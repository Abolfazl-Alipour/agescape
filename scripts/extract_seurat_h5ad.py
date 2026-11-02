#!/usr/bin/env python3

import argparse
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        parsed = rdata.parser.parse_file(args.input)
        slots = rdata.conversion.convert(parsed.object.attributes)
    metadata = slots["meta.data"].copy()
    metadata.columns = [str(column) for column in metadata.columns]
    metadata.index = metadata.index.astype(str)
    counts = slots["assays"]["RNA"].counts
    matrix = sp.csc_matrix(
        (counts.x, counts.i, counts.p), shape=tuple(counts.Dim), dtype=np.float32
    )
    genes = pd.Index(counts.Dimnames[0].astype(str))
    cells = pd.Index(counts.Dimnames[1].astype(str))
    if not cells.equals(metadata.index):
        raise ValueError("Count-matrix cells and Seurat metadata cells do not match")
    result = ad.AnnData(X=matrix.T.tocsr(), obs=metadata, var=pd.DataFrame(index=genes))
    result.var_names_make_unique()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    ad.settings.allow_write_nullable_strings = True
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.unlink(missing_ok=True)
    result.write_h5ad(temporary, compression="gzip")
    temporary.replace(args.output)
    print(f"Wrote {args.output}: {result.shape[0]} cells, {result.shape[1]} genes")
    print(f"Metadata columns: {list(result.obs.columns)}")


if __name__ == "__main__":
    main()
