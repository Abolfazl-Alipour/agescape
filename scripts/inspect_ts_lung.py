import anndata as ad
import pandas as pd

print("Opening TS_Lung.h5ad...")
a = ad.read_h5ad("data/raw/TS_Lung.h5ad", backed="r")
print("Shape:", a.shape)
print("\nDonors:")
print(a.obs.groupby(["donor", "gender", "method"], observed=True).size())
print("\nTop cell ontology classes:")
print(a.obs["cell_ontology_class"].value_counts().head(20))
print("\nVar names sample:", a.var_names[:10].tolist())
print("Var columns:", a.var.columns.tolist())

