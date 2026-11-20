# Block 4: Reverse signature transport

## Status

INCONCLUSIVE: blocked before analysis because a required input is missing.

## Missing required input

Block 4 requires 16 Tabula Muris Senis mouse target files matching:

`data/interim/m1-qc-cleaned-*.h5ad`

Direct inspection of `data/interim/` found zero matching files. The expected per-tissue mouse targets are therefore unavailable at the stated paths.

The shared rules require stopping the block when a stated input is missing and prohibit substituting another input. Human signature derivation, ortholog mapping, TMS scoring, matched random-set tests, and transport verdict calculations were not run. No result CSVs were generated.

## Required resolution

Restore the 16 per-tissue TMS files at the stated `data/interim/m1-qc-cleaned-*.h5ad` paths, then rerun Block 4 exactly as specified. No other block was started.
