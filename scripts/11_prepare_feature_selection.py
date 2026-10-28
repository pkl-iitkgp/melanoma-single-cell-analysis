from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path.home() / "single_cell_project"
INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
OUTPUT = ROOT / "results/qc/feature_selection_precheck.txt"

print("Scanpy version:", version("scanpy"))
print("Reading:", INPUT)

adata = sc.read_h5ad(INPUT)

print("\nObject summary")
print(adata)
print("X type:", type(adata.X))
print("X dtype:", adata.X.dtype)
print("Shape:", adata.shape)

X = adata.X

if hasattr(X, "toarray"):
    # The matrix is expected to be dense for this dataset.
    # Avoid converting a sparse matrix to a dense matrix unnecessarily.
    print("Sparse matrix:", True)
    values_min = float(X.min())
    values_max = float(X.max())
else:
    print("Sparse matrix:", False)
    values_min = float(np.min(X))
    values_max = float(np.max(X))

print("Minimum expression:", values_min)
print("Maximum expression:", values_max)
print("Contains non-finite values:", not np.isfinite(X).all())

print("\nMetadata columns:")
print(adata.obs.columns.tolist())

print("\nTumor samples:", adata.obs["tumor_sample"].nunique())
print("Malignancy categories:")
print(adata.obs["malignant_status"].value_counts(dropna=False).to_string())

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w") as f:
    f.write(f"Scanpy version: {version('scanpy')}\n")
    f.write(f"Shape: {adata.shape}\n")
    f.write(f"X dtype: {adata.X.dtype}\n")
    f.write(f"Minimum expression: {values_min}\n")
    f.write(f"Maximum expression: {values_max}\n")
    f.write(
        "Contains non-finite values: "
        f"{not np.isfinite(adata.X).all()}\n"
    )
    f.write(f"Tumor samples: {adata.obs['tumor_sample'].nunique()}\n")
    f.write("\nMalignancy categories:\n")
    f.write(adata.obs["malignant_status"].value_counts(
        dropna=False
    ).to_string())
    f.write("\n")

print("\nSaved diagnostic report:", OUTPUT)
