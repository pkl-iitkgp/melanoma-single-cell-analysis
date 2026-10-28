from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

INPUT = Path("data/processed/GSE72056_melanoma.h5ad")
OUT = Path("results/qc")
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT)

print("\n=== DATASET DIMENSIONS ===")
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

X = adata.X
if hasattr(X, "toarray"):
    X = X.toarray()

X = np.asarray(X)
print("\n=== EXPRESSION VALUE DIAGNOSTICS ===")
print("Minimum:", float(np.min(X)))
print("Maximum:", float(np.max(X)))
print("Mean:", float(np.mean(X)))
print("Median:", float(np.median(X)))
print("Fraction exactly zero:", float(np.mean(X == 0)))

quantiles = np.quantile(X, [0, 0.25, 0.5, 0.75, 0.9, 0.99, 1])
print("Quantiles [0, 25, 50, 75, 90, 99, 100%]:")
print(quantiles)

print("\n=== ANNData OBSERVATION METADATA ===")
print(adata.obs.dtypes.to_string())

print("\n=== MALIGNANCY STATUS COUNTS ===")
print(adata.obs["malignant_status"].value_counts(dropna=False).to_string())

print("\n=== PUBLISHED CELL TYPE COUNTS ===")
print(adata.obs["published_cell_type"].value_counts(dropna=False).to_string())

print("\n=== TUMOR SAMPLE COUNTS ===")
print(adata.obs["tumor_sample"].value_counts().sort_index().to_string())

print("\n=== SAMPLE–MALIGNANCY CROSS-TABULATION ===")
cross = pd.crosstab(
    adata.obs["tumor_sample"],
    adata.obs["malignant_status"],
)
print(cross.to_string())
cross.to_csv(OUT / "sample_by_malignancy_status.csv")

print("\n=== GENE SYMBOL EXAMPLES ===")
print(adata.var["gene_symbol"].head(15).to_string())

report = [
    f"Cells: {adata.n_obs}",
    f"Genes: {adata.n_vars}",
    f"Minimum expression: {np.min(X)}",
    f"Maximum expression: {np.max(X)}",
    f"Mean expression: {np.mean(X)}",
    f"Median expression: {np.median(X)}",
    f"Fraction exactly zero: {np.mean(X == 0)}",
    "",
    "Expression quantiles:",
    str(quantiles),
    "",
    "Observation columns:",
    str(list(adata.obs.columns)),
    "",
    "Original gene-symbol column exists: "
    + str("gene_symbol" in adata.var.columns),
]

(OUT / "expression_metadata_diagnostics.txt").write_text(
    "\n".join(report) + "\n"
)

print("\nDiagnostic report saved to:")
print(OUT / "expression_metadata_diagnostics.txt")
print("Sample cross-tabulation saved to:")
print(OUT / "sample_by_malignancy_status.csv")
