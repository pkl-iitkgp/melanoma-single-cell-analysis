from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

INPUT = Path("data/processed/GSE72056_melanoma.h5ad")
OUT = Path("results/qc")
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT)

# The input contains already-transformed expression values.
# Do not apply normalize_total() or log1p() here.
X = adata.X

if hasattr(X, "toarray"):
    X = X.toarray()

X = np.asarray(X, dtype=np.float32)

print("Calculating gene-level variability...")
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

gene_mean = X.mean(axis=0)
gene_variance = X.var(axis=0, ddof=1)
gene_detection = (X > 0).mean(axis=0)

summary = pd.DataFrame({
    "gene_symbol": adata.var["gene_symbol"].astype(str).to_numpy(),
    "mean_expression": gene_mean,
    "variance": gene_variance,
    "fraction_cells_detected": gene_detection,
})

summary = summary.sort_values("variance", ascending=False)

summary.to_csv(
    OUT / "gene_variability_all.csv",
    index=False
)

summary.head(100).to_csv(
    OUT / "top_100_variable_genes.csv",
    index=False
)

print("\nTop 20 genes by variance on the published expression scale:")
print(summary.head(20).to_string(index=False))

print("\nSaved:")
print(OUT / "gene_variability_all.csv")
print(OUT / "top_100_variable_genes.csv")

print("\nNo genes or cells were removed.")
print("The original AnnData file was not modified.")
