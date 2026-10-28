from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = Path(".")
INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
OUT = ROOT / "results/qc"
FIG = ROOT / "figures/qc"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT)
X = np.asarray(adata.X, dtype=np.float32)

# The values are already transformed; calculate variability directly.
gene_mean = X.mean(axis=0)
gene_var = X.var(axis=0, ddof=1)
detection_fraction = (X > 0).mean(axis=0)

results = pd.DataFrame({
    "gene_symbol": adata.var["gene_symbol"].astype(str).values,
    "gene_mean": gene_mean,
    "gene_variance": gene_var,
    "detection_fraction": detection_fraction,
})

results = results.sort_values(
    "gene_variance", ascending=False
).reset_index(drop=True)

results.insert(0, "variance_rank", np.arange(1, len(results) + 1))

results.to_csv(OUT / "gene_variability_ranked.csv", index=False)
results.head(1000).to_csv(OUT / "top_1000_variable_genes.csv", index=False)

# Visualize the variance distribution.
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(results["gene_variance"], bins=100)
ax.set_xlabel("Gene variance on published expression scale")
ax.set_ylabel("Number of genes")
ax.set_title("Gene variability across all melanoma cells")
fig.tight_layout()
fig.savefig(FIG / "gene_variance_distribution.png", dpi=300)
plt.close(fig)

# Variance vs mean helps identify mean-dependent variability.
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(
    results["gene_mean"],
    results["gene_variance"],
    s=5,
    alpha=0.35,
)
ax.set_xlabel("Mean expression")
ax.set_ylabel("Expression variance")
ax.set_title("Mean–variance relationship")
fig.tight_layout()
fig.savefig(FIG / "gene_mean_variance.png", dpi=300)
plt.close(fig)

print("Dimensions:", adata.n_obs, "cells x", adata.n_vars, "genes")
print("Saved:", OUT / "gene_variability_ranked.csv")
print("Saved:", OUT / "top_1000_variable_genes.csv")
print("Saved figures:")
print(" ", FIG / "gene_variance_distribution.png")
print(" ", FIG / "gene_mean_variance.png")

print("\nTop 25 variable genes:")
print(results.head(25).to_string(index=False))

print("\nVariance summary:")
print(results["gene_variance"].describe().to_string())

print("\nGenes detected in at least 10% of cells:",
      int((results["detection_fraction"] >= 0.10).sum()))
