from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

INPUT = Path("data/processed/GSE72056_melanoma.h5ad")
OUT = Path("results/qc")
FIG = Path("figures/qc")

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

print("Step 1: Loading the validated dataset...")
adata = sc.read_h5ad(INPUT)
X = adata.X

print("Dimensions:", adata.n_obs, "cells x", adata.n_vars, "gene rows")

print("\nStep 2: Calculating expression-based QC metrics...")

# The matrix contains processed log-transformed expression values.
# Therefore, these are NOT raw sequencing library-size metrics.
adata.obs["n_genes_detected"] = np.count_nonzero(X > 0, axis=1)
adata.obs["transformed_expression_sum"] = np.asarray(X.sum(axis=1)).ravel()
adata.obs["mean_transformed_expression"] = np.asarray(X.mean(axis=1)).ravel()

qc = adata.obs[
    [
        "tumor_sample",
        "malignant_status",
        "published_cell_type",
        "n_genes_detected",
        "transformed_expression_sum",
        "mean_transformed_expression",
    ]
].copy()

qc.to_csv(OUT / "cell_qc_metrics.csv")

sample_summary = qc.groupby("tumor_sample", observed=True).agg(
    n_cells=("n_genes_detected", "size"),
    median_detected_genes=("n_genes_detected", "median"),
    mean_detected_genes=("n_genes_detected", "mean"),
    min_detected_genes=("n_genes_detected", "min"),
    max_detected_genes=("n_genes_detected", "max"),
)

sample_summary = sample_summary.sort_values("n_cells", ascending=False)
sample_summary.to_csv(OUT / "sample_qc_summary.csv")

print("\nStep 3: Overall detected-gene summary")
print(qc["n_genes_detected"].describe().round(2).to_string())

print("\nStep 4: Cell counts and median detected genes by sample")
print(sample_summary.round(2).to_string())

print("\nStep 5: Creating QC figures...")

# Figure 1: Distribution across all cells.
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(qc["n_genes_detected"], bins=40, edgecolor="black")
ax.set_xlabel("Number of genes with nonzero expression")
ax.set_ylabel("Number of cells")
ax.set_title("GSE72056: detected genes per cell")
fig.tight_layout()
fig.savefig(FIG / "detected_genes_distribution.png", dpi=300)
plt.close(fig)

# Figure 2: Detected genes across tumor sample codes.
plot_data = [
    group["n_genes_detected"].to_numpy()
    for _, group in qc.groupby("tumor_sample", sort=True)
]
sample_labels = [
    str(sample)
    for sample, _ in qc.groupby("tumor_sample", sort=True)
]

fig, ax = plt.subplots(figsize=(12, 6))
ax.boxplot(plot_data, showfliers=False)
ax.set_xticks(range(1, len(sample_labels) + 1))
ax.set_xticklabels(sample_labels, rotation=45, ha="right")
ax.set_xlabel("Tumor sample code")
ax.set_ylabel("Genes with nonzero expression")
ax.set_title("GSE72056: detected genes by tumor sample")
fig.tight_layout()
fig.savefig(FIG / "detected_genes_by_sample.png", dpi=300)
plt.close(fig)

print("\nQC assessment complete.")
print("Cell metrics:", OUT / "cell_qc_metrics.csv")
print("Sample summary:", OUT / "sample_qc_summary.csv")
print("Figures saved in:", FIG)
print("\nNo cells were filtered in this step.")
