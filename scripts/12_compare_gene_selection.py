from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path.home() / "single_cell_project"
INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
RANKED = ROOT / "results/qc/gene_variability_ranked.csv"
OUTDIR = ROOT / "results/qc"
FIGDIR = ROOT / "figures/qc"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

print("Loading expression object...")
adata = sc.read_h5ad(INPUT)
X = adata.X

# Calculate gene-wise mean and variance without changing X.
gene_mean = np.asarray(X.mean(axis=0)).ravel()
gene_var = np.asarray(X.var(axis=0, ddof=1)).ravel()
detection_fraction = np.asarray((X > 0).mean(axis=0)).ravel()

stats = pd.DataFrame({
    "gene_id": adata.var_names.astype(str),
    "gene_symbol": adata.var["gene_symbol"].astype(str).to_numpy(),
    "mean_expression": gene_mean,
    "variance": gene_var,
    "detection_fraction": detection_fraction,
})

# Match the original variance ranking.
stats = stats.sort_values(
    "variance", ascending=False
).reset_index(drop=True)
stats["variance_rank"] = np.arange(1, len(stats) + 1)

# Exploratory mean-adjusted variability:
# fit a smooth trend to log(variance) versus log(mean) for genes
# with positive means and variances. Residuals indicate variance
# above/below the fitted trend. This is a ranking, not a formal
# statistical test or a count-based HVG procedure.
valid = (
    (stats["mean_expression"] > 0)
    & (stats["variance"] > 0)
    & np.isfinite(stats["mean_expression"])
    & np.isfinite(stats["variance"])
)

x = np.log1p(stats.loc[valid, "mean_expression"].to_numpy())
y = np.log1p(stats.loc[valid, "variance"].to_numpy())

# Fit a polynomial trend to avoid requiring additional packages.
# Degree 3 captures a broad nonlinear relationship.
coefficients = np.polyfit(x, y, deg=3)
predicted = np.polyval(coefficients, x)
residuals = y - predicted

stats["mean_adjusted_variability"] = np.nan
stats.loc[valid, "mean_adjusted_variability"] = residuals

# Genes with zero variance or zero mean are not informative for
# this residual ranking.
adjusted = stats.loc[valid].sort_values(
    "mean_adjusted_variability", ascending=False
).copy()

top_variance = stats.head(1000).copy()
top_adjusted = adjusted.head(1000).copy()

top_variance.to_csv(
    OUTDIR / "candidate_top_1000_by_variance.csv", index=False
)
top_adjusted.to_csv(
    OUTDIR / "candidate_top_1000_mean_adjusted.csv", index=False
)
stats.to_csv(
    OUTDIR / "gene_selection_comparison_all.csv", index=False
)

variance_ids = set(top_variance["gene_id"])
adjusted_ids = set(top_adjusted["gene_id"])
overlap = variance_ids & adjusted_ids

summary = [
    f"Cells: {adata.n_obs}",
    f"Genes: {adata.n_vars}",
    f"Top variance candidates: {len(variance_ids)}",
    f"Top mean-adjusted candidates: {len(adjusted_ids)}",
    f"Overlap: {len(overlap)}",
    f"Jaccard similarity: {len(overlap) / len(variance_ids | adjusted_ids):.4f}",
    "",
    "Top 20 by variance:",
    top_variance[["gene_symbol", "mean_expression", "variance"]]
        .head(20).to_string(index=False),
    "",
    "Top 20 by mean-adjusted variability:",
    top_adjusted[[
        "gene_symbol", "mean_expression", "variance",
        "mean_adjusted_variability"
    ]].head(20).to_string(index=False),
]

report = "\n".join(summary)
print(report)
(OUTDIR / "gene_selection_comparison_summary.txt").write_text(report)

# Visualize the mean-variance relationship and highlight the
# mean-adjusted candidates.
selected_ids = set(top_adjusted["gene_id"])
selected = stats["gene_id"].isin(selected_ids)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(
    stats["mean_expression"], stats["variance"],
    s=5, alpha=0.20, label="All genes"
)
ax.scatter(
    stats.loc[selected, "mean_expression"],
    stats.loc[selected, "variance"],
    s=8, alpha=0.65, label="Top mean-adjusted candidates"
)
ax.set_xlabel("Mean expression")
ax.set_ylabel("Variance")
ax.set_title("Exploratory mean-adjusted gene variability")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(
    FIGDIR / "mean_adjusted_gene_candidates.png",
    dpi=200
)
plt.close(fig)

print("\nSaved comparison tables and figure.")
print("Original expression object was not modified.")
