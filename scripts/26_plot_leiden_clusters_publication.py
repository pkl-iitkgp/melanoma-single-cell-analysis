from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/clustering/A_all_2000_HVGs_leiden_r05.h5ad"
OUT = ROOT / "figures/clustering/A_all_2000_HVGs_umap_leiden_r05_clean.png"

if not INPUT.exists():
    raise FileNotFoundError(f"Missing clustered AnnData: {INPUT}")

adata = sc.read_h5ad(INPUT)

if "X_umap" not in adata.obsm:
    raise ValueError("Saved UMAP coordinates are missing.")
if "leiden_r05" not in adata.obs:
    raise ValueError("Leiden cluster labels are missing.")

OUT.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    adata,
    color="leiden_r05",
    title="Melanoma single-cell clusters (Leiden resolution 0.5)",
    legend_loc="on data",
    legend_fontsize=9,
    legend_fontweight="bold",
    frameon=False,
    size=9,
    alpha=0.85,
    show=False,
)

fig = plt.gcf()
fig.set_size_inches(9, 7)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Saved:", OUT)
print("Cells:", adata.n_obs)
print("Clusters:", adata.obs["leiden_r05"].nunique())
print("Existing UMAP coordinates reused; no reclustering performed.")
