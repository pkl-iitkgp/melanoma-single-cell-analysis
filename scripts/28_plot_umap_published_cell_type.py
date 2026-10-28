from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/umap/A_all_2000_HVGs_umap.h5ad"
OUT = ROOT / "figures/umap/A_all_2000_HVGs_umap_published_cell_type_clean.png"

if not INPUT.exists():
    raise FileNotFoundError(f"Missing UMAP object: {INPUT}")

adata = sc.read_h5ad(INPUT)

if "X_umap" not in adata.obsm:
    raise ValueError("Saved UMAP coordinates are missing.")
if "published_cell_type" not in adata.obs:
    raise ValueError("Published cell-type annotations are missing.")

print("Published cell-type counts:")
print(
    adata.obs["published_cell_type"]
    .value_counts(dropna=False)
    .to_string()
)

OUT.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    adata,
    color="published_cell_type",
    title="Melanoma cells by published cell-type annotation",
    frameon=False,
    size=10,
    alpha=0.8,
    legend_loc="right margin",
    show=False,
)

fig = plt.gcf()
fig.set_size_inches(11, 7)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("\nSaved:", OUT)
print("Cells:", adata.n_obs)
print("Existing UMAP coordinates reused; no embedding recomputed.")
