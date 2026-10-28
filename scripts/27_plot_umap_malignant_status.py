from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/umap/A_all_2000_HVGs_umap.h5ad"
OUT = ROOT / "figures/umap/A_all_2000_HVGs_umap_malignant_status_clean.png"

if not INPUT.exists():
    raise FileNotFoundError(f"Missing UMAP object: {INPUT}")

adata = sc.read_h5ad(INPUT)

required = ["X_umap"]
if "X_umap" not in adata.obsm:
    raise ValueError("Saved UMAP coordinates are missing.")
if "malignant_status" not in adata.obs:
    raise ValueError("malignant_status metadata is missing.")

print("Malignant-status categories:")
print(adata.obs["malignant_status"].value_counts(dropna=False).to_string())

OUT.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    adata,
    color="malignant_status",
    title="Melanoma cells by malignant status",
    frameon=False,
    size=10,
    alpha=0.8,
    show=False,
)

fig = plt.gcf()
fig.set_size_inches(9, 7)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("\nSaved:", OUT)
print("Cells:", adata.n_obs)
print("Existing UMAP coordinates reused; no embedding recomputed.")
