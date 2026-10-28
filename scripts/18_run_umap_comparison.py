from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

INPUTS = {
    "A_all_2000_HVGs":
        "results/pca/A_all_2000_HVGs_pca.h5ad",
    "B_HVGs_detected_in_at_least_5pct":
        "results/pca/B_HVGs_detected_in_at_least_5pct_pca.h5ad",
}

OUT = Path("results/umap")
FIG = Path("figures/umap")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

summary = [
    "UMAP COMPARISON",
    "===============",
    "Input: saved PCA objects",
    "Neighbors: n_neighbors=15, n_pcs=30",
    "UMAP: min_dist=0.3, random_state=0",
    "No additional normalization or log transformation",
    "",
]

for name, input_path in INPUTS.items():
    print(f"\n===== Processing {name} =====")

    adata = sc.read_h5ad(input_path)

    if "X_pca" not in adata.obsm:
        raise ValueError(f"PCA coordinates missing in {input_path}")

    if adata.obsm["X_pca"].shape[1] < 30:
        raise ValueError(f"Fewer than 30 PCs available in {input_path}")

    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=30,
        random_state=0,
    )

    sc.tl.umap(
        adata,
        min_dist=0.3,
        random_state=0,
    )

    coords = np.asarray(adata.obsm["X_umap"])
    if not np.isfinite(coords).all():
        raise ValueError("UMAP coordinates contain non-finite values")

    for column in [
        "malignant_status",
        "published_cell_type",
        "tumor_sample",
    ]:
        if column not in adata.obs:
            raise ValueError(f"Missing metadata column: {column}")

        sc.pl.umap(
            adata,
            color=column,
            show=False,
            title=f"{name}: {column}",
        )
        fig = plt.gcf()
        fig.savefig(
            FIG / f"{name}_umap_{column}.png",
            dpi=160,
            bbox_inches="tight",
        )
        plt.close(fig)

    adata.write_h5ad(OUT / f"{name}_umap.h5ad")

    summary.extend([
        f"Feature set: {name}",
        f"Cells: {adata.n_obs}",
        f"Genes: {adata.n_vars}",
        f"UMAP coordinates finite: {np.isfinite(coords).all()}",
        f"Saved object: {OUT / f'{name}_umap.h5ad'}",
        "",
    ])

    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"UMAP coordinates shape: {coords.shape}")
    print("Saved three UMAP plots and the UMAP AnnData object.")

(OUT / "umap_comparison_summary.txt").write_text(
    "\n".join(summary) + "\n"
)

print("\nUMAP comparison completed.")
print(f"Summary: {OUT / 'umap_comparison_summary.txt'}")
