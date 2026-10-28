from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd

OUT = Path("results/pca")
FIG = Path("figures/pca")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

feature_sets = {
    "A_all_2000_HVGs":
        "data/processed/GSE72056_melanoma_A_all_2000_HVGs.h5ad",
    "B_HVGs_detected_in_at_least_5pct":
        "data/processed/GSE72056_melanoma_B_HVGs_detected_in_at_least_5pct.h5ad",
}

summary_lines = [
    "PCA COMPARISON",
    "==============",
    "Input values: processed GEO expression values",
    "Transformation: no additional log transform",
    "Scaling: gene-wise centering/scaling; max_value=10",
    "Principal components requested: 30",
    "",
]

for name, path in feature_sets.items():
    print(f"\n===== Processing {name} =====")
    adata = sc.read_h5ad(path)

    # Work only on this in-memory copy; source files remain unchanged.
    sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.tl.pca(
        adata,
        n_comps=30,
        zero_center=True,
        svd_solver="arpack",
        random_state=0,
    )

    ratios = np.asarray(adata.uns["pca"]["variance_ratio"])
    variances = np.asarray(adata.uns["pca"]["variance"])

    pca_table = pd.DataFrame({
        "PC": [f"PC{i}" for i in range(1, len(ratios) + 1)],
        "variance": variances,
        "variance_ratio": ratios,
        "cumulative_variance_ratio": np.cumsum(ratios),
    })
    pca_table.to_csv(OUT / f"{name}_pca_variance.csv", index=False)

    sc.pl.pca(
        adata,
        color="malignant_status",
        show=False,
        save=None,
    )
    fig = plt.gcf()
    fig.savefig(FIG / f"{name}_pca_malignant_status.png",
                dpi=160, bbox_inches="tight")
    plt.close(fig)

    sc.pl.pca(
        adata,
        color="published_cell_type",
        show=False,
    )
    fig = plt.gcf()
    fig.savefig(FIG / f"{name}_pca_published_cell_type.png",
                dpi=160, bbox_inches="tight")
    plt.close(fig)

    sc.pl.pca(
        adata,
        color="tumor_sample",
        show=False,
    )
    fig = plt.gcf()
    fig.savefig(FIG / f"{name}_pca_tumor_sample.png",
                dpi=160, bbox_inches="tight")
    plt.close(fig)

    adata.write_h5ad(OUT / f"{name}_pca.h5ad")

    summary_lines.extend([
        f"Feature set: {name}",
        f"Cells: {adata.n_obs}",
        f"Genes: {adata.n_vars}",
        f"PC1 variance explained: {ratios[0]:.4%}",
        f"PC2 variance explained: {ratios[1]:.4%}",
        f"First 10 PCs cumulative variance: {ratios[:10].sum():.4%}",
        f"First 30 PCs cumulative variance: {ratios.sum():.4%}",
        f"PCA object: {OUT / f'{name}_pca.h5ad'}",
        "",
    ])

    print(f"Genes: {adata.n_vars}")
    print(f"PC1 variance explained: {ratios[0]:.4%}")
    print(f"PC2 variance explained: {ratios[1]:.4%}")
    print(f"First 10 PCs cumulative variance: {ratios[:10].sum():.4%}")
    print(f"First 30 PCs cumulative variance: {ratios.sum():.4%}")

(OUT / "pca_comparison_summary.txt").write_text(
    "\n".join(summary_lines) + "\n"
)
print("\nSaved PCA objects, variance tables, and plots.")
