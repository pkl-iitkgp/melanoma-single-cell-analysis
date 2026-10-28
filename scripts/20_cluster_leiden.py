from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

INPUTS = {
    "A_all_2000_HVGs":
        "results/umap/A_all_2000_HVGs_umap.h5ad",
    "B_HVGs_detected_in_at_least_5pct":
        "results/umap/B_HVGs_detected_in_at_least_5pct_umap.h5ad",
}

OUT = Path("results/clustering")
FIG = Path("figures/clustering")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

summary = [
    "LEIDEN CLUSTERING COMPARISON",
    "============================",
    "Input: saved UMAP objects with existing neighbor graphs",
    "Method: Leiden",
    "Resolution: 0.5",
    "Random seed: 0",
    "",
]

for name, input_path in INPUTS.items():
    print(f"\n===== Processing {name} =====")
    adata = sc.read_h5ad(input_path)

    if "connectivities" not in adata.obsp:
        raise ValueError(f"Neighbor connectivities missing: {input_path}")

    # Store cluster labels as a separate annotation.
    sc.tl.leiden(
        adata,
        resolution=0.5,
        random_state=0,
        key_added="leiden_r05",
        flavor="leidenalg",
        directed=False,
    )

    counts = adata.obs["leiden_r05"].value_counts().sort_index()
    counts_df = counts.rename_axis("cluster").reset_index(name="n_cells")
    counts_df["fraction_of_cells"] = counts_df["n_cells"] / adata.n_obs
    counts_df.to_csv(OUT / f"{name}_cluster_counts.csv", index=False)

    # Save the cluster annotation on the UMAP.
    sc.pl.umap(
        adata,
        color="leiden_r05",
        show=False,
        title=f"{name}: Leiden resolution 0.5",
        legend_loc="on data",
    )
    fig = plt.gcf()
    fig.savefig(
        FIG / f"{name}_umap_leiden_r05.png",
        dpi=160,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Compare clusters with existing metadata annotations.
    for column in ["malignant_status", "published_cell_type", "tumor_sample"]:
        table = pd.crosstab(
            adata.obs["leiden_r05"],
            adata.obs[column],
            dropna=False,
        )
        table.to_csv(OUT / f"{name}_cluster_by_{column}.csv")

    adata.write_h5ad(OUT / f"{name}_leiden_r05.h5ad")

    summary.extend([
        f"Feature set: {name}",
        f"Cells: {adata.n_obs}",
        f"Genes: {adata.n_vars}",
        f"Number of Leiden clusters: {len(counts)}",
        "Cluster sizes:",
        counts.to_string(),
        "",
    ])

    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Number of clusters: {len(counts)}")
    print("Cluster sizes:")
    print(counts.to_string())
    print("Saved cluster counts, metadata cross-tabs, plot, and AnnData.")

(OUT / "leiden_comparison_summary.txt").write_text(
    "\n".join(summary) + "\n"
)

print("\nLeiden clustering comparison completed.")
print(f"Summary: {OUT / 'leiden_comparison_summary.txt'}")
