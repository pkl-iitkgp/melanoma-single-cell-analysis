from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("results/qc")
FIG = Path("figures/qc")

INPUT = OUT / "cell_qc_metrics.csv"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

print("Loading previously calculated QC metrics...")
qc = pd.read_csv(INPUT, index_col=0)

required = [
    "malignant_status",
    "published_cell_type",
    "tumor_sample",
    "n_genes_detected",
]

missing = [column for column in required if column not in qc.columns]
if missing:
    raise ValueError(f"Missing required columns in {INPUT}: {missing}")

print("Cells:", len(qc))

status_summary = qc.groupby(
    "malignant_status", observed=True
).agg(
    n_cells=("n_genes_detected", "size"),
    median_genes=("n_genes_detected", "median"),
    mean_genes=("n_genes_detected", "mean"),
    min_genes=("n_genes_detected", "min"),
    max_genes=("n_genes_detected", "max"),
).round(2)

status_summary.to_csv(OUT / "qc_by_malignancy_status.csv")

type_summary = qc.groupby(
    "published_cell_type", observed=True
).agg(
    n_cells=("n_genes_detected", "size"),
    median_genes=("n_genes_detected", "median"),
    mean_genes=("n_genes_detected", "mean"),
    min_genes=("n_genes_detected", "min"),
    max_genes=("n_genes_detected", "max"),
).round(2).sort_values("n_cells", ascending=False)

type_summary.to_csv(OUT / "qc_by_published_cell_type.csv")

print("\nDetected genes by malignancy status:")
print(status_summary.to_string())

print("\nDetected genes by published cell type:")
print(type_summary.to_string())

def save_boxplot(group_column, order, output_name, xlabel, title):
    values = [
        qc.loc[qc[group_column] == group, "n_genes_detected"].to_numpy()
        for group in order
    ]

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 1.1), 5))
    ax.boxplot(values, tick_labels=order, showfliers=False)
    ax.set_ylabel("Genes with nonzero expression")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(FIG / output_name, dpi=300)
    plt.close(fig)

status_order = [
    x for x in ["malignant", "non-malignant", "unresolved"]
    if x in qc["malignant_status"].unique()
]

type_order = type_summary.index.tolist()

save_boxplot(
    "malignant_status",
    status_order,
    "detected_genes_by_malignancy_status.png",
    "Published malignancy status",
    "Detected genes by malignancy status",
)

save_boxplot(
    "published_cell_type",
    type_order,
    "detected_genes_by_published_cell_type.png",
    "Published cell type",
    "Detected genes by published cell type",
)

print("\nQC comparison complete.")
print("Tables saved in:", OUT)
print("Figures saved in:", FIG)
print("No cells were removed.")
