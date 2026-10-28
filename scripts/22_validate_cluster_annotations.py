from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]

EXPR = ROOT / "data/processed/GSE72056_melanoma.h5ad"
CLUSTERS = ROOT / "results/clustering/A_all_2000_HVGs_leiden_r05.h5ad"
MARKERS = ROOT / "results/markers/leiden_r05_all_cluster_markers.csv"
OUT = ROOT / "results/markers/cluster_annotation_validation.csv"

adata = sc.read_h5ad(EXPR)
clustered = sc.read_h5ad(CLUSTERS)

labels = clustered.obs["leiden_r05"].reindex(adata.obs_names)
if labels.isna().any():
    raise ValueError("Some cells could not be matched to cluster labels.")

adata.obs["leiden_r05"] = labels.astype(str)

# Canonical marker panels used as supporting evidence, not rigid rules.
panels = {
    "T_cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "TCF7"],
    "NK_cytotoxic": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "NCR1"],
    "B_cell": ["MS4A1", "CD79A", "CD19", "CD22", "PAX5"],
    "Plasma_cell": ["MZB1", "SDC1", "JCHAIN", "IGJ", "TNFRSF17"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "LYZ", "CSF1R", "CD14", "CD163"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "TCF4", "GZMB"],
    "Mast_cell": ["KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Fibroblast_CAF": ["COL1A1", "COL1A2", "DCN", "LUM", "DPT", "COL3A1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "PLVAP"],
    "Melanocytic_melanoma": ["MLANA", "PMEL", "TYR", "TYRP1", "DCT", "MITF"],
    "Cell_cycle": ["MKI67", "TOP2A", "UBE2C", "CDK1", "CCNB1"],
}

# Use gene_symbol mapping because the original var_names are numeric IDs.
gene_symbols = adata.var["gene_symbol"].astype(str)
symbol_to_positions = {}
for pos, symbol in enumerate(gene_symbols):
    if symbol and symbol != "nan":
        symbol_to_positions.setdefault(symbol, []).append(pos)

X = adata.X
rows = []

for cluster in sorted(adata.obs["leiden_r05"].unique(), key=int):
    cell_mask = (adata.obs["leiden_r05"].to_numpy() == cluster)
    cluster_n = int(cell_mask.sum())

    row = {
        "cluster": int(cluster),
        "n_cells": cluster_n,
    }

    # Percentages based on existing published labels.
    for field in ["malignant_status", "published_cell_type"]:
        if field not in adata.obs.columns:
            continue
        counts = adata.obs.loc[cell_mask, field].value_counts()
        for label, count in counts.items():
            safe_label = str(label).replace(" ", "_").replace("/", "_")
            row[f"{field}__{safe_label}_pct"] = round(100 * count / cluster_n, 2)

    # Mean expression and fraction detected for canonical markers.
    for panel, genes in panels.items():
        panel_means = []
        panel_detection = []

        for gene in genes:
            positions = symbol_to_positions.get(gene, [])
            if not positions:
                continue

            values = np.asarray(X[cell_mask, positions]).reshape(cluster_n, -1)
            panel_means.extend(values.mean(axis=0).tolist())
            panel_detection.extend((values > 0).mean(axis=0).tolist())

        row[f"{panel}__mean_expression"] = (
            round(float(np.mean(panel_means)), 4) if panel_means else np.nan
        )
        row[f"{panel}__mean_fraction_detected"] = (
            round(float(np.mean(panel_detection)), 4)
            if panel_detection else np.nan
        )

    rows.append(row)

result = pd.DataFrame(rows)
result.to_csv(OUT, index=False)

print("Saved:", OUT)
print("Clusters:", len(result))
print("\nCanonical-marker panel summary:")
summary_cols = [
    "cluster", "n_cells",
    "T_cell__mean_fraction_detected",
    "NK_cytotoxic__mean_fraction_detected",
    "B_cell__mean_fraction_detected",
    "Plasma_cell__mean_fraction_detected",
    "Myeloid__mean_fraction_detected",
    "pDC__mean_fraction_detected",
    "Mast_cell__mean_fraction_detected",
    "Fibroblast_CAF__mean_fraction_detected",
    "Endothelial__mean_fraction_detected",
    "Melanocytic_melanoma__mean_fraction_detected",
    "Cell_cycle__mean_fraction_detected",
]
print(result[summary_cols].round(3).to_string(index=False))
