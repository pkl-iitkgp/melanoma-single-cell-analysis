from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

markers_path = ROOT / "results/markers/leiden_r05_all_cluster_markers.csv"
validation_path = ROOT / "results/markers/cluster_annotation_validation.csv"
clustered_path = ROOT / "results/clustering/A_all_2000_HVGs_leiden_r05.h5ad"

markers = pd.read_csv(markers_path)
validation = pd.read_csv(validation_path)

# Review clear lineage candidates and clusters requiring extra scrutiny.
review_clusters = [3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19]

print("TOP 15 MARKERS FOR SELECTED CLUSTERS")
print("=" * 60)

for cluster in review_clusters:
    subset = markers[markers["group"].astype(str) == str(cluster)].copy()
    subset = subset.sort_values(["pvals_adj", "scores"], ascending=[True, False])
    genes = subset["names"].head(15).astype(str).tolist()
    print(f"\nCluster {cluster}: {', '.join(genes)}")

print("\n\nPUBLISHED ANNOTATION AND MALIGNANT STATUS")
print("=" * 60)

cols = ["cluster", "n_cells"]
for col in validation.columns:
    if col.startswith("malignant_status__") or col.startswith("published_cell_type__"):
        cols.append(col)

print(validation.loc[validation["cluster"].isin(review_clusters), cols]
        .fillna(0).round(2).to_string(index=False))
