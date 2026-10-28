from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

counts_path = ROOT / "results/clustering/A_all_2000_HVGs_cluster_counts.csv"
annotations_path = ROOT / "results/markers/cluster_annotations_working.csv"
output_path = ROOT / "results/markers/cluster_annotations_with_counts.csv"

counts = pd.read_csv(counts_path)
annotations = pd.read_csv(annotations_path)

print("Count-file columns:", counts.columns.tolist())
print("Annotation-file columns:", annotations.columns.tolist())

# Identify the cluster-ID column in the count file.
possible_cluster_cols = ["cluster", "leiden", "leiden_r05"]
cluster_col = next(
    (col for col in possible_cluster_cols if col in counts.columns),
    None
)

if cluster_col is None:
    raise ValueError(
        "Could not identify the cluster-ID column. "
        f"Available columns: {counts.columns.tolist()}"
    )

counts = counts.rename(columns={cluster_col: "cluster"})
counts["cluster"] = counts["cluster"].astype(str)
annotations["cluster"] = annotations["cluster"].astype(str)

# Merge the original counts with the working biological annotations.
merged = annotations.merge(
    counts,
    on="cluster",
    how="left",
    validate="one_to_one",
    suffixes=("", "_count")
)

if merged["working_annotation"].isna().any():
    raise ValueError("Some clusters are missing working annotations.")

if merged.isna().all(axis=1).any():
    raise ValueError("At least one merged row contains only missing values.")

merged.to_csv(output_path, index=False)

print("\nMerged table saved:", output_path)
print("Rows:", len(merged))
print("Columns:", merged.columns.tolist())
print("\nPreview:")
print(merged.head(6).to_string(index=False))
