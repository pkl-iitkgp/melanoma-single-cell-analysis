from pathlib import Path
import pandas as pd
import scanpy as sc

root = Path(".")
input_file = root / "results/clustering/A_all_2000_HVGs_leiden_r05.h5ad"
output_file = root / "results/annotation/cluster_sample_composition.csv"

adata = sc.read_h5ad(input_file)

required = {"leiden_r05", "tumor_sample"}
missing = required - set(adata.obs.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

rows = []
for cluster, cells in adata.obs.groupby("leiden_r05", observed=True):
    counts = cells["tumor_sample"].value_counts()
    rows.append({
        "cluster": str(cluster),
        "n_cells": len(cells),
        "n_tumor_samples": cells["tumor_sample"].nunique(),
        "dominant_tumor_sample": str(counts.index[0]),
        "dominant_sample_cells": int(counts.iloc[0]),
        "dominant_sample_fraction": counts.iloc[0] / len(cells),
    })

result = pd.DataFrame(rows)
result["cluster_num"] = pd.to_numeric(result["cluster"])
result = result.sort_values("cluster_num").drop(columns="cluster_num")
result.to_csv(output_file, index=False)

print(result.to_string(index=False, formatters={
    "dominant_sample_fraction": "{:.1%}".format
}))
print(f"\nSaved: {output_file}")
