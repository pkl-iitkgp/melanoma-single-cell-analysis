from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
CLUSTERED = ROOT / "results/clustering/A_all_2000_HVGs_leiden_r05.h5ad"

OUT_DIR = ROOT / "results/markers"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_CELLS_DETECTED = 3
TOP_N = 25

print("Loading original expression data...")
adata = sc.read_h5ad(INPUT)
clustered = sc.read_h5ad(CLUSTERED)

# Transfer cluster labels by cell ID, not by row position.
labels = clustered.obs["leiden_r05"].reindex(adata.obs_names)

if labels.isna().any():
    raise ValueError("Some cells have no matching Leiden cluster label.")

adata.obs["leiden_r05"] = labels.astype("category")

# Confirm gene symbols are available.
if "gene_symbol" not in adata.var.columns:
    raise ValueError("Expected adata.var['gene_symbol'] was not found.")

symbols = adata.var["gene_symbol"].astype(str)
valid_symbols = symbols.notna() & symbols.ne("") & (symbols != "nan")

# Filter genes detected in at least three cells.
detected_per_gene = np.asarray((adata.X > 0).sum(axis=0)).ravel()
keep = (detected_per_gene >= MIN_CELLS_DETECTED) & valid_symbols.to_numpy()

n_genes_before = adata.n_vars
adata = adata[:, keep].copy()

# Use readable gene symbols for reporting. Retain unique names when symbols
# are duplicated, so distinct source rows are not silently merged.
adata.var["source_gene_id"] = adata.var_names.astype(str)
adata.var_names = pd.Index(adata.var["gene_symbol"].astype(str))
adata.var_names_make_unique(join="-dup")

print(f"Cells: {adata.n_obs}")
print(f"Genes before filtering: {n_genes_before}")
print(f"Genes retained: {adata.n_vars}")
print(f"Genes removed: {n_genes_before - adata.n_vars}")
print(f"Clusters: {adata.obs['leiden_r05'].nunique()}")
print("Running Wilcoxon marker analysis...")

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always", RuntimeWarning)
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden_r05",
        method="wilcoxon",
        use_raw=False,
        pts=True,
        tie_correct=True,
        corr_method="benjamini-hochberg",
    )

runtime_warnings = [
    str(w.message) for w in caught
    if issubclass(w.category, RuntimeWarning)
]

markers = sc.get.rank_genes_groups_df(adata, group=None)

# Remove rows with invalid numerical statistics before ranking.
numeric_cols = [
    c for c in ["scores", "logfoldchanges", "pvals", "pvals_adj"]
    if c in markers.columns
]
if numeric_cols:
    finite = np.isfinite(markers[numeric_cols].to_numpy(dtype=float)).all(axis=1)
    n_invalid = int((~finite).sum())
    markers = markers.loc[finite].copy()
else:
    n_invalid = 0

# P-values may underflow to zero for very strong results; preserve them.
# Sort within each cluster by adjusted p-value and then by test score.
markers = markers.sort_values(
    ["group", "pvals_adj", "scores"],
    ascending=[True, True, False],
    kind="stable",
)

markers.to_csv(OUT_DIR / "leiden_r05_all_cluster_markers.csv", index=False)

top_markers = (
    markers.groupby("group", observed=True, sort=True)
    .head(TOP_N)
    .copy()
)
top_markers.to_csv(
    OUT_DIR / f"leiden_r05_top{TOP_N}_markers_per_cluster.csv",
    index=False,
)

summary_path = OUT_DIR / "leiden_r05_marker_analysis_summary.txt"
with summary_path.open("w") as f:
    f.write("Leiden cluster marker analysis\n")
    f.write("==============================\n")
    f.write(f"Input: {INPUT.relative_to(ROOT)}\n")
    f.write(f"Cluster labels: {CLUSTERED.relative_to(ROOT)}\n")
    f.write("Test: Scanpy Wilcoxon rank-sum test\n")
    f.write("Multiple testing correction: Benjamini-Hochberg\n")
    f.write(f"Minimum detected cells per gene: {MIN_CELLS_DETECTED}\n")
    f.write(f"Cells analysed: {adata.n_obs}\n")
    f.write(f"Genes before filtering: {n_genes_before}\n")
    f.write(f"Genes after filtering: {adata.n_vars}\n")
    f.write(f"Invalid-statistic rows excluded: {n_invalid}\n")
    f.write(f"Marker rows retained: {len(markers)}\n")
    f.write(f"Clusters: {adata.obs['leiden_r05'].nunique()}\n")
    f.write(f"Top markers per cluster exported: {TOP_N}\n")
    f.write(f"Runtime warnings recorded: {len(runtime_warnings)}\n")
    for message in sorted(set(runtime_warnings)):
        f.write(f"Warning: {message}\n")
    f.write("\nInterpretation notes:\n")
    f.write("- Markers are exploratory cluster-associated genes.\n")
    f.write("- Cluster marker results do not by themselves prove cell identity.\n")
    f.write("- The original expression object was not modified on disk.\n")

print("\nMarker analysis completed.")
print("All markers:", OUT_DIR / "leiden_r05_all_cluster_markers.csv")
print("Top markers:", OUT_DIR / f"leiden_r05_top{TOP_N}_markers_per_cluster.csv")
print("Summary:", summary_path)
print("Rows retained:", len(markers))
print("Rows excluded for invalid statistics:", n_invalid)
print("Runtime warnings recorded:", len(runtime_warnings))

print("\nPreview of top markers:")
print(
    top_markers[
        ["group", "names", "scores", "logfoldchanges", "pvals_adj"]
    ].head(15).to_string(index=False)
)
