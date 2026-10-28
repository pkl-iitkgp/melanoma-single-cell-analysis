import numpy as np
import pandas as pd
import scanpy as sc

path = "data/processed/GSE72056_melanoma_pca_input.h5ad"
adata = sc.read_h5ad(path)

X = np.asarray(adata.X)
detection = (X > 0).mean(axis=0)
variance = X.var(axis=0)
mean = X.mean(axis=0)

summary = pd.DataFrame({
    "gene_symbol": adata.var["gene_symbol"].astype(str).to_numpy(),
    "mean_expression": mean,
    "variance": variance,
    "detection_fraction": detection,
})

print("===== FEATURE SPARSITY =====")
for threshold in [0.01, 0.05, 0.10, 0.25, 0.50]:
    count = int((detection >= threshold).sum())
    print(f"Detected in >= {threshold:.0%} of cells: {count} / {adata.n_vars}")

print("\n===== LOWEST-DETECTION FEATURES =====")
print(
    summary.sort_values("detection_fraction")
    .head(15).round(4).to_string(index=False)
)

print("\n===== HIGHEST-VARIANCE FEATURES =====")
print(
    summary.sort_values("variance", ascending=False)
    .head(15).round(4).to_string(index=False)
)

total_variance = variance.sum()
ranked = np.sort(variance)[::-1]
print("\n===== VARIANCE CONCENTRATION =====")
for n in [10, 50, 100, 500]:
    fraction = ranked[:n].sum() / total_variance
    print(f"Top {n} genes account for {fraction:.2%} of summed gene variance")

print("\nDiagnostic complete; no files were modified.")
