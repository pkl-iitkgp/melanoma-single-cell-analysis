from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

INPUT = Path("data/processed/GSE72056_melanoma_pca_input.h5ad")
OUTDIR = Path("data/processed")
REPORT = Path("results/pca/feature_set_comparison.txt")

adata = sc.read_h5ad(INPUT)
X = np.asarray(adata.X)

detection_fraction = (X > 0).mean(axis=0)
variance = X.var(axis=0)

# Set A: all selected HVGs.
mask_all = np.ones(adata.n_vars, dtype=bool)

# Set B: HVGs detected in at least 5% of all cells.
mask_filtered = detection_fraction >= 0.05

sets = {
    "A_all_2000_HVGs": mask_all,
    "B_HVGs_detected_in_at_least_5pct": mask_filtered,
}

lines = [
    "PCA FEATURE SET COMPARISON",
    "==========================",
    f"Input: {INPUT}",
    f"Cells: {adata.n_obs}",
    f"Original selected HVGs: {adata.n_vars}",
    "",
]

for name, mask in sets.items():
    subset = adata[:, mask].copy()
    output = OUTDIR / f"GSE72056_melanoma_{name}.h5ad"

    subset.uns["pca_feature_set"] = {
        "name": name,
        "source_file": str(INPUT),
        "detection_threshold": (
            "none" if name == "A_all_2000_HVGs"
            else "detected in at least 5% of all cells"
        ),
        "expression_transform": "unchanged from input",
        "scaling": "not applied",
    }

    subset.write_h5ad(output)

    gene_variances = variance[mask]
    lines.extend([
        f"Set: {name}",
        f"Genes: {int(mask.sum())}",
        f"Output: {output}",
        f"All expression values finite: {bool(np.isfinite(subset.X).all())}",
        f"Any negative values: {bool((np.asarray(subset.X) < 0).any())}",
        f"Mean gene variance: {float(gene_variances.mean()):.6f}",
        f"Median gene variance: {float(np.median(gene_variances)):.6f}",
        "",
    ])

summary = pd.DataFrame({
    "gene_symbol": adata.var["gene_symbol"].astype(str).to_numpy(),
    "detection_fraction": detection_fraction,
    "variance": variance,
    "retained_in_set_B": mask_filtered,
})
summary.to_csv("results/pca/feature_set_gene_audit.csv", index=False)

REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("\n".join(lines) + "\n")

print("\n".join(lines))
print("Gene-level audit: results/pca/feature_set_gene_audit.csv")
print("Both feature sets retain the original expression values.")
