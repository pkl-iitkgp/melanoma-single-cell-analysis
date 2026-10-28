from pathlib import Path

import numpy as np
import scanpy as sc

ROOT = Path(".")
INPUT = ROOT / "data/processed/GSE72056_melanoma_hvg_annotated.h5ad"
OUTPUT = ROOT / "data/processed/GSE72056_melanoma_pca_input.h5ad"
REPORT = ROOT / "results/pca/pca_input_summary.txt"

adata = sc.read_h5ad(INPUT)

# Keep only the HVGs selected in the previous step.
hvg_mask = adata.var["highly_variable"].to_numpy()
pca_adata = adata[:, hvg_mask].copy()

# Record provenance and analysis decisions.
pca_adata.uns["pca_preparation"] = {
    "source_file": str(INPUT),
    "feature_selection": "Scanpy highly_variable flag from flavor='seurat'",
    "selected_feature_count": int(hvg_mask.sum()),
    "expression_values": (
        "Original processed GEO expression values retained; "
        "no additional normalization or log transformation applied."
    ),
    "scaling": "Not applied in this preparation step.",
}

# Validate the object before saving.
X = np.asarray(pca_adata.X)
assert X.shape == (adata.n_obs, int(hvg_mask.sum()))
assert np.isfinite(X).all()
assert (X >= 0).all()
assert "malignant_status" in pca_adata.obs
assert "published_cell_type" in pca_adata.obs
assert "tumor_sample" in pca_adata.obs

pca_adata.write_h5ad(OUTPUT)

report = [
    "PCA INPUT PREPARATION SUMMARY",
    "=============================",
    f"Input file: {INPUT}",
    f"Output file: {OUTPUT}",
    f"Cells: {pca_adata.n_obs}",
    f"Selected genes: {pca_adata.n_vars}",
    f"Expression dtype: {X.dtype}",
    f"Expression minimum: {float(X.min()):.6f}",
    f"Expression maximum: {float(X.max()):.6f}",
    f"All values finite: {bool(np.isfinite(X).all())}",
    f"Any negative values: {bool((X < 0).any())}",
    f"Zero fraction: {float((X == 0).mean()):.6f}",
    "Additional normalization/log transformation: none",
    "Centering/scaling: not applied yet",
    f"Malignant status categories: {pca_adata.obs['malignant_status'].value_counts().to_dict()}",
    f"Published cell types: {pca_adata.obs['published_cell_type'].value_counts().to_dict()}",
    f"Tumor samples: {pca_adata.obs['tumor_sample'].nunique()}",
]

REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("\n".join(report) + "\n")

print("\n".join(report))
print("\nSaved PCA input object and summary report.")
