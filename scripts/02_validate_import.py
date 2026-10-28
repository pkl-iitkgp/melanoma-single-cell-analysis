from pathlib import Path

import numpy as np
import scanpy as sc

INPUT = Path("data/processed/GSE72056_melanoma.h5ad")
REPORT = Path("results/qc/import_validation_summary.txt")

print("Loading saved AnnData...")
adata = sc.read_h5ad(INPUT)

checks = {
    "Expected cells (4645)": adata.n_obs == 4645,
    "Expected gene rows (23686)": adata.n_vars == 23686,
    "Unique cell identifiers": adata.obs_names.is_unique,
    "Unique AnnData gene names": adata.var_names.is_unique,
    "All expression values finite": bool(np.isfinite(adata.X).all()),
    "No negative expression values": bool((adata.X >= 0).all()),
    "Required annotations present": all(
        col in adata.obs.columns
        for col in [
            "tumor_sample",
            "malignant_status",
            "published_cell_type",
        ]
    ),
}

lines = [
    "GSE72056 IMPORT VALIDATION",
    "=" * 32,
    f"Dimensions: {adata.n_obs} cells x {adata.n_vars} gene rows",
    f"Expression matrix type: {type(adata.X).__name__}",
    f"Expression dtype: {adata.X.dtype}",
    f"Minimum expression: {float(adata.X.min()):.6f}",
    f"Maximum expression: {float(adata.X.max()):.6f}",
    "",
    "Validation checks:",
]

for name, passed in checks.items():
    lines.append(f"[{'PASS' if passed else 'FAIL'}] {name}")

lines.extend([
    "",
    "Cells by malignant status:",
    adata.obs["malignant_status"].value_counts().to_string(),
    "",
    "Cells by tumor sample (first 10):",
    adata.obs["tumor_sample"].value_counts().head(10).to_string(),
    "",
    "Published cell-type annotations:",
    adata.obs["published_cell_type"].value_counts().to_string(),
    "",
    f"Unique original gene symbols: {adata.var['gene_symbol'].nunique()}",
    f"Repeated original gene-symbol rows: "
    f"{adata.var['gene_symbol'].duplicated().sum()}",
])

REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nReport saved to: {REPORT}")

if not all(checks.values()):
    raise SystemExit("Validation failed. Do not proceed to QC yet.")

print("\nVALIDATION SUCCESSFUL")
