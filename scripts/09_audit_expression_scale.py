from pathlib import Path
from datetime import datetime

import numpy as np
import scanpy as sc

ROOT = Path(".")
INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
OUT = ROOT / "results/qc"
OUT.mkdir(parents=True, exist_ok=True)

REPORT = OUT / "expression_scale_audit.txt"

print(f"Loading: {INPUT}")
adata = sc.read_h5ad(INPUT)

X = adata.X

# Compute diagnostics without changing the matrix.
values = X.data if hasattr(X, "data") and not isinstance(X, np.ndarray) else X
values = np.asarray(values).ravel()

finite = np.isfinite(values)
nonzero = values[values != 0]

lines = [
    "GSE72056 EXPRESSION-SCALE AUDIT",
    "=" * 40,
    f"Audit date: {datetime.now().isoformat(timespec='seconds')}",
    f"Input: {INPUT}",
    f"Scanpy version: {sc.__version__}",
    "",
    "DATA DIMENSIONS",
    f"Cells: {adata.n_obs}",
    f"Genes: {adata.n_vars}",
    f"Matrix type: {type(X).__name__}",
    f"Data dtype: {X.dtype}",
    "",
    "NUMERICAL DIAGNOSTICS",
    f"Total matrix entries: {values.size}",
    f"Finite entries: {int(finite.sum())}",
    f"Non-finite entries: {int((~finite).sum())}",
    f"Minimum: {float(values[finite].min()) if finite.any() else 'NA'}",
    f"Maximum: {float(values[finite].max()) if finite.any() else 'NA'}",
    f"Zero fraction: {float(np.mean(values == 0))}",
    f"Negative entries: {int(np.sum(values < 0))}",
    f"Nonzero-entry median: {float(np.median(nonzero)) if nonzero.size else 'NA'}",
    "",
    "PROVENANCE AND INTERPRETATION",
    "GSE72056 is a Smart-seq2 melanoma single-cell RNA-seq dataset.",
    "Published methods describe expression as E = log2(TPM/10 + 1).",
    "The supplied matrix is treated as already transformed expression,",
    "not as raw UMI counts or raw gene counts.",
    "No normalization, log transformation, filtering, or scaling was applied.",
    "",
    "LIMITATIONS",
    "Numerical diagnostics alone cannot prove the transformation formula.",
    "The formula is based on published dataset methods and must be cited.",
    "The original matrix and AnnData object were not modified by this script.",
]

if not finite.all():
    lines.append("WARNING: Non-finite values detected; investigate before analysis.")

if finite.any() and values[finite].min() < 0:
    lines.append("WARNING: Negative values detected; investigate before analysis.")

if adata.n_obs != 4645:
    lines.append("WARNING: Cell count differs from the expected 4645.")

if adata.n_vars != 23686:
    lines.append("WARNING: Gene count differs from the expected 23686.")

REPORT.write_text("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nSaved report: {REPORT}")
