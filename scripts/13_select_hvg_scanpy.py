from pathlib import Path

import scanpy as sc
import pandas as pd

ROOT = Path(".")
INPUT = ROOT / "data/processed/GSE72056_melanoma.h5ad"
OUTPUT_H5AD = ROOT / "data/processed/GSE72056_melanoma_hvg_annotated.h5ad"
OUTPUT_CSV = ROOT / "results/qc/scanpy_hvg_top2000.csv"
OUTPUT_REPORT = ROOT / "results/qc/scanpy_hvg_summary.txt"

OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

print(f"Loading original object: {INPUT}")
adata = sc.read_h5ad(INPUT)

print("\n===== INPUT CHECK =====")
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)
print("Expression dtype:", adata.X.dtype)
print("Expression note:", adata.uns.get("expression_note", "Not recorded"))

if adata.n_obs != 4645 or adata.n_vars != 23686:
    raise ValueError("Unexpected input dimensions; stopping before analysis.")

if not adata.var_names.is_unique:
    raise ValueError("Gene identifiers are not unique; stopping.")

print("\n===== SELECTING HIGHLY VARIABLE GENES =====")
print("Method: Scanpy dispersion-based Seurat flavor")
print("Requested genes: 2000")
print("No additional normalization or log transformation will be applied.")

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000,
    subset=False,
    inplace=True,
)

hvg = adata.var.loc[adata.var["highly_variable"]].copy()
hvg.insert(0, "gene_id", hvg.index.astype(str))

columns = [
    col for col in [
        "gene_id",
        "gene_symbol",
        "means",
        "dispersions",
        "dispersions_norm",
        "highly_variable",
    ]
    if col in hvg.columns
]

hvg = hvg[columns].sort_values(
    "dispersions_norm", ascending=False
)

hvg.to_csv(OUTPUT_CSV, index=False)
adata.write_h5ad(OUTPUT_H5AD, compression="gzip")

report = [
    "SCANPY HIGHLY VARIABLE GENE SELECTION",
    "=" * 42,
    f"Scanpy version: {sc.__version__}",
    f"Input: {INPUT}",
    f"Cells: {adata.n_obs}",
    f"Genes before selection: {adata.n_vars}",
    f"Selected HVGs: {int(adata.var['highly_variable'].sum())}",
    "Method: flavor='seurat'",
    "Requested number: 2000",
    "Expression transformation: original values retained",
    "Additional normalization/log transformation: none",
    "",
    "IMPORTANT INTERPRETATION",
    "HVGs reflect variation across the included cell population.",
    "Variation can reflect cell identity, malignant state, tumor of origin,",
    "and technical effects; HVG status alone does not establish biology.",
    "",
    f"HVG table: {OUTPUT_CSV}",
    f"Annotated object: {OUTPUT_H5AD}",
]

OUTPUT_REPORT.write_text("\n".join(report) + "\n")

print("\n===== RESULTS =====")
print("Selected HVGs:", int(adata.var["highly_variable"].sum()))
print("\nTop 20 HVGs by normalized dispersion:")
print(hvg.head(20).to_string(index=False))
print("\nSaved:", OUTPUT_CSV)
print("Saved:", OUTPUT_H5AD)
print("Saved:", OUTPUT_REPORT)
print("\nOriginal input object was not modified.")
