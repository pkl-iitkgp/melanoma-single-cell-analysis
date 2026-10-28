import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

INPUT = Path("data/raw/GSE72056_melanoma_single_cell_revised_v2.txt.gz")
OUTPUT = Path("data/processed/GSE72056_melanoma.h5ad")

if not INPUT.exists():
    raise FileNotFoundError(f"Dataset not found: {INPUT}")

print("Step 1: Inspecting dataset...")

with gzip.open(INPUT, "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    tumor_row = next(reader)
    malignant_row = next(reader)
    cell_type_row = next(reader)

    n_cells = len(header) - 1
    n_genes = sum(1 for _ in reader)

print(f"Cells: {n_cells:,}")
print(f"Gene rows: {n_genes:,}")

print("\nStep 2: Reading expression values...")

gene_symbols = []
expression = np.empty((n_genes, n_cells), dtype=np.float32)

with gzip.open(INPUT, "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")

    for _ in range(4):
        next(reader)

    for i, row in enumerate(reader):
        if i >= n_genes:
            raise ValueError("More gene rows than expected.")

        if len(row) != n_cells + 1:
            raise ValueError(f"Unexpected column count on gene row {i + 1}.")

        gene_symbols.append(row[0])
        expression[i, :] = np.asarray(row[1:], dtype=np.float32)

        if (i + 1) % 5000 == 0:
            print(f"  Loaded {i + 1:,} / {n_genes:,} genes")

if len(gene_symbols) != n_genes:
    raise ValueError("Gene count changed between reading passes.")

if not np.isfinite(expression).all():
    raise ValueError("Expression contains NaN or infinite values.")

print("\nStep 3: Checking values...")
print("Minimum:", float(expression.min()))
print("Maximum:", float(expression.max()))
print("Fraction of zero values:", round(float(np.mean(expression == 0)), 4))

print("\nStep 4: Building cell metadata...")

obs = pd.DataFrame(
    {
        "tumor_sample": tumor_row[1:],
        "malignant_code": malignant_row[1:],
        "cell_type_code": cell_type_row[1:],
    },
    index=header[1:],
)

malignant_map = {
    "0": "unresolved",
    "1": "non-malignant",
    "2": "malignant",
}

cell_type_map = {
    "0": "not assigned",
    "1": "T cell",
    "2": "B cell",
    "3": "macrophage",
    "4": "endothelial cell",
    "5": "CAF",
    "6": "NK cell",
}

obs["malignant_status"] = (
    obs["malignant_code"].map(malignant_map).fillna("unknown")
)
obs["published_cell_type"] = (
    obs["cell_type_code"].map(cell_type_map).fillna("unknown")
)

var = pd.DataFrame({"gene_symbol": gene_symbols})

print("\nStep 5: Creating AnnData...")

adata = sc.AnnData(
    X=expression.T,
    obs=obs,
    var=var,
)

adata.obs_names_make_unique()
adata.var_names_make_unique(join="-dup")

adata.uns["dataset"] = "GEO GSE72056"
adata.uns["expression_note"] = (
    "Processed GEO expression matrix. Original values retained; "
    "no additional normalization or log transformation applied."
)

print("\nStep 6: Saving AnnData...")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUTPUT, compression="gzip")

print("\nIMPORT COMPLETE")
print(adata)
print("\nMalignant status:")
print(adata.obs["malignant_status"].value_counts().to_string())
print("\nPublished cell types:")
print(adata.obs["published_cell_type"].value_counts().to_string())
print("\nSaved file:", OUTPUT)
print("File size (MB):", round(OUTPUT.stat().st_size / 1_000_000, 2))
