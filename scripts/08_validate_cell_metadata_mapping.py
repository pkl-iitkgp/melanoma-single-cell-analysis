from pathlib import Path
import re

import pandas as pd
import scanpy as sc

META = Path("data/metadata/GSE72056_sample_metadata.csv")
H5AD = Path("data/processed/GSE72056_melanoma.h5ad")
OUT = Path("data/metadata")
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(META, dtype=str).fillna("")
adata = sc.read_h5ad(H5AD)

# GEO titles should correspond to the original cell identifiers.
meta["cell_id"] = meta["title"].str.strip()

print("GEO metadata records:", len(meta))
print("AnnData cells:", adata.n_obs)
print("Unique GEO cell titles:", meta["cell_id"].nunique())
print("Unique AnnData cell IDs:", adata.obs_names.nunique())

if meta["cell_id"].duplicated().any():
    print("\nWARNING: duplicated GEO cell titles detected.")
    print(meta.loc[meta["cell_id"].duplicated(False), "cell_id"].head(20))
else:
    print("\nGEO cell titles are unique.")

geo_ids = set(meta["cell_id"])
adata_ids = set(adata.obs_names.astype(str))

print("\nCell IDs present in both:", len(geo_ids & adata_ids))
print("GEO-only IDs:", len(geo_ids - adata_ids))
print("AnnData-only IDs:", len(adata_ids - geo_ids))

# Extract tumor-of-origin from the GEO title, then use the
# explicit GEO "tumor index" characteristic as a fallback.
def extract_tumor(title, characteristics):
    match = re.match(r"(?i)^CY(\d+)", str(title).strip())
    if match:
        return match.group(1)

    match = re.search(
        r"(?i)tumor\s*index\s*:\s*(\d+)",
        str(characteristics)
    )
    return match.group(1) if match else ""

meta["geo_tumor_code"] = meta.apply(
    lambda row: extract_tumor(row["cell_id"], row["characteristics"]),
    axis=1
)

# Preserve cell IDs and compare GEO-derived tumor codes with existing labels.
obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["adata_tumor_code"] = obs["tumor_sample"].astype(str).str.replace(
    r"\\.0$", "", regex=True
)

merged = obs.reset_index(drop=True).merge(
    meta[["cell_id", "geo_accession", "title", "source_name",
          "characteristics", "geo_tumor_code"]],
    on="cell_id",
    how="left",
    validate="one_to_one",
    indicator=True,
)

merged["tumor_code_match"] = (
    merged["adata_tumor_code"] == merged["geo_tumor_code"]
)

print("\nMetadata merge status:")
print(merged["_merge"].value_counts().to_string())

print("\nGEO-derived tumor codes:")
print(merged["geo_tumor_code"].value_counts(dropna=False).sort_index().to_string())

print("\nTumor code agreement:")
print(merged["tumor_code_match"].value_counts(dropna=False).to_string())

print("\nExamples of mismatches:")
cols = ["cell_id", "adata_tumor_code", "geo_tumor_code"]
print(
    merged.loc[~merged["tumor_code_match"], cols]
    .head(20)
    .to_string(index=False)
)

merged.to_csv(OUT / "GSE72056_cell_metadata_merged.csv", index=False)

print("\nSaved merged metadata:")
print(OUT / "GSE72056_cell_metadata_merged.csv")
