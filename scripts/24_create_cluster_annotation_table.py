from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "results/clustering/A_all_2000_HVGs_leiden_r05_cluster_counts.csv"
OUT = ROOT / "results/markers/cluster_annotations_working.csv"

annotations = [
    (0, "Cytotoxic T-cell-enriched", "High", "CD8A, CD8B, NKG7, PRF1"),
    (1, "CD4/T-cell-enriched", "High", "CD4, IL7R, TCF7, CD3D"),
    (2, "B cells", "High", "MS4A1, CD79A, CD19, PAX5"),
    (3, "Malignant melanoma-associated state", "High", "92.84% malignant; XAGE family"),
    (4, "Malignant melanocytic-state candidate", "Moderate", "80.59% malignant; DCT, TYRP1"),
    (5, "Myeloid/macrophage-enriched", "High", "CSF1R, CD163, FCGR1A, TYROBP"),
    (6, "T-cell-enriched; subtype unresolved", "Moderate", "95.06% non-malignant; mixed marker profile"),
    (7, "Malignant melanoma-associated state candidate", "Moderate", "85.81% malignant; MAGEA4, TERT"),
    (8, "Malignant-associated; state unresolved", "Low", "94.03% malignant; OCA2, KIT"),
    (9, "Malignant-associated; state unresolved", "Low", "89.92% malignant; LY6K"),
    (10, "Malignant melanoma-associated state candidate", "Moderate", "98.41% malignant; MAGEC1, MAGEC2"),
    (11, "Fibroblast/stromal-enriched; mixed", "Moderate", "DPT, LUM, DCN; 37.07% published CAF"),
    (12, "NK-cell-enriched", "High", "GNLY, NCR1, FGFBP2, GZMB"),
    (13, "Proliferating T-cell-enriched", "Moderate", "MKI67, TOP2A, CDK1; published T-cell-enriched"),
    (14, "T-cell-enriched; subtype unresolved", "Moderate", "Published labels predominantly T cells; SNAR-rich"),
    (15, "Endothelial cells", "High", "CLDN5, CDH5, ERG, TIE1"),
    (16, "Malignant-associated; state unresolved", "Low", "88.89% malignant; unusual marker profile"),
    (17, "Endothelial cells", "High", "CDH5, ESAM, ROBO4, PLVAP"),
    (18, "Plasma-cell/B-lineage-enriched", "Moderate", "MZB1, SDC1, IGJ, TNFRSF17"),
    (19, "pDC-like candidate", "Low", "LILRA4, PLD4, FLT3, DNASE1L3; only 19 cells"),
]

columns = [
    "cluster",
    "working_annotation",
    "confidence",
    "supporting_evidence",
]

annotation_df = pd.DataFrame(annotations, columns=columns)

if INPUT.exists():
    counts = pd.read_csv(INPUT)
    print("Cluster-count file columns:", counts.columns.tolist())
    print("Cluster-count file found; inspect its schema before automatic merging.")
else:
    print("Cluster-count file not found; creating annotation table without cell counts.")

annotation_df.to_csv(OUT, index=False)

print("\nSaved:", OUT)
print("\nWorking annotations:")
print(annotation_df.to_string(index=False))
