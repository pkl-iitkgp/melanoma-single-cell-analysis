from pathlib import Path

import numpy as np
import scanpy as sc

inputs = {
    "Set A": "results/umap/A_all_2000_HVGs_umap.h5ad",
    "Set B": "results/umap/B_HVGs_detected_in_at_least_5pct_umap.h5ad",
}

for label, path in inputs.items():
    print(f"\n===== {label} =====")
    adata = sc.read_h5ad(path)

    print("Shape:", adata.shape)
    print("UMAP coordinates:", adata.obsm["X_umap"].shape)
    print(
        "UMAP coordinates finite:",
        bool(np.isfinite(adata.obsm["X_umap"]).all()),
    )
    print("Neighbor metadata:", adata.uns.get("neighbors", {}))
    print("UMAP metadata:", adata.uns.get("umap", {}))

    graph = adata.obsp.get("connectivities")
    if graph is not None:
        print("Graph shape:", graph.shape)
        print("Graph nonzero entries:", graph.nnz)
        print("Graph symmetry:", (graph != graph.T).nnz == 0)
        print("Graph values finite:", bool(np.isfinite(graph.data).all()))
