"""Explore immune-cell populations in the public 10x PBMC 3k dataset.

Run: python analysis.py
Results are written under results/. No cell-type labels are assigned automatically.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=110, facecolor="white")
sc.settings.figdir = OUT

# PBMC 3k: one healthy donor; not a tumor sample.
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
print(f"Raw: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Basic filters before computing quality metrics.
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

def save_qc(plot_name: str) -> None:
    """Save the distributions without assuming any particular cell count."""
    cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
    for axis, col in zip(axes, cols):
        axis.hist(adata.obs[col], bins=45)
        axis.set_xlabel(col)
        axis.set_ylabel("Cells")
    fig.tight_layout()
    fig.savefig(OUT / plot_name, dpi=160)
    plt.close(fig)

save_qc("qc_before.png")

# These thresholds are starting points, not universal biological rules.
adata = adata[(adata.obs["n_genes_by_counts"] < 2500) & (adata.obs["pct_counts_mt"] < 5)].copy()
if adata.n_obs < 50:
    raise RuntimeError("QC thresholds removed too many cells. Inspect QC plots.")
save_qc("qc_after.png")
print(f"After QC: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# Preserve raw counts to make follow-up work reproducible.
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()  # full normalized genes for marker visualization and DE
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)

sc.pl.umap(adata, color="leiden", show=False, return_fig=False)
plt.savefig(OUT / "umap_clusters.png", dpi=170, bbox_inches="tight")
plt.close("all")

markers = ["CD3D", "MS4A1", "NKG7", "LYZ", "FCGR3A", "PPBP"]
available = [gene for gene in markers if gene in adata.raw.var_names]
if available:
    sc.pl.umap(adata, color=available, use_raw=True, ncols=3, show=False)
    plt.savefig(OUT / "umap_markers.png", dpi=170, bbox_inches="tight")
    plt.close("all")

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)
marker_table = sc.get.rank_genes_groups_df(adata, group=None)
marker_table.to_csv(OUT / "cluster_markers.csv", index=False)
counts = adata.obs["leiden"].value_counts().sort_index()
counts.rename_axis("cluster").rename("n_cells").to_csv(OUT / "cluster_sizes.csv")
adata.write_h5ad(OUT / "processed_pbmc3k.h5ad", compression="gzip")
print("Saved results in:", OUT)
print(counts.to_string())
