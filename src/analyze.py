"""Reproducible exploratory scRNA-seq analysis of 10x PBMC 3k counts.

Run: python -m src.analyze
Note: Results are exploratory; automated cell labels must be manually reviewed.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc

OUT = Path('results')
FIG = OUT / 'figures'
MARKERS = {
    'T cells': ['CD3D', 'CD3E'],
    'B cells': ['MS4A1', 'CD79A'],
    'NK cells': ['NKG7', 'GNLY'],
    'Monocytes': ['LYZ', 'S100A8'],
    'Dendritic cells': ['FCER1A', 'CST3'],
    'Platelets': ['PPBP', 'PF4'],
}

def qc_mask(adata, min_genes=200, max_mt=20.0):
    """Return cell-retention mask; percent mitochondrial counts on 0-100 scale."""
    return (adata.obs['n_genes_by_counts'] >= min_genes) & (adata.obs['pct_counts_mt'] < max_mt)

def save_plot(path):
    plt.savefig(path, dpi=180, bbox_inches='tight')
    plt.close('all')

def main():
    FIG.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white')
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    cells_before, genes_before = adata.shape
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.25, multi_panel=True, show=False)
    save_plot(FIG / '01_qc_violin.png')
    adata = adata[qc_mask(adata)].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_obs < 10 or adata.n_vars < 100:
        raise ValueError('Too few cells/genes after QC; inspect filtering thresholds')
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    adata = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden', flavor='igraph', n_iterations=2)
    sc.pl.umap(adata, color='leiden', show=False)
    save_plot(FIG / '02_umap_clusters.png')
    markers = [g for genes in MARKERS.values() for g in genes if g in adata.raw.var_names]
    if markers:
        sc.pl.dotplot(adata, markers, groupby='leiden', use_raw=True, show=False)
        save_plot(FIG / '03_marker_dotplot.png')
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True)
    de = sc.get.rank_genes_groups_df(adata, group=None)
    de.to_csv(OUT / 'cluster_marker_candidates.csv', index=False)
    adata.obs[['leiden']].to_csv(OUT / 'cell_clusters.csv')
    adata.write_h5ad(OUT / 'pbmc3k_processed.h5ad', compression='gzip')
    summary = {
        'dataset': 'Scanpy PBMC 3k (public 10x raw counts)',
        'cells_before_qc': int(cells_before),
        'genes_before_qc': int(genes_before),
        'cells_after_qc': int(adata.n_obs),
        'genes_hvg': int(adata.n_vars),
        'n_clusters': int(adata.obs['leiden'].nunique()),
        'note': 'Exploratory clusters; check QC, doublets, batch effects, annotations, and biological replicates before drawing conclusions.',
    }
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
