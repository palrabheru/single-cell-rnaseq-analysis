"""Test QC threshold behavior using a tiny synthetic AnnData object."""
import anndata as ad
import pandas as pd
import numpy as np
from src.analyze import qc_mask

def test_qc_mask():
    a = ad.AnnData(np.ones((3, 2)))
    a.obs = pd.DataFrame({'n_genes_by_counts': [250, 150, 300], 'pct_counts_mt': [5.0, 1.0, 25.0]}, index=a.obs_names)
    assert qc_mask(a).tolist() == [True, False, False]
