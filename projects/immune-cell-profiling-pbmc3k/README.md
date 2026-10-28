# Immune Cell Profiling with Single-Cell RNA-seq

A separate mini-project that uses **Python + Scanpy** to explore the immune-cell populations in a public PBMC dataset. The goal is to go from raw gene counts to interpretable clusters, without pretending the computer can identify cell types perfectly.

## Question

Can we find distinct immune-cell populations from thousands of individual cells, and which genes help tell them apart?

## Data

- **Source:** 10x Genomics PBMC 3k, available from Scanpy's built-in dataset loader.
- **Sample:** Peripheral blood mononuclear cells (PBMCs) from one healthy donor.
- **Input:** Around 2,700 cells and 32,738 genes before filtering.
- **Important limitation:** This is **not a cancer dataset**, and this analysis cannot establish cancer biomarkers or clinical conclusions.

Dataset documentation: https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html

## What the pipeline does

1. Loads raw single-cell RNA counts.
2. Filters cells/genes with very little information and calculates mitochondrial RNA fraction.
3. Plots quality-control (QC) metrics and filters high-mitochondrial cells.
4. Normalizes counts, log-transforms them, and selects highly variable genes.
5. Performs PCA, builds a nearest-neighbor graph, creates a UMAP, and runs Leiden clustering.
6. Finds genes enriched in each cluster and plots immune-cell marker expression.
7. Saves a cluster summary, marker gene table, PNG figures, and processed AnnData file.

## Run on Windows (PowerShell)

Requires Python 3.10–3.12 (a clean Python 3.11 environment is a good choice).

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
python analysis.py
```

The first run downloads the public data. Generated results appear under `results/`. This analysis has not been executed in this repository yet, so don't treat example cluster names or counts as confirmed findings.

## Outputs

| File | Meaning |
| --- | --- |
| `results/qc_before.png` | QC distributions before filtering |
| `results/qc_after.png` | QC distributions after filtering |
| `results/umap_clusters.png` | Leiden-defined cell clusters |
| `results/umap_markers.png` | Expression of selected immune markers |
| `results/cluster_sizes.csv` | Number of cells in each cluster |
| `results/cluster_markers.csv` | Ranked candidate marker genes |
| `results/processed_pbmc3k.h5ad` | Processed single-cell data |

## How to interpret results

Clusters reflect similarity in gene-expression profiles; they are not automatically validated cell types. Use multiple markers (e.g. **CD3D**, **MS4A1**, **NKG7**, **LYZ**, **FCGR3A**) together, check for mixed signatures, and compare to published references before giving clusters biological labels.

## Next steps

- Manually annotate each Leiden cluster using a marker-gene panel.
- Examine how different QC settings change clusters.
- Compare these healthy immune populations with tumor-infiltrating immune cells from a cancer dataset.
- Try reference-based cell-type annotation and check confidence.

## Credits

Data: 10x Genomics PBMC 3k. Analysis follows widely used Scanpy workflows, with custom outputs and explicit interpretation limits. Scanpy tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html
