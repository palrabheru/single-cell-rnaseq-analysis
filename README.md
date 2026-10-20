# Single-Cell RNA-seq Analysis | Python + Scanpy

An end-to-end, reproducible introductory single-cell RNA-seq workflow with a planned neuro-oncology extension. **The runnable first phase uses public PBMC 3k raw count data**, *not* glioblastoma cells. No cancer-specific claims can be made from the PBMC demo.

## Research questions
1. Which immune cell populations emerge from an unsupervised scRNA-seq analysis?
2. Which marker genes characterize the inferred clusters?
3. How can this workflow be adapted to study immune-cell diversity in glioblastoma?

## Analysis pipeline
Raw 10x PBMC counts → per-cell QC → filtering → library-size normalization → log1p → HVGs → PCA → neighbors → Leiden → UMAP → cluster marker candidates → CSV/plots/H5AD.

## Reproducibility and limitations
- QC thresholds (`>=200` detected genes; `<20%` mitochondrial counts) are illustrative and should be tuned after examining QC distributions. A doublet detection step is not yet implemented.
- The original PBMC data contain counts; normalization is appropriate for this demo. Do **not** apply raw-count QC metrics or count normalization blindly to already normalized TPM matrices.
- Clusters and ranked markers are exploratory. They are not automatically validated cell type assignments. Differential expression between conditions requires donor/sample-aware design, with appropriate biological replicates and multiple-testing controls.
- Files under `data/` and generated `results/` are intentionally ignored in Git.

## Next phase: Glioblastoma
Explore GEO **GSE131928**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928

This study reports 24,131 single cells from 28 glioblastoma patients across Smart-seq2 and 10x assays. Its GEO 10x processed matrix contains **TPM**, not raw UMI counts. Before extending the pipeline, plan memory-efficient loading and TPM-appropriate methods, inspect patient/tumor metadata, distinguish malignant and non-malignant cells with multiple lines of evidence, and account for patient effects. Tumor microenvironment findings must not be inferred from the PBMC example.

## Project status
- [x] Initial Python pipeline and basic automated QC test
- [ ] Run full workflow and review output plots
- [ ] Curate immune-cell annotations
- [ ] Add glioblastoma dataset ingestion and patient-aware analysis
- [ ] Write results narrative with limitations and references

## References
- Scanpy preprocessing and clustering: https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html
- Scanpy installation: https://scanpy.readthedocs.io/en/stable/installation.html
- Glioblastoma series GSE131928: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928
