# ZFP36L2 — Organoid scRNA-seq Analysis

Computational workflow for the **organoid** single-cell RNA-seq analyses in
*ZFP36L2 orchestrates stress-adaptive plasticity during intestinal regeneration
and metastasis* (Jiang et al., 2026). This repository holds the organoid arm of
the study; see the [lab index repo](https://github.com/joechanlab/Jiang_etal_2026_scRNAseq)
for the other components (including the [HTAN patient analyses](https://github.com/mlallo/Jiang_etal_2026_HTAN)).

## Repository structure

The pipeline is organized by analysis stage. Files are grouped as-is (no renaming);
each stage folder keeps its cluster submission scripts (`bsub_*`) alongside the
scripts they launch.

| Folder | Stage |
|---|---|
| `01_qc_preprocessing/` | Batch setup, ambient-RNA removal (CellBender), doublet detection, concatenation, scran normalization, hashtag removal, sample manifests |
| `02_integration_scvi/` | scVI integration and post-processing |
| `03_stratify/` | Stratification and reprocessing of subsets |
| `04_metacells/` | Metacell construction and imputation |
| `05_differential_expression/` | Differential expression (MAST, DESeq2, limma) and DEG plotting |
| `06_gsea_enrichment/` | GSEA (gseapy), Enrichr, leading-edge, gene-network, gene-set (`.gmt`) utilities |
| `07_phenograph_classification/` | PhenoGraph classification, marker genes, ternary/archetype-fraction plots |
| `08_local_correlation_hotspot/` | Hotspot, autocorrelation, local-correlation, kNN similarity |
| `09_abundance/` | Differential abundance (Milo, MELD, Dirichlet regression) |
| `10_archetypal_analysis/` | Archetypal analysis and circular projections |
| `11_cell_scoring_annotation/` | Cell-cycle / gene-set scoring, CellTypist, CellRank, marker dotplots |
| `12_plotting_misc/` | Post-processing plots and exploratory notebooks |

`scrna.yml` is the conda environment specification. `.RData` is a saved R session
used by some of the R-based analyses.

## Environment

```bash
conda env create -f scrna.yml
conda activate scrna
```

## Large files (Git LFS)

Three files exceed GitHub's 100 MB limit and are stored via [Git LFS](https://git-lfs.github.com):

- `.RData`
- `07_phenograph_classification/11b_plotting_phenographclassification.ipynb`
- `09_abundance/differential_abundance_meld.ipynb`

To retrieve their actual contents (rather than small pointer files), install Git LFS
before cloning:

```bash
git lfs install
git clone https://github.com/britneyf/Jiang_etal_2026_Organoid.git
```

## Notes

The workflow was developed and run on an LSF cluster; the `bsub_*` shell scripts are
the corresponding job-submission wrappers. Notebooks and scripts reference paths that
may need to be adapted to your environment.
