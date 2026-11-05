# RAPHE_Code

Scripts for the analysis of data published in:

**Serotonergic innervation and cortical progenitor regulation in human brain assembloids**
Pérez Fernández, Gasparotto, et al. (2026)

This repository contains the R, Python and Jupyter code used for all computational analyses in the study. Scripts are numbered in the order in which they were run.

## Repository structure

```
RAPHE_Code/
├── Sequencing/
│   ├── scripts/
│   │   ├── 001_QC.Rmd                                  QC, filtering and merging of the CO/RO control and fused libraries (Seurat); export for scVI
│   │   ├── 002_Integration.ipynb                       scVI integration and CellTypist annotation (Developing_Human_Brain model)
│   │   ├── 003_Annotation_analysis.Rmd                 Cell-type annotation, composition, DEG and GO enrichment, cell cycle (tricycle),
│   │   │                                               CellChat cell-cell communication, pseudotime (monocle3)
│   │   ├── 004_Quantifications.Rmd                     qPCR and immunostaining quantifications (Figures 1 and 4)
│   │   └── 005_PH3_Quantification_ROCO_Antagonist.Rmd  PH3+ proliferation in fused organoids treated with serotonin receptor antagonists (KET, VOL)
│   ├── helpers/
│   │   ├── Set_environment.R                           packages and colour palettes sourced by the scripts
│   │   └── Set_functions.R                             QC and GO plotting functions sourced by the scripts
│   ├── rawdata/Quantifications/                        qPCR and immunostaining quantification tables read by 004 and 005
│   ├── data/                                           intermediate objects written by the scripts (not tracked)
│   └── Figures/                                        figure output (not tracked)
└── sDarken/                                            sDarken serotonin sensor imaging analysis (see sDarken/README.md)
    ├── code/calcimetrics/                              Python package for ImageJ/Fiji ROI fluorescence exports
    ├── notebooks/                                      acute slice, fused organoid and dissociated neuron analyses
    ├── input_data/                                     place the source data here (see Data availability)
    └── figures/                                        figure output (not tracked)
```

All R scripts use paths relative to `Sequencing/scripts/` and should be run from that directory.

## Data availability

The RNA sequencing data generated in this study have been deposited in the German Human Genome-Phenome Archive (GHGA) under accession code GHGAS76412448666741. Cell Ranger outputs (`filtered_feature_bc_matrix`) for each library are expected in `Sequencing/rawdata/<sample>/` (CO_CTRL, CO_FUSED, RO_CTRL, RO_FUSED).

Source data for all other analyses are provided with the paper.

## Software

- R with Seurat, SeuratWrappers, monocle3, CellChat, tricycle, clusterProfiler, org.Hs.eg.db, ComplexHeatmap, circlize, ggplot2, patchwork, ggpubr, ggalluvial, rstatix, dplyr, tidyr, readxl, openxlsx
- Python with scanpy, scvi-tools, CellTypist (integration); numpy, pandas, scipy, matplotlib, seaborn, SciencePlots (sDarken, installed via `pip install -e sDarken/code/`)

## Contributors

RPF (Raquel Pérez Fernández), MG (Matteo Gasparotto)

## License

MIT, see [LICENSE](LICENSE).
