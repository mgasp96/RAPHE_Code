# Packages and colour palettes shared by the Sequencing/scripts notebooks.
# Minimal version: only what the scripts in this repository use.

# Single-cell analysis
library(Seurat)

# Visualisation
library(ggplot2)
library(patchwork)

# Data manipulation
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)

# Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Palettes
rna_palette <- colorRampPalette(c("#C2E4EF", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"))
rna_palette_2 <- colorRampPalette(c("#cdcdcd", "#f7bbb1", "#cf796f", "#9f3032", "#6e1720", "#480310"))
diverging_palette <- colorRampPalette(c("#0057B7", "#4FA8F0", "#F0F0F0", "#F58B3C", "#D65A00"))
diverging_palette_4 <- colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", ""))
