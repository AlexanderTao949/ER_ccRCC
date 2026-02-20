#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(data.table)
library(ggplot2)
library(grid)
library(scTookit)

#### Load data and process ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)

#### Calculate OR ####
test <- CalculateOR(seu, metacol = "ERScore_group", cellcol = "mps_celltype", metacol_values = c("High_ERScore", "Low_ERScore"))
qs::qsave(test, file = "analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/04_OR/01_OR_test.qs")

p <- PlotOR(test)
print(p)
ggsave(filename = "analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/04_OR/01_OR_test.pdf", p, width = 3.89, height = 4.56)
