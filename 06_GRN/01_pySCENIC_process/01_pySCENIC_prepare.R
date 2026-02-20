#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(arrow)
library(SCopeLoomR)

#### Load data and process ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
cisdb <- arrow::read_feather("analysis/resource/cisTargetDB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

counts <- as.matrix(seu@assays$RNA@counts)
genes.use <- intersect(colnames(cisdb), rownames(counts))
mc.mat <- counts[genes.use, ]

loom <- SCopeLoomR::build_loom(
  file.name         = "analysis/data/06_GRN/01_pySCENIC_process/mc_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()
