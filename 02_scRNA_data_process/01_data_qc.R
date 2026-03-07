#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(decontX)
library(ggplot2)
library(scTookit)
library(SingleCellExperiment)
library(reticulate)
use_condaenv("scrublet")

#### Load data and process ####
dir <- "./scRNA_upstream/cellranger_results/"
files <- list.files(dir)
seulist <- pbmcapply::pbmclapply(files, function(x){
  path <- file.path(dir, x, "outs/filtered_feature_bc_matrix")
  scdata <- Read10X(path)
  colnames(scdata) <- paste(x, colnames(scdata), sep = "_")
  seu <- CreateSeuratObject(scdata, project = x, min.cells=3, min.features=200)
  return(seu)
})
names(seulist) <- files
qs::qsave(seulist, file = "analysis/data/02_scRNA_data_process/01_data_qc/01_raw_seulist.qs", nthreads = 50)

#### RNA contamination, doublet, QC indicator calculation ####
seulist_qc <- list()
for (i in files){
  seu <- seulist[[i]]
  seu <- CalcRNAContam(seu)
  seu <- MarkDoublets(seu)
  seu <- calcQCMetrics(seu)
  seu <- subset(seu, subset = percent_mito_by_sample < 30)
  seulist_qc[[i]] <- seu
}
qs::qsave(seulist_qc, file = "analysis/data/02_scRNA_data_process/01_data_qc/02_seulist_qc.qs", nthreads = 50)

#### Overall view ####
metrics <- c("nCount_RNA", "nFeature_RNA", "contamination_score", 
             "doublet_score", "percent_mito_by_sample", "percent_ribo_by_sample", 
             "percent_hb_by_sample", "percent_hsp_by_sample", 
             "percent_dissociation_by_sample")
plots_list <- lapply(metrics, function(x){p <- PlotMultiQC(seulist_qc, x, y_limits = NULL)})
combined_plot <- patchwork::wrap_plots(plots_list, ncol = 3)
ggsave("analysis/figure/02_scRNA_data_process/01_data_qc/01_qc_plot_all.pdf", combined_plot, width = 15, height = 15)

seulist_qc <- seulist_qc[setdiff(names(seulist_qc), c("RCC101", "RCC115"))]
qs::qsave(seulist_qc, file = "analysis/data/02_scRNA_data_process/01_data_qc/02_seulist_qc_filter.qs", nthreads = 50)

#### Filter ####
seu_new_list <- NULL
for(i in names(seulist_qc)){
  cat(i, "\n")
  seu <- seulist_qc[[i]]
  seu <- sample_qc(seu, plotdir = "analysis/figure/02_scRNA_data_process/01_data_qc/02_qc_samples_plot", filename = i)
  seu_new_list[[i]] <- seu
}
qs::qsave(seu_new_list, file = "analysis/data/02_scRNA_data_process/01_data_qc/03_seulist_clean.qs", nthreads = 50)

#### Merge data ####
seu <- merge(x = seu_new_list[[1]], y = seu_new_list[-1])
qs::qsave(seu, file = "analysis/data/02_scRNA_data_process/01_data_qc/04_seu_merge.qs", nthreads = 50)
