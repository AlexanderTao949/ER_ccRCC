#### Load packages and functions ####
library(tidyverse)
library(AUCell)
library(Seurat)
source("resource/R/mcAUCell.R")
regulon2list <- function(regulons){
  rg.names <- unique(regulons$term)
  regulon.list <- lapply(rg.names, function(rg) {
    subset(regulons, term == rg)$gene
  })
  names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
  return(regulon.list)
}


#### Load data and process ####
seu <- qs::qread("data/07_Mono_subcluster_annotation/03_mono_seu_anno.qs", nthreads = 25)
regulons <- clusterProfiler::read.gmt("data/11_GRN_Mono/01_pySCENIC_process/pyscenic_TF.regulons.gmt")

regulons_list <- regulon2list(regulons)
qs::qsave(regulons_list, file = "data/11_GRN_Mono/01_pySCENIC_process/regulons_list.qs")


#### AUCell ####
matrix <- seu@assays$RNA@data %>% as.matrix()
res <- Run_mc_AUCell(matrix, regulons_list)
seu[["SCENIC"]] <- CreateAssayObject(data = res)
qs::qsave(seu, "data/11_GRN_Mono/01_pySCENIC_process/mono_seu_pyscenic.qs")
