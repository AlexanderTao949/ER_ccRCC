#### Load packages and functions ####
library(tidyverse)
library(AUCell)
library(Seurat)
Run_mc_AUCell <- function(expr.matrix, gene.sets, batch.size = 500, ncores = 10){
  if (!is.list(gene.sets)){
    stop("gene.sets must be a list!")
  }
  
  n.cells <- ncol(expr.matrix)
  batches <- floor((1:n.cells - 1) / batch.size)
  batch.levels <- unique(batches)
  
  fun <- function(x){
    tmp <- expr.matrix[, batches == x]
    au <- AUCell::AUCell_buildRankings(tmp, nCores=1, plotStats=F, verbose = T)
    cr <- AUCell::AUCell_calcAUC(gene.sets, au, nCores=1, verbose = T)
    dat <- AUCell::getAUC(cr)
  }
  result <- pbmcapply::pbmclapply(batch.levels, fun, mc.cores = ncores)
  table <- do.call(cbind, result)
}
regulon2list <- function(regulons){
  rg.names <- unique(regulons$term)
  regulon.list <- lapply(rg.names, function(rg) {
    subset(regulons, term == rg)$gene
  })
  names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
  return(regulon.list)
}


#### Load data and process ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 25)
regulons <- clusterProfiler::read.gmt("analysis/data/06_GRN/01_pySCENIC_process/pyscenic_TF.regulons.gmt")

regulons_list <- regulon2list(regulons)
qs::qsave(regulons_list, file = "analysis/data/06_GRN/01_pySCENIC_process/regulons_list.qs")

#### AUCell ####
matrix <- seu@assays$RNA@data %>% as.matrix()
res <- Run_mc_AUCell(matrix, regulons_list)
seu[["SCENIC"]] <- CreateAssayObject(data = res)
qs::qsave(res, "analysis/data/06_GRN/01_pySCENIC_process/pyscenic_auc.qs", nthreads = 50)
qs::qsave(seu, "analysis/data/06_GRN/01_pySCENIC_process/mps_seu_pyscenic.qs", nthreads = 50)
