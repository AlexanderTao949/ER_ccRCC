#### Load packages ####
library(Seurat)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)

#### Transfer to Annadata ####
DefaultAssay(seu) <- 'RNA'
reticulate::use_condaenv("sceasy")
sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts", 
                      outFile = "analysis/data/07_mps_trajectory/03_scTour/01_mps_adata.h5ad")
