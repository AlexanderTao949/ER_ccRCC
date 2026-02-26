#### Load packages ####
library(Seurat)

#### Choose Early Cell ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
root.cell <- CellSelector(DimPlot(seu, reduction = "umap", group.by = "mps_celltype"))
root.cell <- "RCC114_CGGCAGTTCGCGTGCA-1"
DimPlot(seu, reduction = "umap", cells.highlight = root.cell)

#### Transfer to Annadata ####
DefaultAssay(seu) <- 'RNA'
seu2 <- seu[rownames(seu[["RNA"]]@scale.data), ]
reticulate::use_condaenv("sceasy")
sceasy::convertFormat(seu2, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = "analysis/data/07_mps_trajectory/02_Palantir/01_mps_scaled_adata.h5ad")

