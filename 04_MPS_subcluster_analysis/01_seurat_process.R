#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(harmony)
library(clustree)
library(ggplot2)
library(scTookit)

#### Load data and processing ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthread = 50)
seu <- subset(seu, subset = major_celltype == "Mononuclear_Phagocytes")

#### Seurat process ####
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()
ndim <- pcs_determine(seu)
seu <- seu %>% 
  RunHarmony(reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony", dims.use = 1:ndim) %>% 
  RunUMAP(reduction = "harmony", dims = 1:ndim, reduction.name = "umap")

#### Find clusters ####
seu <- seu %>% 
  FindNeighbors(reduction = "harmony", dims = 1:ndim, graph.name = "MPS_RNA_snn") %>% 
  FindClusters(resolution = seq(0.1, 1, 0.05), graph.name = "MPS_RNA_snn")
p <- clustree(seu, prefix = "MPS_RNA_snn_res.")
ggsave(filename = "analysis/figure/04_MPS_subcluster_analysis/01_seurat_process/01_clustree.pdf", p, width = 15, height = 30)

seu$MPS_subclusters <- seu$MPS_RNA_snn_res.0.35
seu$MPS_subclusters <- factor(seu$MPS_subclusters, levels = as.character(0:8))
qs::qsave(seu, file = "analysis/data/04_MPS_subcluster_analysis/01_seurat_process/01_mps_subclusters", nthread = 50)

p <- DimPlot(seu, reduction = "umap", group.by = "major_celltype", pt.size = 0.05, label = T, raster = FALSE)+NoLegend()+
  ggsci::scale_color_d3("category20")
ggsave("analysis/figure/04_MPS_subcluster_analysis/01_seurat_process/02_mps_major_celltype_umap.png", p, width = 6, height = 6, dpi = 500)

p <- DimPlot(seu, reduction = "umap", group.by = "MPS_subclusters", pt.size = 0.05, label = T, raster = FALSE)+NoLegend()+
  ggsci::scale_color_d3("category20")
ggsave("analysis/figure/04_MPS_subcluster_analysis/01_seurat_process/03_mps_subcluster_umap.png", p, width = 6, height = 6, dpi = 500)

