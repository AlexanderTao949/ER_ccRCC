#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(scTookit)
library(presto)

#### Load data and processe ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/02_seurat_process/01_seu_clusters.qs", nthreads = 50)
Idents(seu) <- seu$seurat_clusters
metadata <- seu@meta.data

#### ACT ####
markers <- wilcoxauc(seu,
                     group_by = "seurat_clusters",
                     seurat_assay = "RNA",
                     assay = "data",
                     verbose = TRUE)
markers <- markers %>% dplyr::filter(logFC > 0)
markers$pct_diff <- markers$pct_in - markers$pct_out
qs::qsave(markers,file = "analysis/data/02_scRNA_data_process/03_annotation/01_findmarkers.qs")

topmarkers <- markers %>% 
  group_by(group) %>% 
  dplyr::filter(auc > 0.7, padj < 0.05, pct_diff > 30) %>% 
  dplyr::arrange(desc(logFC), .by_group = TRUE) %>% 
  slice_head(n = 30) %>% 
  ungroup()
topmarkers$group <- factor(topmarkers$group, levels = as.character(0:12))

act_input <- topmarkers %>% 
  group_by(group) %>% 
  summarise(
    marker = paste(feature, collapse = ","),
    .groups = "drop"
  ) %>% 
  mutate(line = paste0("cluster_", group, ":", marker))
writeLines(act_input$line, con = "analysis/data/02_scRNA_data_process/03_annotation/02_act_input.txt")

cluster_annotation <- c(
  "0" = "Epithelial_Cells",
  "1" = "CD4_T_Cells",
  "2" = "Natural_Killer_Cells",
  "3" = "CD8_T_Cells",
  "4" = "Epithelial_Cells",
  "5" = "Mononuclear_Phagocytes",
  "6" = "Fibroblast",
  "7" = "Endothelial_Cells",
  "8" = "Mononuclear_Phagocytes",
  "9" = "unknown",
  "10" = "B_Cells",
  "11" = "Endothelial_Cells",
  "12" = "Mast_Cells"
)
metadata$act_annotation <- cluster_annotation[as.character(metadata$seurat_clusters)]
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "act_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/01_act_annotation.png",p ,width = 6, height = 6, dpi = 500)

#### First level annotation ####
markerlist <- list("Immune" = "PTPRC",
                   "Epithelial" = c("EPCAM", "KRT19", "KRT18", "KRT8", "PROM1", "ALDH1A1", "CD24"),
                   "Fibroblasts" = c("COL1A1", "COL3A1", "ACAT2", "PDGFRB", "RGS5", "FGF7", "MME", "DCN", "LUM", "GSN"),
                   "Endothelial" = c("PECAM1", "VWF", "PLVAP"))

p <- MarkDotplot(markerlist, seu, facet = TRUE)
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/02_first_level_annotation_dotplot.pdf", p, width = 12, height = 4)

cluster_annotation <- c(
  "0" = "Epithelial_Cells",
  "1" = "Immune_Cells",
  "2" = "Immune_Cells",
  "3" = "Immune_Cells",
  "4" = "Epithelial_Cells",
  "5" = "Immune_Cells",
  "6" = "Fibroblast",
  "7" = "Endothelial_Cells",
  "8" = "Immune_Cells",
  "9" = "Immune_Cells",
  "10" = "Immune_Cells",
  "11" = "Endothelial_Cells",
  "12" = "Immune_Cells"
)
metadata$first_level_annotation <- cluster_annotation[as.character(metadata$seurat_clusters)]
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "first_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/02_first_level_annotation_umap.png",p ,width = 6, height = 6, dpi = 500)

#### Second level annotation ####
seu_immune <- subset(seu, subset = first_level_annotation == "Immune_Cells")
Idents(seu_immune) <- seu_immune$seurat_clusters

markerlist <- list("CD4_T" = c("CD4", "IL7R"),
                   "CD8_T" = c("CD8A", "CD8B"),
                   "B" = c("MS4A1", "CD19", "CD79A", "IGHG3"),
                   "NK" = c("NKG7", "GNLY", "KLRD1", "FGFBP2", "CX3CR1"),
                   "Mono" = c("S100A8", "S100A9", "LYZ", "CD14","ITGAL","ITGAX","ITGB2"),
                   "Macro" = c("APOE", "C1QA", "C1QB", "FCGR1A","CD68", "CD163","MRC1","ITGAM"),
                   "DC" = c("CD1C", "CD1E", "FCER1A", "LILRA4", "TPM2"),
                   "Mast"= c("TPSAB1", "TPSB2", "CPA3"))

p <- MarkDotplot(markerlist, seu_immune, facet = T)
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/03_second_level_annotation_dotplot.pdf", p, width = 12, height = 4)

cluster_annotation <- c(
  "0" = "Epithelial_Cells",
  "1" = "CD4_T_Cells",
  "2" = "Natural_Killer_Cells",
  "3" = "CD8_T_Cells",
  "4" = "Epithelial_Cells",
  "5" = "Mononuclear_Phagocytes",
  "6" = "Fibroblast",
  "7" = "Endothelial_Cells",
  "8" = "Mononuclear_Phagocytes",
  "9" = "unknown",
  "10" = "B_Cells",
  "11" = "Endothelial_Cells",
  "12" = "Mast_Cells"
)
metadata$second_level_annotation <- cluster_annotation[as.character(metadata$seurat_clusters)]
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "second_level_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/03_second_level_annotation_umap.png",p ,width = 6, height = 6, dpi = 500)


seu <- subset(seu, subset = second_level_annotation == "unknown", invert = T)
seu$major_celltype <- seu$second_level_annotation
seu$major_celltype[seu$major_celltype == "Epithelial_Cells"] <- "Tumor_Cells"
p <- DimPlot(seu, reduction = "umap", group.by = "major_celltype", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/02_scRNA_data_process/03_annotation/04_major_celltype.png", p, width = 6, height = 6, dpi = 500)

qs::qsave(seu, file = "analysis/data/02_scRNA_data_process/03_annotation/03_seu_anno.qs", nthreads = 50)

