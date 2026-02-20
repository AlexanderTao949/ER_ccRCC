#### Load packages ####
library(Seurat)
library(tidyverse)
library(scTookit)
library(presto)

#### load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/01_seurat_process/01_mps_subclusters", nthreads = 50)
Idents(seu) <- seu$MPS_subclusters
metadata <- seu@meta.data

#### ACT ####
markers <- wilcoxauc(seu,
                     group_by = "MPS_subclusters",
                     seurat_assay = "RNA",
                     assay = "data",
                     verbose = TRUE)
markers <- markers %>% dplyr::filter(logFC > 0)
markers$pct_diff <- markers$pct_in - markers$pct_out
qs::qsave(markers,file = "analysis/data/04_MPS_subcluster_analysis/02_annotation/01_findmarkers.qs")

topmarkers <- markers %>% 
  group_by(group) %>% 
  dplyr::filter(auc > 0.6, padj < 0.05, pct_diff > 10) %>% 
  dplyr::arrange(desc(logFC), .by_group = TRUE) %>% 
  slice_head(n = 30) %>% 
  ungroup()
topmarkers$group <- factor(topmarkers$group, levels = as.character(0:8))

act_input <- topmarkers %>% 
  group_by(group) %>% 
  summarise(
    marker = paste(feature, collapse = ","),
    .groups = "drop"
  ) %>% 
  mutate(line = paste0("cluster_", group, ":", marker))
writeLines(act_input$line, con = "analysis/data/04_MPS_subcluster_analysis/02_annotation/02_act_input.txt")

# 0: Macrophages
# 1: Macrophages
# 2: Monocytes
# 3: Monocytes
# 4: Dendritic_cell
# 5: unknown
# 6: unknown
# 7: Dendritic_cell
# 8: unknown

cluster_annotation <- c(
  "0" = "Macrophages",
  "1" = "Macrophages",
  "2" = "Monocytes",
  "3" = "Monocytes",
  "4" = "Dendritic_cell",
  "5" = "unknown",
  "6" = "unknown",
  "7" = "Dendritic_cell",
  "8" = "unknown"
)

metadata$mps_act_annotation <- cluster_annotation[as.character(metadata$MPS_subclusters)]
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "mps_act_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/04_MPS_subcluster_analysis/02_annotation/01_act_annotation.png",p ,width = 6, height = 6, dpi = 500)

#### Classical Markers ####
markerlist <- list(
  Mono_CD14 = c("CD14","CSF3R","S100A9","S100A8","S100A12","FCN1"),
  Mono_CD16 = c("FCGR3A", "LST1", "LILRB2","CX3CR1","CDKN1C","CSF1R","ITGAL"),
  IFN_TAMs  = c("CXCL10", "CD274", "ISG15", "CD86", "MHCII"),
  Reg_TAMs  = c("ARG1", "MRC1", "CX3CR1"),
  Inflam_TAMs = c("IL1B", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CCL3", "CCL3L1"),
  Angio_TAMs = c("VEGFA", "SPP1", "VCAN", "FCN1", "THBS1"),
  LA_TAMs = c("APOC1", "APOE", "ACP5", "FABP5"),
  Prolif_TAMs = c("MKI67", "CDK1", "CDC45"),
  RTM_TAMs = c("LYVE1", "HES1", "FOLR2"),
  pDC = c("LILRA4", "GZMB", "IL3RA", "TCF4", "TCL1A", "CLEC4C", "CLIC3"),
  cDC1 = c("CLEC9A", "FLT3", "IDO1", "IRF8", "XCR1", "BATF3", "C1orf54"),
  cDC2 = c("CD1A","CD1C", "CD1E", "FCER1A", "HLA-DQA1", "CLEC10A", "TIMP1", "CEBPD", "LST1"),
  cDC3 = c("LAMP3", "CCR7", "FSCN1", "CCL19", "CCL22", "BIRC3")
)

p <- MarkDotplot(markerlist, seu, facet = TRUE)
print(p)
ggsave("analysis/figure/04_MPS_subcluster_analysis/02_annotation/02_classical_markers_dotplot.pdf", p, width = 20, height = 4.7)

cluster_annotation <- c(
  "0" = "Maph_Inflam_CCL3",
  "1" = "Maph_LA_APOE",
  "2" = "Mono_CD14_S100A8",
  "3" = "Mono_CD16_TCF7L2",
  "4" = "cDC2_FCER1A",
  "5" = "unknown",
  "6" = "unknown",
  "7" = "cDC1_CCSER1",
  "8" = "unknown"
)

metadata$classical_marker_annotation <- cluster_annotation[as.character(metadata$MPS_subclusters)]
seu@meta.data <- metadata
p <- DimPlot(seu, reduction = "umap", group.by = "classical_marker_annotation", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/04_MPS_subcluster_analysis/02_annotation/02_classical_marker_annotation.png",p ,width = 6, height = 6, dpi = 500)


#### seu filter ####
seu <- subset(seu, subset = classical_marker_annotation == "unknown", invert = TRUE)
seu$mps_celltype <- seu$classical_marker_annotation
qs::qsave(seu, file = "analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)

p <- DimPlot(seu, reduction = "umap", group.by = "mps_celltype", pt.size = 0.05, label = T, raster = FALSE)+
  ggsci::scale_color_d3("category20")+ NoLegend()
ggsave("analysis/figure/04_MPS_subcluster_analysis/02_annotation/03_mps_celltype.png",p ,width = 6, height = 6, dpi = 500)
