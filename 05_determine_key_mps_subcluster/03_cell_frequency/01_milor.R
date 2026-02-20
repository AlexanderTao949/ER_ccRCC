#### Load pacakges ####
library(miloR)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(ggplot2)
library(ggbeeswarm)

#### Load data and process ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
meta.data <- seu@meta.data[,c("orig.ident", "ERScore_group", "mps_celltype")]
colnames(meta.data) = c("sample", "group","celltype")
pca.embeddings <- Embeddings(seu, reduction = "pca")
harmony.embeddings <- Embeddings(seu, reduction = "harmony")
umap.embeddings <- Embeddings(seu, reduction = "umap")

#### Run milor ####
sce <- SingleCellExperiment(assays=list(counts=seu[["RNA"]]@counts),
                            reducedDims=SimpleList(PCA=pca.embeddings,
                                                   HARMONY=harmony.embeddings,
                                                   UMAP=umap.embeddings))
milo <- Milo(sce)
colData(milo) <- DataFrame(meta.data)
reducedDimNames(milo)
milo <- buildGraph(milo, k = 30, d = 20, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.05, k = 30, d=20, refined = TRUE, reduced_dims = "HARMONY")
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sample")

milo_design <- as.data.frame(colData(milo))
milo_design <- distinct(milo_design[,-3])
milo_design$group = factor(milo_design$group,levels = c("Low_ERScore","High_ERScore"))
rownames(milo_design) <- milo_design$sample
milo_design
milo <- calcNhoodDistance(milo, d=20, reduced.dim = "HARMONY")
da_results <- testNhoods(milo, design = ~ group, design.df = milo_design, reduced.dim = "HARMONY")

milo <- buildNhoodGraph(milo)
qs::qsave(milo, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/01_milo.qs")
## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "UMAP", colour_by="group", text_by = "celltype",
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=1)+
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 8),  # 标签文本大小
    legend.title = element_text(size = 10, face = "bold")  # 标题文本大小及加粗
  )

umap_plot <- umap_pl + nh_graph_pl + 
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(guide = guide_legend(ncol = 2))
print(umap_plot)
ggsave(filename = "analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/01_umap.pdf", umap_plot, height = 7.56, width = 7.65)

da_results <- annotateNhoods(milo, da_results, coldata_col = "celltype")
qs::qsave(da_results, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/02_da_results.qs")

da_results <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/02_da_results.qs")
p <- ggplot(da_results, aes(y = celltype, x = logFC, color = logFC)) +
  geom_quasirandom(alpha = 1, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_color_gradient2(high = "#D8383A", mid = "#7d7f82", low = "dodgerblue4", midpoint = 0) +
  labs(x = "log(fold change)", y = NULL)+
  theme_classic()+
  scale_x_continuous(limits = c(-max(abs(da_results$logFC)), max(abs(da_results$logFC))))+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )
print(p)
ggsave(filename = "analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/02_bee_plot.pdf", p, width = 6.10, height = 2.63)
