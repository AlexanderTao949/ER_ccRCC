#### Load packages ####
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scater)
library(scTookit)
library(ggplot2)
library(mgcv)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs",nthreads = 50)
seu <- subset(seu, subset = mps_celltype %in% c("Maph_Inflam_CCL3", "Mono_CD14_S100A8"))

#### slingshot
scale.data <- GetAssayData(seu, slot = "scale.data", assay = "RNA")
scale.gene <- rownames(seu)
counts <- GetAssayData(seu, slot = "counts", assay = "RNA")
counts <- counts[scale.gene, ]
sim <- SingleCellExperiment(assays = List(counts = counts))#输入为counts格式
umap <- seu@reductions$umap@cell.embeddings
colnames(umap) <- c('UMAP-1', 'UMAP-2')
reducedDims(sim) <- SimpleList(UMAP = umap)
colData(sim) <- DataFrame(seu@meta.data)
sim <- slingshot(sim, 
                 clusterLabels = 'mps_celltype',  
                 reducedDim = 'UMAP',  
                 start.clus = "Mono_CD14_S100A8",  
                 end.clus = "Maph_Inflam_CCL3")  
qs::qsave(sim, file = "analysis/data/07_mps_trajectory/04_slingshot/01_slingshot_res.qs")


#### All lineages plot ####
celltype <- colData(sim)$mps_celltype 
plotcol <- brewer.pal(n = length(unique(celltype)), name = "Set1")[as.factor(celltype)]
pdf("analysis/figure/07_mps_trajectory/04_slingshot/01_celltype_lineages.pdf", width = 6.5, height = 6.5)
plot(reducedDims(sim)$UMAP, 
     col = plotcol, 
     pch = 16, 
     asp = 1,
     main = "Cell Type")
lines(SlingshotDataSet(sim), lwd=2, col="black")
legend("bottomright",
       legend = levels(as.factor(celltype)),
       col = brewer.pal(n = length(unique(celltype)), name = "Set1"),
       pch = 16,
       cex = 0.7,     # 文字大小
       pt.cex = 0.7)  # 点的大小
dev.off()

ptime <- colData(sim)$slingPseudotime_1
seu$slingshot_ptime <- ptime
metadata <- seu@meta.data
qs::qsave(seu, file = "analysis/data/07_mps_trajectory/04_slingshot/02_seu_slingshot_ptime.qs", nthreads = 50)

#### ptime umap ####
seu <- qs::qread("analysis/data/07_mps_trajectory/04_slingshot/02_seu_slingshot_ptime.qs", nthreads = 50)
p <- ThreshPlot(seu, 
                feature = "slingshot_ptime",
                threshold = 0,
                dot.size = 1,
                alpha = 1,
                reduction = "umap",
                color = RColorBrewer::brewer.pal(9, "Blues"),
                legend.position = "right",
                legend.title = "pseudotime")
print(p)
ggsave("analysis/figure/07_mps_trajectory/04_slingshot/02_ptime.pdf", p, width = 3.69, height = 2.51)












