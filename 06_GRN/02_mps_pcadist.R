#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(RColorBrewer)

#### Load data and process ####
seu <- qs::qread("analysis/data/06_GRN/01_pySCENIC_process/mps_seu_pyscenic.qs", nthreads = 50)

DefaultAssay(seu) <- "SCENIC"
seu <- ScaleData(seu)
seu <- RunPCA(object = seu,
              features = rownames(seu),
              reduction.name = "pcaRAS",
              reduction.key = "pcaRAS_")

####
groups <- c("High_ERScore", "Low_ERScore")
celltypes <- unique(seu$mps_celltype)
sample.size <- 100
npcs <- 20
N = 1000

res <- lapply(celltypes, function(ct) {
  dists <- sapply(1:N, function(i) {
    cells1 <- rownames(subset(seu@meta.data, mps_celltype == ct & ERScore_group == groups[1]))
    cells2 <- rownames(subset(seu@meta.data, mps_celltype == ct & ERScore_group == groups[2]))
    set.seed(i)
    cells1 <- sample(cells1, size = sample.size)
    cells2 <- sample(cells2, size = sample.size)
    emb1 <- Embeddings(seu, reduction = "pcaRAS")[cells1, 1:npcs]
    emb2 <- Embeddings(seu, reduction = "pcaRAS")[cells2, 1:npcs]
    dists <- proxy::dist(emb1, emb2)
    median(as.numeric(dists))
  })
  data.frame(
    celltype = ct,
    dist = dists
  )
}) %>% do.call(rbind, .)

ct.levels <- res %>% group_by(celltype) %>% 
  reframe(mean.dist = median(dist)) %>% 
  arrange(desc(mean.dist)) %>% 
  dplyr::pull(celltype)

mps_colors <- c(
  "Mono_CD14_S100A8" = "#4DBBD5", 
  "Mono_CD16_TCF7L2" = "#00798C", 
  "Maph_Inflam_CCL3" = "#3CB371", 
  "Maph_LA_APOE"     = "#195F43", 
  "cDC2_FCER1A"      = "#A2C865", 
  "cDC1_CCSER1"      = "#6E8B3D"  
)

p <- ggplot(res, aes(factor(celltype, levels = ct.levels), dist, fill = celltype)) + 
  geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, outlier.alpha = 0.6) +
  scale_fill_brewer(palette = "Set3", guide = "none") +  # 使用调色板
  labs(y = "Distance",
       x = NULL) +
  scale_fill_manual(values = mps_colors, guide = "none") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 去掉x轴主网格线
    panel.grid.minor = element_blank(),    # 去掉次网格线
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # x轴标签倾斜
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )
print(p)
ggsave(filename = "figure/13_pyscenic_mps/03_celltype_TFs_pcadist.pdf", p, height = 3.66, width = 5.41)
