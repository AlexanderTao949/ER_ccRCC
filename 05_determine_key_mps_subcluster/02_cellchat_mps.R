#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(CellChat)
library(ggplot2)
Run_CellChat <- function(seu){
  cellchat <- createCellChat(object = seu, group.by = "cell_type_mps", assay = "RNA")
  
  CellChatDB <- CellChatDB.human 
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat) 
  options(future.globals.maxSize = 100 * 1024^3) 
  future::plan("multisession", workers = 8) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.human)
  cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

#### Load data and processing ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthreads = 50)
seu_mps <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
meta <- seu@meta.data
meta_mps <- seu_mps@meta.data

matched_indices <- match(rownames(meta), rownames(meta_mps))
meta$cell_type_mps <- meta$major_celltype
meta$cell_type_mps[!is.na(matched_indices)] <- meta_mps$mps_celltype[matched_indices[!is.na(matched_indices)]]
seu@meta.data <- meta
seu <- subset(seu, subset = cell_type_mps == "Mononuclear_Phagocytes", invert = TRUE)
seu$samples <- seu$orig.ident

seu_high <- subset(seu, subset = ERScore_group == "High_ERScore")
seu_low <- subset(seu, subset = ERScore_group == "Low_ERScore")

#### Run CellChat ####
cellchat_high <- Run_CellChat(seu_high)
cellchat_low <- Run_CellChat(seu_low)
object.list <- list(High_ERScore = cellchat_high, Low_ERScore = cellchat_low)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
qs::qsave(cellchat_high, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/02_cellchat_mps/01_cellchat_high.qs")
qs::qsave(cellchat_low, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/02_cellchat_mps/01_cellchat_low.qs")
qs::qsave(cellchat, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/02_cellchat_mps/01_cellchat_merge.qs")


cellchat <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/02_cellchat_mps/01_cellchat_merge.qs", nthreads = 25)

# all_sources <- c("B_Cells", "CD4_T_Cells", "CD8_T_Cells", 
#                  "Endothelial_Cells", "Fibroblast", "Mast_Cells", 
#                  "Mononuclear_Phagocytes", "Natural_Killer_Cells")
all_sources <- c("Maph_Inflam_CCL3", "Maph_LA_APOE", "Mono_CD14_S100A8", "Mono_CD16_TCF7L2", "cDC1_CCSER1", "cDC2_FCER1A",
                 "B_Cells", "CD4_T_Cells", "CD8_T_Cells", "Endothelial_Cells", "Fibroblast", "Mast_Cells", "Natural_Killer_Cells")

df.net <- subsetCommunication(cellchat)
df.net.high <- df.net$High_ERScore
df.net.low <- df.net$Low_ERScore

df.net.high <- df.net.high %>% 
  filter(pathway_name == "TNF", source %in% all_sources, target == "Tumor_Cells") %>% 
  arrange(desc(prob))
df.net.low <- df.net.low %>% 
  filter(pathway_name == "TNF", source %in% all_sources, target == "Tumor_Cells") %>% 
  arrange(desc(prob))


df.plot <- bind_rows(
  High_ERScore = df.net.high,
  Low_ERScore = df.net.low,
  .id = "group"
) %>%
  mutate(
    source = factor(source, levels = all_sources),
    group = factor(group, levels = c("High_ERScore", "Low_ERScore")),
    LR_pair = paste0(ligand, " - ", receptor)
  ) %>%
  complete(
    source = all_sources,
    group = c("High_ERScore", "Low_ERScore"),
    fill = list(prob = NA, pval = 1, LR_pair = "TNF - TNFRSF1A")
  )

p <- ggplot(df.plot, aes(x = source, y = LR_pair)) +
  geom_point(aes(size = prob, color = pval), 
             na.rm = TRUE, alpha = 0.9) +
  
  facet_wrap(~group, ncol = 2, scales = "free_x") +
  
  # CellChat 默认配色和尺寸
  scale_size_continuous(
    range = c(2, 8),
    name = "Interaction\nProbability",
    breaks = c(0.0001, 0.001),
    labels = expression(10^-4, 10^-3)
  ) +
  scale_color_gradient(
    low = "red", high = "grey90",
    name = "P-value",
    limits = c(0, 0.05),
    na.value = "white"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(colour = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  scale_x_discrete(drop = FALSE) 
print(p)
ggsave("manuscript_figure/Figure3/cellchat_allcells_mps_dotplot.png", width = 9.40, height = 3.64, dpi = 500)
