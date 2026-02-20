#### Load packages and functions ####
library(Seurat)
library(tidyverse)
library(CellChat)
library(ggplot2)
Run_CellChat <- function(seu){
  cellchat <- createCellChat(object = seu, group.by = "major_celltype", assay = "RNA")
  
  CellChatDB <- CellChatDB.human 
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat) 
  options(future.globals.maxSize = 100 * 1024^3) 
  future::plan("multisession", workers = 10) 
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
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process_1/05_pseudobulk_groups/04_seu_group.qs", nthread = 50)
seu_high <- subset(seu, subset = ERScore_group == "High_ERScore")
seu_low <- subset(seu, subset = ERScore_group == "Low_ERScore")

#### Run CellChat ####
cellchat_high <- Run_CellChat(seu_high)
cellchat_low <- Run_CellChat(seu_low)
object.list <- list(High_ERScore = cellchat_high, Low_ERScore = cellchat_low)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
qs::qsave(cellchat_high, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/03_determine_key_major_celltype/02_CellChat_major_celltype/01_cellchat_high.qs")
qs::qsave(cellchat_low, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/03_determine_key_major_celltype/02_CellChat_major_celltype/01_cellchat_low.qs")
qs::qsave(cellchat, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/03_determine_key_major_celltype/02_CellChat_major_celltype/01_cellchat_merge.qs")


# cellchat_high <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/data/05_cellchat_all/01_cellchat_high.qs", nthreads = 25)
# 
# p <- netVisual_bubble(cellchat_high, 
#                       sources.use = 1:9,            
#                       targets.use = 10, 
#                       remove.isolate = FALSE,      
#                       signaling = "TNF",    
#                       grid.on = T,                 
#                       color.grid = "black")
# 
# range(p$data$pval, na.rm = T)  # 查看通信概率(prob)的范围
# summary(p$data$prob, na.rm = T) # 查看概率的统计摘要
# 
# p1 <- p + 
#   scale_size_continuous(range = c(3, 9),
#                         breaks = c(1, 3),
#                         labels = c("p > 0.05", "p < 0.01"),
#                         name = "P-Value")+  
#   scale_color_gradientn(
#     colours = c("#2760a9", "white", "#e50f20"),
#     values = scales::rescale(c(0.1058550, 0.1307231, 0.1555912)),  
#     breaks = c(0.1058550, 0.1555912),  
#     labels = c("min", "max"),     
#     name = "Commun. Prob."
#   ) +
#   coord_fixed(ratio = 1) +
#   xlab(NULL) + 
#   ylab(NULL) +
#   theme(
#     axis.text.x = element_text(
#       size = 14, 
#       angle = 90,
#       hjust = 1
#     ),
#     axis.text.y = element_text(
#       size = 14, 
#       face = "italic"
#     ),
#     panel.border = element_rect(
#       color = "black", 
#       fill = NA, 
#       size = 1
#     ),
#     legend.key.size = unit(1, 'cm'),
#     legend.spacing.x = unit(0.5, "cm"),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(
#       size = 13, 
#       face = "bold",
#       margin = margin(b = 10)
#     ),
#     legend.box = "horizontal",
#     legend.direction = "vertical",
#     legend.position = "right"
#   )
# 
# print(p1)
# ggsave(filename = "figure/05_CellChat_all_cells/01_dotplot.pdf", p, width = 10.51, height = 6.43)
# 
# ggsave(filename = "figure/05_CellChat_all_cells/01_dotplot.pdf", p, height = 6.69, width = 6.83)
