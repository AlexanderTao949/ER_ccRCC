#### Load pacakges ####
library(tidyverse)
library(scTookit)
library(Seurat)
library(ggplot2)
library(gghalves)
library(ggbeeswarm)
library(ggpubr)
ViolinPlot <- function(plotdata, cellchoose){
  plotdata <- plotdata %>% dplyr::filter(major_celltype == cellchoose & TNF > 0)
  
  p <- ggplot(data = plotdata, aes(x = ERScore_group,y = TNF, color = ERScore_group)) + 
    geom_violin(trim = FALSE,position = position_dodge(0.9)) +  # size：加粗小提琴的边框,trim = FALSE：确保小提琴图的尾部不会被裁剪
    geom_signif(comparisons = list(c("High_ERScore", "Low_ERScore")), # 添加显著性，y_position控制标记的y轴位置，textsize控制文本大小，size控制线的粗细
                test = "t.test",colour = "black",
                test.args = list(var.equal = TRUE, alternative = "two.sided"),
                map_signif_level = TRUE, textsize = 6, tip_length = c(0, 0),
                y_position = max(plotdata[,"TNF"]),
                vjust=0
    )+  
    stat_summary(aes(group = ERScore_group,color = ERScore_group), position = position_dodge(0.9),
                 width = 0.2, fun = mean, geom = "errorbar", 
                 fun.min = function(x) quantile(x, 0.25),
                 fun.max = function(x) quantile(x, 0.75)) + 
    stat_summary( aes(group = ERScore_group,color = ERScore_group), position = position_dodge(0.9),
                  width = 0.2,  show.legend = FALSE,
                  fun = mean, geom = "crossbar") + 
    scale_color_manual(values = c("Low_ERScore"="#667a57","High_ERScore"="#dfadab")) + 
    scale_y_continuous(limits = c(0, max(plotdata[,"TNF"])*1.1)) + 
    labs(title = cellchoose,
         x = "",
         y = paste0("TNF expression levels")) + 
    theme_classic() + 
    theme (legend.position = "none",
           axis.line = element_line(color = "black"), 
           axis.title = element_text(size = 10,,face = "bold"),
           axis.text = element_text(size = 9),
           plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  return(p)
}

#### Load data ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthreads = 50)
exp_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/03_ERScore_exp_list.qs")
meta_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")

#### TNF expression in Macrophages and Monocytes ####
plotdata <- FetchData(seu, var = c("major_celltype", "ERScore_group", "TNF"))
p <- ViolinPlot(plotdata, "Mononuclear_Phagocytes")
print(p)
ggsave("analysis/figure/03_determine_key_major_celltype/03_ligand_receptor_expression/01_LR_exp_violin_plot.pdf", p, width = 2.64, height = 3.20)

#### Correlation between TNF expression and markers in scRNA-seq ####
plotdata <- FetchData(seu, var = c("major_celltype", "TNF", "CD68", "ITGAM", "FCGR1A", "CD14", "LYZ", "ITGB2"))
data.use <- plotdata %>% dplyr::filter(major_celltype == "Mononuclear_Phagocytes" & TNF > 0 & CD68 > 0 & ITGAM > 0 & FCGR1A > 0)
plist <- lapply(c("CD68", "ITGAM", "FCGR1A"), function(x){
  p <- ScatterPlot(data.use, x_var = x, y_var = "TNF", Xlab = paste("Expression Level of", x), Ylab = "Expression Level of TNF", title = NULL)
})
p <- patchwork::wrap_plots(plist, ncol = 3)
print(p)
ggsave("analysis/figure/03_determine_key_major_celltype/03_ligand_receptor_expression/02_ligand_markers_scRNA_cor.png", p, width = 15, height = 5, dpi = 500)

#### Correlation between TNF/TNFRSF1A expression and ERScore in scRNA-seq ####
plist <- lapply(names(exp_list), function(x){
  data.use <- exp_list[[x]]
  p1 <- ScatterPlot(data.use, x_var = "TNF", y_var = "Everolimus_Response_Score", 
                   Xlab = "TNF", Ylab = "ERScore", title = x)
  p2 <- ScatterPlot(data.use, x_var = "TNFRSF1A", y_var = "Everolimus_Response_Score", 
                    Xlab = "TNFRSF1A", Ylab = "ERScore", title = x)
  p <- patchwork::wrap_plots(list(p1, p2), ncol = 2)
})
p <- patchwork::wrap_plots(plist, ncol = 1)
print(p)                
ggsave("analysis/figure/03_determine_key_major_celltype/03_ligand_receptor_expression/03_LR_ERScore_cor_cohorts.pdf",p, width = 6.27, height = 10)


#### Correlation between TNF and Mononuclear_Phagocytes markers (cohorts) ####
plist <- lapply(names(exp_list), function(x){
  data.use <- exp_list[[x]]
  p_list <- list()
  for(i in c("CD68", "ITGAM", "FCGR1A", # Maph
             "CD14", "LYZ", "ITGB2",    # Mono
             "CLEC9A", "XCR1", "BATF3", # cDC1
             "CD1C", "FCER1A", "CLEC10A")){ # cDC2
    p <- ScatterPlot(data.use, x_var = "TNF", y_var = i, 
                     Xlab = "TNF", Ylab = i, title = x)
    p_list[[i]] <- p
  }
  pp <- patchwork::wrap_plots(p_list, ncol = 1)
  return(pp)
})
