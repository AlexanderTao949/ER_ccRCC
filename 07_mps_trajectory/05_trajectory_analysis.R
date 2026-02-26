#### Load packages ####
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scater)
library(scTookit)
library(ggplot2)
library(mgcv)
pseu_line <- function(seu, feature = NULL){
  data.use <- FetchData(seu, var = c("slingshot_ptime", feature, "ERScore_group"))
  data.use <- data.use %>% dplyr::arrange(slingshot_ptime)
  p <- ggplot(data.use, aes(x = slingshot_ptime, y = .data[[feature]], color = ERScore_group)) +
    geom_smooth(
      method = "gam",           
      formula = y ~ s(x, bs = "cs"),       
      linewidth = 1.2,          
      se = FALSE                
    ) +
    scale_color_manual(
      name = NULL,           # 图例标题
      values = c(
        "High_ERScore" = "#E41A1C",  # 红色
        "Low_ERScore" = "#377EB8"    # 蓝色
      ),
      breaks = c("High_ERScore", "Low_ERScore"),  # 确保图例顺序
      labels = c("High_ERScore", "Low_ERScore")   # 图例显示标签（带空格）
    ) +
    labs(
      x = "Mono_CD14_S100A8 -> Maph_Inflam_CCL3",
      y = paste(feature, "Expression Level")
    ) +
    theme_classic(base_size = 15) +
    theme(
      # legend.position = c(0.02, 0.98),  
      legend.position = "top",
      legend.justification = c("center", "top"),
      legend.background = element_blank(),
      axis.text.x = element_text(colour = "black")
    )
  return(p)
}

#### Load data ####
seu <- qs::qread("analysis/data/07_mps_trajectory/04_slingshot/02_seu_slingshot_ptime.qs", nthreads = 50)
pyscenic_auc <- qs::qread("analysis/data/06_GRN/01_pySCENIC_process/pyscenic_auc.qs")
pyscenic_auc <- pyscenic_auc[, colnames(seu)]

#### Density Plot ####
plot_data <- FetchData(seu, var = c("slingshot_ptime", "ERScore_group"))

plot_data <- plot_data %>% 
  filter(!is.na(slingshot_ptime)) %>%
  filter(!is.na(ERScore_group))

p <- ggplot(plot_data, aes(x = slingshot_ptime, color = ERScore_group, fill = ERScore_group)) +
  geom_density(alpha = 0.4, size = 1) + 
  
  scale_color_manual(values = c("High_ERScore" = "#E41A1C", 
                                "Low_ERScore" = "#377EB8")) +
  scale_fill_manual(values = c("High_ERScore" = "white", "Low_ERScore" = "white")) +
  
  theme_classic() +
  
  labs(x = "Pseudotime", 
       y = "Cell Density", 
       title = "Pseudotime Distribution") +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.8, 0.9),
    legend.title = element_blank(), 
    legend.background = element_rect(fill = "transparent")
  )

print(p)

ggsave("analysis/figure/07_mps_trajectory/05_trajectory_analysis/01_pseudotime_Density.pdf", width = 5, height = 4)

#### marker pseudotime ####
data.use <- FetchData(seu, var = c("S100A8", "CCL3", "TNF", "ERScore_group", 
                                   "slingshot_ptime", "mps_celltype"))
data.use <- data.use %>% dplyr::arrange(slingshot_ptime)

plist <- lapply(c("S100A8", "CCL3", "TNF"), function(x){p <- pseu_line(seu, feature = x)})
p <- patchwork::wrap_plots(plist, ncol = 3)+ 
  patchwork::plot_layout(guides = "collect")&
  theme(legend.position = "top")
p
ggsave("analysis/figure/07_mps_trajectory/05_trajectory_analysis/02_marker_ptime_exp.pdf", p, width = 13.68, height = 3.76)


#### TF activity pseudotime ####
seu[["SCENIC"]] <- CreateAssayObject(data = pyscenic_auc)
DefaultAssay(seu) <- "SCENIC"
data.use <- FetchData(seu, var = c("NFKB1(+)", "NFKB2(+)", "RELB(+)", # NFkB
                                   "FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)", 
                                   "CEBPB(+)", 
                                   "ERScore_group", "slingshot_ptime"))
data.use <- data.use %>% dplyr::arrange(slingshot_ptime)

plist <- lapply(c("NFKB1(+)", "NFKB2(+)", "RELB(+)", # NFkB
                    "FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)", 
                    "CEBPB(+)"), function(x){p <- pseu_line(seu, feature = x)})
p <- patchwork::wrap_plots(plist, ncol = 4)+ 
  patchwork::plot_layout(guides = "collect")&
  theme(legend.position = "top")
p
ggsave("analysis/figure/07_mps_trajectory/05_trajectory_analysis/02_marker_ptime_exp.pdf", p, width = 13.68, height = 3.76)
