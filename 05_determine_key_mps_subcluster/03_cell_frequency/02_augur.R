#### Load packages and functions ####
library(Augur)
library(tidyverse)
library(Seurat)
library(ggplot2)

#### Load data and process ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
augur <- calculate_auc(seu,
                       cell_type_col = "mps_celltype", 
                       label_col = "ERScore_group",
                       n_threads = 50)
qs::qsave(augur, file= "~/projects/Everolimus_Resistance_ccRCC/analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/02_augur/01_augur.qs")

#### UMAP ####
auc <- augur$AUC
metadata <- seu@meta.data
metadata$augur_auc <- "unknown"
metadata$augur_auc <- auc$auc[match(metadata$mps_celltype, auc$cell_type)]
seu@meta.data <- metadata
p <- ThreshPlot(seu,feature = "augur_auc", threshold = 0, legend.position = "right")
print(p)
ggsave("analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/02_augur/01_auc_umap.png", p, width = 2.47, height = 2.65, dpi = 500)

#### bar plot ####
aucs = augur$AUC
size_sm = 6
size_lg = 7
range = range(aucs$auc)
expand = abs(diff(range)) * 0.1
cellColors <- c(
  "Maph_Inflam_CCL3" = "#E41A1C",  
  "Maph_LA_APOE"     = "#FF7F00",  
  "cDC2_FCER1A"      = "#4DAF4A",  
  "cDC1_CCSER1"      = "#377EB8",  
  "Mono_CD14_S100A8" = "#984EA3",  
  "Mono_CD16_TCF7L2" = "#41B6C4"   
)
p <- aucs %>%
  # mutate(auc = ifelse(auc < 0.5, 0.5, auc)) %>%
  ggplot(aes(x = reorder(cell_type, auc), y = auc, colour = cell_type)) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted', size = 0.3) +
  geom_point(size = 3) +
  geom_text(aes(label = format(auc, digits = 3),
                y = ifelse(auc < 0.5, 0.5, auc)), size = 4,
            nudge_y = expand,
            hjust = 0.5) +
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 1.2) +  # 增加 size 参数
  scale_y_continuous('AUC', limits = c(min(range[1] - expand, 0.5),
                                       range[2] + expand * 1.5)) +
  scale_colour_manual(values = cellColors)+
  coord_flip() +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  labs(x = NULL)

print(p)
ggsave(filename = "analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/02_augur/02_line_plot.pdf", p, width = 6.16, height = 2.08)
