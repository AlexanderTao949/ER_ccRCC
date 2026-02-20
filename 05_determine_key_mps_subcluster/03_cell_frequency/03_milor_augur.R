#### Load packages and functions ####
library(dplyr)
library(ggplot2)

#### Load data and process ####
milo.res <- qs::qread("analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/01_milor/02_da_results.qs")
augur.res <- qs::qread("analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/02_augur/01_augur.qs")

data1 <- milo.res %>% group_by(celltype) %>% summarise(logFC_mean = mean(logFC),xerr=sd(logFC))
colnames(data1) <- c("cell_type", "logFC", "xerr")

data2 <- subset(augur.res$results, metric == "roc_auc")
data2 <- data2 %>% group_by(cell_type) %>% summarise(AUC = mean(estimate),yerr=sd(estimate))
colnames(data2) <- c("cell_type", "auc", "yerr")

data_plot <- inner_join(data1, data2, by = "cell_type")
cellColors <- c(
  "Maph_Inflam_CCL3" = "#E41A1C",  
  "Maph_LA_APOE"     = "#FF7F00",  
  "cDC2_FCER1A"      = "#4DAF4A",  
  "cDC1_CCSER1"      = "#377EB8",  
  "Mono_CD14_S100A8" = "#984EA3",  
  "Mono_CD16_TCF7L2" = "#41B6C4"   
)
p <- ggplot(data_plot, aes(x = logFC, y = auc, color = cell_type, label = cell_type)) +
  geom_point(size = 4) +  
  # geom_errorbar(aes(ymin = auc - yerr, ymax = auc + yerr), width = 0.1,alpha=0.3) +  
  # geom_errorbarh(aes(xmin = logFC- xerr, xmax = logFC + xerr), height = 0.1,alpha=0.3) +  
  geom_hline(yintercept = 0.5, linetype="dashed", color='black') + 
  geom_vline(xintercept = 0, linetype="dashed", color='black') + 
  ggrepel::geom_text_repel(show.legend = F,color="black") + 
  scale_colour_manual(values = cellColors)+
  theme_bw(base_size = 15)+
  coord_flip() +
  labs(x = "logFC by miloR",y = "AUC score by Augur") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        aspect.ratio = 1)
print(p)
ggsave(file = "analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/03_milor_augur/03_milor_augur.pdf", p, width = 4.71, height = 4.26)
