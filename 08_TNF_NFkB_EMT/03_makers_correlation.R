#### load packages ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

#### Load data and process ####
ERScore_exp_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/03_ERScore_exp_list.qs")
TCGA_dat <- ERScore_exp_list[[1]]
ever_vsd <- qs::qread("RNA_seq_upstream/GSE99875/02_annotated_expression_matrices/vsd.qs")
metadata <- qs::qread("RNA_seq_upstream/GSE99875/resource/metadata.qs")

#### RNA-seq 786 ####
data.use <- ever_vsd[c("RELA","SNAI2", "ZEB1", "ZEB2", "FOXC2"), ] %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
data.use <- data.use %>% dplyr::inner_join(metadata, by = "sample") %>% column_to_rownames("sample")

pf <- function(gene){
  ggscatter(data.use, 
            x = "RELA",
            y = gene, 
            palette = c("#377EB8", "#E41A1C"),
            add = "reg.line",                    
            conf.int = TRUE, 
            color = "group"
  )+
    labs(title = NULL)+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 15),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1
    )+
    stat_cor(aes(color = group))
}
plist <- lapply(c("SNAI2", "ZEB1", "ZEB2", "FOXC2"), function(x){p <- pf(x)})
p <- patchwork::wrap_plots(plist, ncol = 2, nrow = 2)
print(p)
ggsave(filename = "figure/10_tumor_cells/04_genes_TFs_express/01_Everolimus_treatment_786/03_P65_EMT_TFs_correlation.pdf", p, height = 6.60, width = 6.60)

#### TCGA ####
data.use <- TCGA_dat[, c("Everolimus_Response_Score", "RELA","SNAI2", "ZEB1", "ZEB2", "FOXC2")]
# data.use$group <- ifelse(data.use$Everolimus_Response_Score > median(data.use$Everolimus_Response_Score), "High_ERScore", "Low_ERscore")
data.use <- data.use %>%
  mutate(group = case_when(
    Everolimus_Response_Score >= quantile(Everolimus_Response_Score, 0.95) ~ "High_ERScore",
    Everolimus_Response_Score <= quantile(Everolimus_Response_Score, 0.05) ~ "Low_ERScore",
    TRUE ~ NA_character_  # 或设为"Medium"保留中间50%样本
  ))
data.use <- data.use[!is.na(data.use$group), ]


pf <- function(gene){
  ggscatter(data.use, 
            x = "RELA",
            y = gene, 
            palette = c("#377EB8", "#E41A1C"),
            add = "reg.line",                    
            conf.int = TRUE, 
            color = "group"
  )+
    labs(title = NULL)+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 15),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1
    )+
    stat_cor(aes(color = group))
}
plist <- lapply(c("SNAI2", "ZEB1", "ZEB2", "FOXC2"), function(x){p <- pf(x)})
p <- patchwork::wrap_plots(plist, ncol = 2, nrow = 2)
print(p)
ggsave(filename = "analysis/figure/test.png", p, height = 6.60, width = 6.60, dpi = 500)
