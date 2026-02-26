#### Load packages ####
library(presto)
library(tidyverse)
library(GSVA)
library(bulkTookit)
library(survival)

#### load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
TCGA_expr <- qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/exprSet/vsd01A.qs")
ERscore_meta <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")

#### Find markers ####
markers <- wilcoxauc(seu,
                     group_by = "mps_celltype",
                     assay = "data",
                     seurat_assay = "RNA")
markers <- markers %>% dplyr::filter(logFC > 0)
markers$pct_diff <- markers$pct_in - markers$pct_out
qs::qsave(markers,file = "analysis/data/05_determine_key_mps_subcluster/06_ssgsea_cell_scoring/01_cell_markers.qs")

topmarkers <- markers %>% 
  dplyr::filter(feature %in% rownames(TCGA_expr)) %>% 
  group_by(group) %>% 
  dplyr::filter(auc > 0.5, padj < 0.05, pct_diff > 10) %>% 
  dplyr::arrange(desc(logFC), .by_group = TRUE) %>% 
  slice_head(n = 100) %>% 
  ungroup()

genesets <- lapply(unique(topmarkers$group), function(x){
  dat <- topmarkers %>% dplyr::filter(group == x)
  genes <- dat %>% pull(feature)
  return(genes)
})
names(genesets) <- unique(topmarkers$group)

#### ssgsea ####
TCGA_cell_scoring <- Run_ssgsea(TCGA_expr, genesets)
qs::qsave(TCGA_cell_scoring, file = "analysis/data/05_determine_key_mps_subcluster/06_ssgsea_cell_scoring/02_TCGA_cell_scoring.qs")

TCGA_dat <- cbind(ERscore_meta$TCGA_KIRC, TCGA_cell_scoring)
qs::qsave(TCGA_dat, file = "analysis/data/05_determine_key_mps_subcluster/06_ssgsea_cell_scoring/02_TCGA_cell_scoring_ERScore_meta.qs")

#### survival ####
OSdata <- Find_surv_cutoff(data = TCGA_dat, time = "OS.time", event = "OS", variable = "Maph_Inflam_CCL3")
OSdata$Maph_Inflam_CCL3 <- ifelse(OSdata$group == "high", "Maph_Inflam_CCL3 (Hgih)", "Maph_Inflam_CCL3 (Low)")

OSdata <- Find_surv_cutoff(data = OSdata, time = "OS.time", event = "OS", variable = "Everolimus_Response_Score")
OSdata$ERScore_group <- ifelse(OSdata$group == "high", "ERScore (Hgih)", "ERScore (Low)")

OSdata$group <- paste0(OSdata$Maph_Inflam_CCL3, "; ", OSdata$ERScore_group)

fit <- survfit(Surv(OS.time, OS) ~ group, data = OSdata)
groupLevels <- levels(factor(OSdata$group))
diff <- survdiff(Surv(OS.time, OS) ~ group, data = OSdata)
pVal <- 1 - pchisq(diff$chisq, df = length(groupLevels) - 1)
if (pVal < 0.001) {
  pValLabel <- "p < 0.001"
} else {
  pValLabel <- paste0("p = ", formatC(pVal, format = "f", digits = 3))
  
}

p <- ggsurvplot(
  fit,
  data = OSdata,
  conf.int = FALSE,                         # 不显示置信区间
  pval = pValLabel,                         # 显示 P 值
  pval.size = 5,
  legend.labs = groupLevels,               # 图例标签
  legend.title = "Group",                   # 图例标题
  legend = c(0.8, 0.85),                         # 图例位置
  xlab = "Time (Days)",                    # 横轴标签
  risk.table = FALSE,                       # 不显示风险表（可改 TRUE）
  title = "TCGA-KIRC (OS)"
)
p$plot <- p$plot + 
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1)
print(p)
ggsave("analysis/figure/03_determine_key_major_celltype/04_ligand_resceptor_survival_analysis/TNF_ERScore_surv_OS.pdf", width = 5.66, height = 5.66)
