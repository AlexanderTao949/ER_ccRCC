#### Load packages ####
library(bulkTookit)
library(Seurat)
library(cellGeometry)
library(tidyverse)
library(survival)
library(survminer)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
TCGA_counts <- qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/exprSet/counts01A.qs")
ERScore_meta_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")
TCGA_ERScore_meta <- ERScore_meta_list[["TCGA_KIRC"]]
mat <- seu@assays$RNA@counts
meta <- seu@meta.data
subcl <- meta$mps_celltype

#### deconvolution ####
mk <- cellMarkers(mat, subclass = subcl, bulkdata = TCGA_counts, cores = 10)
fit <- deconvolute(mk, TCGA_counts, use_filter = FALSE)
percentage_matrix <- fit$subclass$percent %>% as.data.frame()
qs::qsave(fit, "analysis/data/05_determine_key_mps_subcluster/06_deconvolution/01_TCGA_deconvolustion_res.qs")
qs::qsave(percentage_matrix, "analysis/data/05_determine_key_mps_subcluster/06_deconvolution/01_TCGA_percentage.qs")

TCGA_dat <- cbind(TCGA_ERScore_meta, percentage_matrix)

#### survival ####






OSdata <- Find_surv_cutoff(data = TCGA_dat, time = "OS.time", event = "OS", variable = "cDC1_CCSER1")
OSdata$cDC1_CCSER1_group <- ifelse(OSdata$group == "high", "cDC1_CCSER1 (Hgih)", "cDC1_CCSER1 (Low)")

OSdata <- Find_surv_cutoff(data = OSdata, time = "OS.time", event = "OS", variable = "Everolimus_Response_Score")
OSdata$ERScore_group <- ifelse(OSdata$group == "high", "ERScore (Hgih)", "ERScore (Low)")

OSdata$group <- paste0(OSdata$cDC1_CCSER1_group, "; ", OSdata$ERScore_group)

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
