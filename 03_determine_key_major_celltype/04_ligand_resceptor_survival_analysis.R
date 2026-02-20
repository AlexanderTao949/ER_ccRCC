#### Load packages ####
library(Seurat)
library(tidyverse)
library(bulkTookit)
library(survival)

#### Load data ####
exp_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/03_ERScore_exp_list.qs")
meta_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/02_ERScore_meta_list.qs")
TCGA_dat <- cbind(meta_list[[1]], exp_list[[1]][,2:ncol(exp_list[[1]])])

#### TNF and ERScore ####
OSdata <- Find_surv_cutoff(data = TCGA_dat, time = "OS.time", event = "OS", variable = "TNF")
OSdata$TNF_group <- ifelse(OSdata$group == "high", "TNF (Hgih)", "TNF (Low)")

OSdata <- Find_surv_cutoff(data = OSdata, time = "OS.time", event = "OS", variable = "Everolimus_Response_Score")
OSdata$ERScore_group <- ifelse(OSdata$group == "high", "ERScore (Hgih)", "ERScore (Low)")

OSdata$group <- paste0(OSdata$TNF_group, "; ", OSdata$ERScore_group)

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



#### TNF and TNFRSF1A ####
OSdata <- Find_surv_cutoff(data = TCGA_dat, time = "OS.time", event = "OS", variable = "TNF")
OSdata$TNF_group <- ifelse(OSdata$group == "high", "TNF (Hgih)", "TNF (Low)")

OSdata <- Find_surv_cutoff(data = OSdata, time = "OS.time", event = "OS", variable = "TNFRSF1A")
OSdata$TNFRSF1A_group <- ifelse(OSdata$group == "high", "TNFRSF1A (Hgih)", "TNFRSF1A (Low)")
OSdata$group <- paste0(OSdata$TNF_group, "; ", OSdata$TNFRSF1A_group)

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
ggsave("analysis/figure/03_determine_key_major_celltype/04_ligand_resceptor_survival_analysis/TNF_TNFRSF1A_surv_OS.pdf", width = 5.66, height = 5.66)