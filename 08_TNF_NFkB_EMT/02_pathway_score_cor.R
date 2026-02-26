#### Load packages ####
library(Seurat)
library(bulkTookit)
library(tidyverse)
library(purrr)          
library(ComplexHeatmap) 
library(circlize)       
library(grid)           

#### Load data ####
kegg_db <- qs::qread("~/public_data/kegg_db/kegg_human_pathway.qs")
ERScore_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/01_ERScore_list.qs")
hallmark <- clusterProfiler::read.gmt("analysis/resource/msigdb/h.all.v2026.1.Hs.symbols.gmt")
expr_list <- list(TCGA_KIRC = qs::qread("analysis/data/ccRCC_bulk_data/TCGA_KIRC/exprSet/vsd01A.qs"),
                  E_MTAB_1980 = qs::qread("analysis/data/ccRCC_bulk_data/E_MTAB_1980/exprSet.qs"),
                  ICGC = qs::qread("analysis/data/ccRCC_bulk_data/RECA-EU/vsd01A.qs"),
                  GSE167573 = qs::qread("analysis/data/ccRCC_bulk_data/GSE167573/exprSet.qs"))

#### pathway ssgsea ####
genesets <- list()
TNF_genes <- kegg_db %>% dplyr::filter(term == "TNF signaling pathway") %>% pull(gene)
NFkB_genes <- kegg_db %>% dplyr::filter(term == "NF-kappa B signaling pathway") %>% pull(gene)
TNF_NFkB_genes <- hallmark %>% dplyr::filter(term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>% pull(gene)
EMT_genes <- hallmark %>% dplyr::filter(term == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% pull(gene)

genesets <- list(TNF_signaling_pathway = TNF_genes,
                 NF_kappa_B_signaling_pathway = NFkB_genes,
                 HALLMARK_TNFA_SIGNALING_VIA_NFKB = TNF_NFkB_genes,
                 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = EMT_genes)

pathway_ssgsea <- lapply(names(expr_list), function(x){
  dat <- expr_list[[x]]
  res <- Run_ssgsea(dat, genesets)
  return(res)
})
names(pathway_ssgsea) <- names(expr_list)
qs::qsave(pathway_ssgsea, file = "analysis/data/08_TNF_NFkB_EMT/02_pathway_score_cor/01_pathway_ssgsea.qs")

pathway_score_ERScore <- lapply(names(ERScore_list), function(x){
  dat1 <- pathway_ssgsea[[x]]
  dat2 <- ERScore_list[[x]]
  dat <- cbind(dat1, dat2)
  return(dat)
})
names(pathway_score_ERScore) <- names(ERScore_list)

qs::qsave(pathway_score_ERScore, file = "analysis/data/08_TNF_NFkB_EMT/02_pathway_score_cor/01_pathway_ssgsea_ERScore.qs")

#### Correlation ####
cor_long <- imap_dfr(pathway_score_ERScore, function(df, cohort) {
  # 取前 4 列（通路）和最后一列（Everolimus_Response_Score）
  df_sub <- df %>% select(1:4, last_col())
  
  # 重命名最后一列，方便后续引用
  names(df_sub)[5] <- "Everolimus_Response_Score"
  
  # 转为长格式：每行是一个样本-通路组合
  df_sub %>%
    pivot_longer(
      cols = 1:4,
      names_to = "pathway",
      values_to = "score"
    ) %>%
    group_by(pathway) %>%
    summarise(
      rho  = cor.test(score, Everolimus_Response_Score, method = "spearman")$estimate,
      pval = cor.test(score, Everolimus_Response_Score, method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(cohort = cohort)   # 记录队列名称
})

# --------------------------------------------------------------
# 3. 转回宽矩阵（供 ComplexHeatmap 使用）
# --------------------------------------------------------------
# 相关系数矩阵
cor_mat <- cor_long %>%
  select(pathway, cohort, rho) %>%
  pivot_wider(names_from = cohort, values_from = rho) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# p 值矩阵（与 cor_mat 同维度）
pval_mat <- cor_long %>%
  select(pathway, cohort, pval) %>%
  pivot_wider(names_from = cohort, values_from = pval) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# --------------------------------------------------------------
# 4. 绘制 ComplexHeatmap
# --------------------------------------------------------------
# 颜色方案：-1（蓝） → 0（白） → 1（红）
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# 行注释（通路名）和列注释（队列名）
row_anno <- rowAnnotation(
  Pathway = anno_text(rownames(cor_mat), gp = gpar(fontsize = 10))
)
col_anno <- columnAnnotation(
  Cohort = anno_text(colnames(cor_mat), gp = gpar(fontsize = 10))
)

# 主体热图
ht <- Heatmap(
  cor_mat,
  name = "Spearman ρ",
  col = col_fun,
  
  # 在每个格子中写入 p 值（保留两位有效数字）
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2g", pval_mat[i, j]), 
      x, y, 
      gp = gpar(fontsize = 8, col = "black")
    )
  },
  
  # 注释
  left_annotation = row_anno,
  top_annotation   = col_anno,
  
  # 关闭默认的行列名（因为已经在注释中显示）
  show_row_names    = FALSE,
  show_column_names = FALSE,
  
  # 标题样式
  row_title    = "Pathways",
  column_title = "Cohorts",
  row_title_gp    = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 12),
  
  # 图例位置
  heatmap_legend_param = list(
    title = "Spearman ρ",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

# 画出热图
draw(ht)
