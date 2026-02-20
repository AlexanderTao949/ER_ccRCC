#### Load packgages ####
library(tidyverse)
library(scTookit)
library(Seurat)
library(presto)
library(bulkTookit)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(SCPA)
library(ggrepel)
library(patchwork)

# #### Load data and process ####
seu <- qs::qread("analysis/data/02_scRNA_data_process/03_annotation/03_seu_anno.qs", nthreads = 50)
metadata <- seu@meta.data
ERGs <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/06_ERGs.qs")
ERGs <- list(Everolimus_Response_Score = ERGs)
kegg_db <- qs::qread("~/public_data/kegg_db/kegg_human.qs")[["gene"]]
kegg_db <- kegg_db %>% dplyr::filter(level1_pathway_name == "Environmental Information Processing")
kegg_db_list <- lapply(unique(kegg_db$pathway_name), function(x){
  dat <- kegg_db[kegg_db$pathway_name == x,]
  dat <- dat %>% dplyr::select(pathway_name, gene_name)
  colnames(dat) <- c("Pathway", "Genes")
  return(dat)
})


#### pseudobulk ####
pseu_mat <- getpseudobulk(seu, group_by = "orig.ident")
vsd_pseu <- vst(as.matrix(pseu_mat)) %>% as.data.frame()
qs::qsave(pseu_mat, file = "analysis/data/02_scRNA_data_process/05_pseudobulk_groups/01_pseudobulk_mat.qs")
qs::qsave(vsd_pseu, file = "analysis/data/02_scRNA_data_process/05_pseudobulk_groups/02_pseudobulk_vsd.qs")

#### ssgsea ####
ERScore_dat <- Run_ssgsea(vsd_pseu, geneset = ERGs)
ERScore_dat$group <- ifelse(ERScore_dat$Everolimus_Response_Score > median(ERScore_dat$Everolimus_Response_Score), "High_ERScore", "Low_ERScore") %>%
  factor(levels = c("Low_ERScore", "High_ERScore"))
for (i in 1:nrow(ERScore_dat)){
  index = rownames(ERScore_dat)[i]
  metadata$ERScore[metadata$orig.ident == index] = ERScore_dat$Everolimus_Response_Score[rownames(ERScore_dat) == index]
  metadata$ERScore_group[metadata$orig.ident == index] = ERScore_dat$group[rownames(ERScore_dat) == index] %>% as.character()
}
seu@meta.data <- metadata
seu$major_celltype_treat <- paste(seu$major_celltype, seu$ERScore_group, sep = "_")
qs::qsave(ERScore_dat, file = "analysis/data/02_scRNA_data_process/05_pseudobulk_groups/03_ERScore_group.qs")
qs::qsave(seu, file = "analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthreads = 50)

ERScore_dat <- ERScore_dat %>% rownames_to_column("sample")
ERScore_dat <- ERScore_dat %>% dplyr::arrange(Everolimus_Response_Score)
ERScore_dat$sample <- factor(ERScore_dat$sample, levels = ERScore_dat$sample)
p <- ggplot(ERScore_dat, aes(x = sample,
                             y = Everolimus_Response_Score,
                             color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = median(ERScore_dat$Everolimus_Response_Score),
             linetype = "dashed",
             color = "black",
             linewidth = 0.6) +
  scale_color_manual(values = c(
    "High_ERScore" = "darkred",
    "Low_ERScore" = "darkblue"
  ), name = "ERScore_Group") +
  labs(x = "Sample",
       y = "Everolimus Response Score",
       title = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)
ggsave("analysis/figure/02_scRNA_data_process/05_pseudobulk_groups/01_ERScore_group_dotplot.pdf", p, width = 3.99, height = 3.83)

#### Cell proportion ####
stat <- table(seu$major_celltype, seu$orig.ident)
patient_prop <- as.data.frame(prop.table(stat,margin = 2))
colnames(patient_prop) <- c("major_celltype", "patient", "proportion")
patient_prop <- patient_prop %>%
  dplyr::left_join(ERScore_dat %>% dplyr::select(sample, group),
            by = c("patient" = "sample"))

celltype_colors <- c(
  "B_Cells" = "#5D8AA8",
  "CD4_T_Cells" = "#7B68EE",
  "CD8_T_Cells" = "#483D8B",  # 比CD4深，体现细胞毒性
  "Natural_Killer_Cells" = "#87CEEB",
  "Mononuclear_Phagocytes" = "#DAA520",     
  "Mast_Cells" = "#CD853F",
  "Endothelial_Cells" = "#2E8B57",
  "Fibroblasts" = "#8FBC8F",
  "Tumor_Cells" = "#8B0000"
)
p <- ggplot(patient_prop, aes(patient, proportion, fill = major_celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ group, scales = "free_x", strip.position = "top") +
  xlab(label = "") +
  ylab(label = "Cell Type Proportion") +
  scale_fill_manual(values = celltype_colors) +
  theme_classic() +
  theme(
    axis.ticks.length = unit(0.2, 'cm'),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 12, face = "plain", color = "black"),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    panel.spacing = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(title = NULL, ncol = 1))
print(p)
ggsave("analysis/figure/02_scRNA_data_process/05_pseudobulk_groups/02_celltype_proportion.pdf", p, width = 8.41, height = 4.11)

#### DESeq2 ####
dds <- DESeqDataSetFromMatrix(countData = pseu_mat,
                              colData = ERScore_dat,
                              design = ~ group)
dds <- DESeq(dds)
contrast = c("group", "High_ERScore", "Low_ERScore")
dd1 <- results(dds, contrast = contrast, alpha = 0.05)
plotMA(dd1, ylim = c(-10, 10))
dd2 <- lfcShrink(dds, contrast = contrast, res = dd1, type = "ashr")
plotMA(dd2, ylim = c(-10, 10))
deg_res <- dd2 %>%
  as.data.frame() %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  rownames_to_column("gene")
qs::qsave(deg_res, file = "analysis/data/02_scRNA_data_process/05_pseudobulk_groups/05_DESeq2_res.qs")

#### FindMarkers ####
de <- wilcoxauc(seu,
                group_by = "ERScore_group",
                groups_use = c("High_ERScore", "Low_ERScore"),
                seurat_assay = "RNA",
                assay = "data",
                verbose = TRUE)
de <- de %>% dplyr::filter(group == "High_ERScore")
qs::qsave(de, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/05_findmarkers_groups_allcells.qs")

de_percells <- do.call(rbind, lapply(unique(seu$major_celltype %>% as.character), function(cell){
  high_group <- paste0(cell, "_High_ERScore")
  low_group <- paste0(cell, "_Low_ERScore")
  de <- wilcoxauc(seu,
                  group_by = "major_celltype_treat",
                  groups_use = c(high_group, low_group),
                  seurat_assay = "RNA",
                  assay = "data",
                  verbose = TRUE)
  result <- de %>% dplyr::filter(group == high_group)
  result$group <- cell
  return(result)
}))

de_percells$pct_diff <- de_percells$pct_in - de_percells$pct_out
qs::qsave(de_percells, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/05_findmarkers_groups_percell.qs")

#### TNF expression in all celltypes ####
plot_data <- de_percells[de_percells$feature == "TNF", ]
p <- ggplot(plot_data, aes(x = pct_diff, y = logFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(color = group, size = avgExpr), alpha = 0.85, stroke = 0.5) +
  geom_text_repel(aes(label = group),
                  size = 4,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.size = 0.2,
                  max.overlaps = 15) +

  # scale_color_brewer(palette = "Set3", name = "Cell Type") +
  scale_color_discrete(name = "Cell Type")+
  scale_size_continuous(range = c(1,10), name = "Average\nExpression") +
  labs(
    x = "Expression Percentageg Difference",
    y = "Log2 Fold Change"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

print(p)
ggsave("analysis/figure/02_scRNA_data_process/05_pseudobulk_groups/03_TNF_percells_group_de.pdf", p, width = 6.67, height = 5.11)

seulist <- SplitObject(seu, split.by = "ERScore_group")
TNFexp <- FetchData(seu, var = "TNF")
thes <- mean(TNFexp$TNF) + 2*sd(TNFexp$TNF)
p1 <- ThreshPlot(seulist$High_ERScore, feature = "TNF", threshold = thes, dot.size = 0.8, title = "High_ERScore", legend.title = "Expression")
p2 <- ThreshPlot(seulist$Low_ERScore, feature = "TNF", threshold = thes, dot.size = 0.8, title = "Low_ERScore", legend.title = "Expression")
p <- patchwork::wrap_plots(list(p1, p2), ncol = 2)+
  patchwork::plot_layout(guides = 'collect')&
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8,
                                face = "bold",
                                family = "Arial",
                                color = "black"),
    legend.text = element_text(size = 6),
    legend.key.size = grid::unit(0.15, "inch"),
    legend.box = "horizontal",
    legend.box.margin = margin(10, 0, 0, 0)
  )
print(p)
ggsave("analysis/figure/02_scRNA_data_process/05_pseudobulk_groups/04_TNF_group_featureplot.png", p, width = 10, height = 6, dpi = 500)

#### SCPA KEGG ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthreads = 50)
seu_high <- seurat_extract(seu, meta1 = "major_celltype_treat", value_meta1 = "Tumor_Cells_High_ERScore")
seu_low <- seurat_extract(seu, meta1 = "major_celltype_treat", value_meta1 = "Tumor_Cells_Low_ERScore")

scpa_kegg <- compare_pathways(samples = list(seu_high, seu_low), 
                              pathways = kegg_db_list,
                              downsample = Inf,
                              min_genes = 0,
                              max_genes = 100000,
                              parallel = TRUE, 
                              cores = 50)
qs::qsave(scpa_kegg, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/06_scpa_kegg_tumor.qs")

scpa_kegg <- scpa_kegg %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'yellow3',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

p <- ggplot(scpa_kegg, aes(FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = scpa_kegg$color, stroke = 0.3) +
  xlab("Enrichment") +
  ylab("Qval") +
  xlim(-max(scpa_kegg$FC), max(scpa_kegg$FC)) +
  geom_text_repel(data=scpa_kegg[scpa_kegg$Pathway %in% c("TNF signaling pathway", "NF-kappa B signaling pathway"), ],
                  aes(label = Pathway) ,col="black", alpha = 0.8)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
print(p)