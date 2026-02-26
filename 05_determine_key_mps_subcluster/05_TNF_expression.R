#### Load packages ####
library(Seurat)
library(tidyverse)
library(scTookit)
library(presto)
library(ggrepel)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)

#### DE analysis ####
de <- wilcoxauc(seu,
                group_by = "ERScore_group",
                groups_use = c("High_ERScore", "Low_ERScore"),
                seurat_assay = "RNA",
                assay = "data",
                verbose = TRUE)
de <- de %>% dplyr::filter(group == "High_ERScore")
de$pct_diff <- de$pct_in - de$pct_out
qs::qsave(de, file = "analysis/data/05_determine_key_mps_subcluster/05_TNF_expression/01_de_group.qs")

seu$mps_treat <- paste(seu$mps_celltype, seu$ERScore_group, sep = "_")
de_percells <- do.call(rbind, lapply(unique(seu$mps_celltype %>% as.character), function(cell){
  high_group <- paste0(cell, "_High_ERScore")
  low_group <- paste0(cell, "_Low_ERScore")
  de <- wilcoxauc(seu,
                  group_by = "mps_treat",
                  groups_use = c(high_group, low_group),
                  seurat_assay = "RNA",
                  assay = "data",
                  verbose = TRUE)
  result <- de %>% dplyr::filter(group == high_group)
  result$group <- cell
  return(result)
}))

de_percells$pct_diff <- de_percells$pct_in - de_percells$pct_out
qs::qsave(de_percells, file = "analysis/data/05_determine_key_mps_subcluster/05_TNF_expression/01_de_percell_group.qs")

#### TNF de result ####
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
ggsave("analysis/figure/05_determine_key_mps_subcluster/05_TNF_expression/01_TNF_percells_group_de.pdf", p, width = 6.67, height = 5.11)

#### TNF FeaturePlot ####
seulist <- SplitObject(seu, split.by = "ERScore_group")
TNFexp <- FetchData(seu, var = "TNF")
thres <- median(TNFexp$TNF[TNFexp$TNF > 0])
p1 <- ThreshPlot(seulist$High_ERScore, feature = "TNF", threshold = thres, dot.size = 0.8, title = "High_ERScore", legend.title = "Expression")
p2 <- ThreshPlot(seulist$Low_ERScore, feature = "TNF", threshold = thres, dot.size = 0.8, title = "Low_ERScore", legend.title = "Expression")
p <- patchwork::wrap_plots(list(p1, p2), ncol = 2)
print(p)
ggsave("analysis/figure/05_determine_key_mps_subcluster/05_TNF_expression/02_TNF_group_featureplot.pdf", p, width = 4.86, height = 2.82)

#### TNF+ cells proportion ####
TNF_meta_exp <- FetchData(seu, var = c("orig.ident", "mps_celltype", "ERScore_group", "TNF"))
TNF_meta_exp$TNF_pos_median <- ifelse(TNF_meta_exp$TNF > thres, "TNF_high", "TNF_low")

df_high_only <- TNF_meta_exp %>%
  filter(TNF_pos_median == "TNF_high") %>%
  group_by(orig.ident, ERScore_group, mps_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident, ERScore_group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p <- ggplot(df_high_only, aes(x = orig.ident, 
                              y = proportion, 
                              fill = mps_celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ERScore_group, scales = "free_x", nrow = 1) +  # 关键代码
  labs(
    x = "Sample",
    y = "Proportion of TNF-high cells",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold", size = 12)             
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0))
print(p)

ggsave("analysis/figure/05_determine_key_mps_subcluster/05_TNF_expression/03_TNF_high_cell_proportion.pdf", p, width = 8.05, height = 4.45)
