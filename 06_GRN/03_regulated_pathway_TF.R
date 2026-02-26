#### Load packages and functions ####
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
source("analysis/resource/R/regulon_specificity.R")
source("analysis/resource/R/pyscenic_IO.R")

#### Load data and process ####
seu <- qs::qread("analysis/data/06_GRN/01_pySCENIC_process/mps_seu_pyscenic.qs", nthreads = 50)
de <- qs::qread("analysis/data/05_determine_key_mps_subcluster/05_TNF_expression/01_de_percell_group.qs")
de <- de %>% dplyr::filter(group == "Maph_Inflam_CCL3")
regulon <- read.gmt("analysis/data/06_GRN/01_pySCENIC_process/pyscenic_TF.regulons.gmt")

importancedata <- LoadpySCENICOutput(regulon.gmt = "analysis/data/06_GRN/01_pySCENIC_process/pyscenic_TF.regulons.gmt",
                           adj.mat.file = "analysis/data/06_GRN/01_pySCENIC_process/pyscenic_step1_adj.tsv")
importancedata <- subset(importancedata, importance > 1)

#### KEGG ####
gene_mapping <- bitr(de$feature,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
colnames(gene_mapping) <- c("feature", "entrez")
gene_df <- de %>% dplyr::inner_join(gene_mapping, by = "feature")

geneList <- gene_df$logFC
names(geneList) <- gene_df$entrez
geneList <- sort(geneList, decreasing = T)

y <- gseKEGG(
  geneList,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)
y <- setReadable(y, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
y@result <- subset(y@result, NES > 0 & p.adjust < 0.05)
yd <- as.data.frame(y)
qs::qsave(y, file = "analysis/data/06_GRN/03_regulated_pathway_TF/01_maph_CCL3_kegg_gsea.qs")

yd$Counts <- lengths(strsplit(yd$core_enrichment, "/"))
yd$Description <- factor(yd$Description, levels = rev(yd$Description))
p <- ggplot(yd, aes(NES, Description)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Counts)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  ##scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL)+
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 0.5,size = 13), 
        axis.ticks.y = element_blank(), ## 删去y轴刻度线
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1,1,1,1), "cm"),#画布边缘距离上(top)、右(right)、下(bottom)、左(left) 
        plot.title = element_text(hjust = 0.5,size =  15),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13), 
        legend.position = "right",
        legend.background = element_rect(fill = 'transparent'))
print(p)
ggsave("analysis/figure/06_GRN/03_regulated_pathway_TF/01_maph_CCL3_kegg_gsea.pdf", p, width = 9.06, height = 4.13)

#### TF enrichment ####
pathways <- "hsa04657"
core_enrich_genes <- lapply(pathways, function(x) {
  unlist(strsplit(yd[x, "core_enrichment"], split = "/"))
}) %>%
  unlist() %>%
  unique()

core_enrich_genes <- strsplit(yd["hsa04657", "core_enrichment"], split = "/") %>% unlist()

regulon_enrich <- enricher(gene = core_enrich_genes,
                           pvalueCutoff = 10,
                           TERM2GENE = regulon,
                           minGSSize = 10,
                           maxGSSize = 3200)
regulon_enrich_df <- as.data.frame(regulon_enrich) %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::arrange(desc(FoldEnrichment))
regulon_enrich_df$Description <- factor(regulon_enrich_df$Description, levels = rev(regulon_enrich_df$Description))

# geneList <- de$logFC
# names(geneList) <- de$feature
# geneList <- sort(geneList, decreasing = T)
# 
# regulon_enrich <- GSEA(
#   geneList,
#   exponent = 1,
#   minGSSize = 10,
#   maxGSSize = 3200,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   TERM2GENE = regulon,
#   verbose = TRUE
# )
# regulon_enrich_df <- as.data.frame(regulon_enrich) %>% dplyr::filter(p.adjust < 0.05, NES > 0) %>% dplyr::arrange(desc(NES))
# regulon_enrich_df$Description <- factor(regulon_enrich_df$Description, levels = rev(regulon_enrich_df$Description))

TF_mark <- list(NFκB_signaling_pathway = c("NFKB1(+)", "NFKB2(+)", "RELB(+)"),
                MAPK_signaling_pathway = c("FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)"))

regulon_enrich_df <- regulon_enrich_df %>%
  mutate(Pathway = case_when(
    Description %in% TF_mark$NFκB_signaling_pathway ~ "NFκB signaling",
    Description %in% TF_mark$MAPK_signaling_pathway ~ "MAPK signaling",
    TRUE ~ "Other TFs"
  ))

# 2. 修改 ggplot 代码
p <- ggplot(data = regulon_enrich_df, 
            aes(x = reorder(Description, FoldEnrichment),  # 建议按值排序
                y = FoldEnrichment, 
                fill = Pathway)) +  # 将 fill 映射到通路分组
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 0.2) + 
  # 3. 手动设置颜色：NFκB用红色，MAPK用绿色，其他保持原蓝色
  scale_fill_manual(values = c("NFκB signaling" = "#e41a1c", 
                               "MAPK signaling" = "#4daf4a",
                               "Other TFs" = "black")) +
  labs(x = NULL, y = "NES", fill = "Pathway") +  # 添加图例标题
  coord_flip() + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 13), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13), 
        legend.position = "right",
        legend.background = element_rect(fill = 'transparent'))

print(p)
qs::qsave(regulon_enrich, "analysis/data/06_GRN/03_regulated_pathway_TF/02_regulon_enrich_res.qs")
ggsave(filename = "analysis/figure/06_GRN/03_regulated_pathway_TF/02_TFs_enrich.pdf", p, height = 5.17, width = 8.18)


#### TFs expression between groups ####
DefaultAssay(seu) <- "SCENIC"
data.use <- FetchData(seu, var = c("NFKB1(+)", "NFKB2(+)", "RELB(+)", # NFkB
                                   "FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)", 
                                   "CEBPB(+)", 
                                   "ERScore_group"))

##### **Plot** #####
plist <- list()
for(x in c("NFKB1(+)", "NFKB2(+)", "RELB(+)", # NFkB
           "FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)")){
  cat("\n", x,"\n")
  plotdata <- data.frame(
    gene_expr = data.use[[x]],  
    ERScore_group = data.use[["ERScore_group"]]
  )
  
  # 修改因子水平：将标签从 High_ERScore/Low_ERScore 改为 High/Low
  plotdata$ERScore_group <- factor(plotdata$ERScore_group, 
                                   levels = c("High_ERScore", "Low_ERScore"),
                                   labels = c("High", "Low"))
  
  g <- gsub("\\(\\+\\)", "", x)  
  p <- ggplot(data = plotdata, aes(x = ERScore_group, y = gene_expr, color = ERScore_group)) + 
    geom_violin(trim = FALSE, position = position_dodge(0.9)) +  
    geom_signif(comparisons = list(c("High", "Low")),  # 使用新的标签名
                test = "wilcox.test", colour = "black",
                test.args = list(var.equal = TRUE, alternative = "two.sided"),
                map_signif_level = TRUE, textsize = 6, tip_length = c(0, 0),
                y_position = max(plotdata[,"gene_expr"]),
                vjust = 0
    )+  
    stat_summary(aes(group = ERScore_group, color = ERScore_group), 
                 position = position_dodge(0.9),
                 width = 0.2, fun = mean, geom = "errorbar", 
                 fun.min = function(x) quantile(x, 0.25),
                 fun.max = function(x) quantile(x, 0.75),
                 show.legend = T) +  # 避免显示errorbar的图例
    stat_summary(aes(group = ERScore_group, color = ERScore_group), 
                 position = position_dodge(0.9),
                 width = 0.2, fun = mean, geom = "crossbar",
                 show.legend = F) +  # 避免显示crossbar的图例
    scale_color_manual(name = "ERScore",  # 图例标题
                       values = c("High" = "#E41A1C",  # 红色对应High
                                  "Low" = "#377EB8"),   # 蓝色对应Low
                       labels = c("High", "Low")) +     # 图例标签
    scale_y_continuous(limits = c(0, max(plotdata[,"gene_expr"])*1.1)) + 
    labs(title = NULL,
         x = "",
         y = paste("Transcription factor", g, "activity")) + 
    theme_classic() + 
    theme(legend.position = "right",  # 显示图例（可改为"top"/"bottom"/"left"）
          axis.line = element_line(color = "black"), 
          axis.title = element_text(size = 10, face = "bold"),  # 修正双逗号错误
          axis.text = element_text(size = 9),
          plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  plist[[g]] <- p
}


p <- patchwork::wrap_plots(plist, ncol = 4)+
  patchwork::plot_layout(guides = "collect")&
  theme(legend.position = "top")
print(p)
ggsave(filename = "analysis/figure/06_GRN/03_regulated_pathway_TF/03_TFs_activity_violin.pdf", p, width = 10.65, height = 8.85)


#### UMAP ####
DefaultAssay(seu) <- "SCENIC"
p1 <- DimPlot2(seu, reduction = "UMAP", group.by = "mps_celltype", group.highlight = "Maph_Inflam_CCL3")
ggsave(filename = "analysis/figure/06_GRN/03_regulated_pathway_TF/04_maph_CCL3L1_umap.pdf", p1, height = 4.30, width = 4.30)

seu_high <- subset(seu, subset = ERScore_group == "High_ERScore")
plist <- lapply(c("NFKB1(+)", "NFKB2(+)", "RELB(+)", # NFkB
                  "FOS(+)", "FOSB(+)", "FOSL2(+)", "JUNB(+)", "JUND(+)"), function(x){p <- DimPlot2(seu_high, reduction = "UMAP", regulon = x)})
p <- patchwork::wrap_plots(plist, ncol = 3)
print(p)
ggsave(filename = "figure/13_pyscenic_mps/03_maph_CCL3L1_TFs_umap.pdf", p, height = 8.6, width = 12.9)
