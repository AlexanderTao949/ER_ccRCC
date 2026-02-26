#### Load data and Process ####
library(tidyverse)
library(Seurat)
library(scTookit)
library(bulkTookit)
library(clusterProfiler)
library(GSVA)
library(GseaVis)
library(ggplot2)
library(RcisTarget)

#### Load data and process ####
ever_degs <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/01_DESeq2_res.qs")
ever_vsd <- qs::qread("RNA_seq_upstream/GSE99875/02_annotated_expression_matrices/vsd.qs")
ever_kegg_gsea <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/03_gsea_kegg.qs")
ever_kegg_gsea_dat <- as.data.frame(ever_kegg_gsea) %>% arrange(desc(NES))
ERscore_expr_list <- qs::qread("analysis/data/01_Everolimus_signature_scoring/02_cohorts_ERScore/01_ssgsea/03_ERScore_exp_list.qs")
kegg_db <- qs::qread("~/public_data/kegg_db/kegg_human_pathway.qs")
hallmark <- clusterProfiler::read.gmt("analysis/resource/msigdb/h.all.v2026.1.Hs.symbols.gmt")
motifRankings <- importRankings("analysis/resource/cisTargetDB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
data(motifAnnotations_hgnc)

#### KEGG ####
##### GSEA #####
ever_kegg_gsea_dat$Counts <- lengths(strsplit(ever_kegg_gsea_dat$core_enrichment, "/"))
ever_kegg_gsea_dat$Description <- factor(ever_kegg_gsea_dat$Description, levels = rev(ever_kegg_gsea_dat$Description))
p <- ggplot(ever_kegg_gsea_dat[1:10, ], aes(NES, Description)) + 
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
ggsave("analysis/figure/08_TNF_NFkB_EMT/01_enrichment/01_kegg_gsea_all.pdf", p, width = 10.32, height = 4.52)

TNFgenes <- c("TNFRSF1A", "TNFRSF1B", "TRAF2", "TRAF5", "TRADD", "RIPK1", "BIRC2", "BIRC3", "TNF",
              "NFKB1", "RELA", "NFKBIA", "IKBKB", "TNFAIP3", "IKBKG",
              "TAB2", "TAB3", "MAP3K7")

p <- gseaNb(object = ever_kegg_gsea,
            geneSetID = 'hsa04668',
            subPlot = 2,
            geneCol = "black",
            geneSize = 4,
            addPval = TRUE,
            addGene = TNFgenes,
            arrowType = 'open',
            kegg = T)
print(p)
ggsave(filename = "analysis/figure/08_TNF_NFkB_EMT/01_enrichment/01_gsea_kegg_TNF.pdf", p, height = 4.94, width = 5.64)

NFKBgenes <- c("TNFRSF1A", 
               "RIPK1", "TRADD", "TRAF2", "TRAF5",
               "BIRC2", "BIRC3", 
               "MAP3K7", "TAB1", "TAB2", "TAB3",
               "IKBKG", "CHUK", "IKBKB",
               "NFKBIA", "NFKB1", "RELA")
p <- gseaNb(object = ever_kegg_gsea,
            geneSetID = 'hsa04064',
            subPlot = 2,
            geneCol = "black",
            geneSize = 4,
            addPval = TRUE,
            addGene = NFKBgenes,
            arrowType = 'open',
            kegg = T)
print(p)
ggsave(filename = "analysis/figure/08_TNF_NFkB_EMT/01_enrichment/01_gsea_kegg_NFKB.pdf", p, height = 4.94, width = 5.64)


#### Hall ####
##### GSEA #####
gene_list <- ever_degs$log2FoldChange
names(gene_list) <- ever_degs$gene
gene_list <- sort(gene_list, decreasing = TRUE)

y_hall <- GSEA(gene_list, 
               TERM2GENE    = hallmark,
               minGSSize    = 0,
               maxGSSize    = 1000000,
               pvalueCutoff = 1,
               eps          = 0,
               verbose      = TRUE)
yd_hall <- as.data.frame(y_hall)
qs::qsave(y_hall, file = "analysis/data/08_TNF_NFkB_EMT/01_enrichment/01_hall_gsea.qs")

yd_hall <- yd_hall %>% arrange(desc(NES))
yd_hall$Counts <- lengths(strsplit(yd_hall$core_enrichment, "/"))
yd_hall$Description <- gsub("HALLMARK_", "", yd_hall$Description)
yd_hall$Description <- gsub("_", " ", yd_hall$Description)
yd_hall$Description <- factor(yd_hall$Description, levels = rev(yd_hall$Description))
p <- ggplot(yd_hall[1:10, ], aes(NES, Description)) + 
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
ggsave(filename = "analysis/figure/08_TNF_NFkB_EMT/01_enrichment/02_gsea_hall_all.pdf", p, height = 4.28, width = 8.94)

y_hall@result$Description <- as.character(y_hall@result$Description)

TNFA_NFKB_genes <- c("TNFRSF1A", "RIPK1", "TRADD", "TRAF2", "TRAF5", "BIRC2", "BIRC3",
                     "MAP3K7", "TAB1", "TAB2", "TAB3",
                     "IKBKG", "CHUK", "IKBKB", "NFKB1", "NFKB2", "RELA", "RELB", "TNFAIP3", "NFKBIA", "TNIP1"
)
p <- gseaNb(object = y_hall,
            geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            subPlot = 2,
            geneSize = 4,
            addPval = TRUE,
            addGene = TNFA_NFKB_genes,
            arrowType = 'open')
print(p)


EMT_genes <- c("FN1", "CDH2", 
               "TGFBR3", "WNT5A", "NOTCH2", 
               "SNAI2", "FOXC2", "TWIST1")
EMT_genes <- c("FN1", "CDH2", "VIM", "CDH1",
               "SNAI2", "FOXC2", "ZEB1", "ZEB2")
p <- gseaNb(object = y_hall,
            geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            subPlot = 2,
            geneSize = 4,
            addPval = TRUE,
            addGene = EMT_genes,
            arrowType = 'open')
print(p)
ggsave(filename = "figure/04_tumor_cells_analysis/01_Everolimus_treat_786/03.03_gsea_hall_EMT.pdf", p, height = 4.94, width = 5.64)


#### Rcistarget ####
geneList <- strsplit(yd_hall["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "core_enrichment"], split = "/") %>% unlist()
motifEnrichmentTable_wGenes <- cisTarget(geneList, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations)
