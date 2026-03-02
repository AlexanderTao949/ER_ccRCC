Sys.setenv(http_proxy = "http://iyun70.com:7890")
Sys.setenv(https_proxy = "http://iyun70.com:7890") 
#### Load packages ####
library(Seurat)
library(tidyverse)
library(scTookit)
library(presto)
library(clusterProfiler)
library(org.Hs.eg.db)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs")
protein_coding_genes <- qs::qread("~/public_data/genome_files/GRCh38_gtf/protein_coding_genes_v49.qs") %>% pull(gene_name)

#### FindMarkers ####
de <- wilcoxauc(seu,
                group_by = "mps_celltype",
                seurat_assay = "RNA",
                assay = "data",
                verbose = TRUE)
de <- de %>% 
  group_by(group) %>% 
  dplyr::arrange(desc(logFC), .by_group = TRUE)
qs::qsave(de, file = "analysis/data/04_MPS_subcluster_analysis/03_function_analysis/01_FindMarkers_celltype.qs")

#### compareCluster KEGG ####
gene_mapping <- bitr(de$feature, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = "org.Hs.eg.db")
colnames(gene_mapping) <- c("feature", "entrez")
gene_mapping <- gene_mapping[gene_mapping$feature %in% protein_coding_genes, ]
gene_df <- de %>% dplyr::inner_join(gene_mapping, by = "feature", relationship = "many-to-many")

#### ORA ####
input_list <- lapply(unique(gene_df$group), function(x){
  dat <- gene_df %>% dplyr::filter(group == x)
  g <- dat %>% 
    dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(desc(logFC)) %>% 
    slice_head(n = 150) %>% 
    pull(entrez)
  return(g)
})
names(input_list) <- unique(gene_df$group)

kegg_ora <- compareCluster(geneClusters = input_list,
                           fun = "enrichKEGG",
                           organism = "hsa",
                           keyType = "ncbi-geneid",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           universe = gene_mapping$entrez,
                           qvalueCutoff = 0.2)
kegg_ora <- setReadable(kegg_ora,
                        OrgDb = org.Hs.eg.db,
                        keyType="ENTREZID")
kegg_ora@compareClusterResult <- subset(kegg_ora@compareClusterResult, p.adjust < 0.05)
kegg_ora_dat <- as.data.frame(kegg_ora)
qs::qsave(kegg_ora, file = "analysis/data/04_MPS_subcluster_analysis/03_function_analysis/02_kegg_ora.qs")

#### GSEA ####
input_list <- lapply(unique(gene_df$group), function(x){
  dat <- gene_df %>% dplyr::filter(group == x)
  entrez_list <- dat$logFC
  names(entrez_list) <- dat$entrez
  entrez_list <- sort(entrez_list, decreasing = TRUE)
  return(entrez_list)
})
names(input_list) <- unique(gene_df$group)

kegg_gsea <- compareCluster(geneCluster = input_list,
                            fun = "gseKEGG",
                            organism = "hsa",
                            keyType = "ncbi-geneid",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            verbose = TRUE,
                            pvalThreshold = 0.05)
kegg_gsea <- setReadable(kegg_gsea,
                         OrgDb = org.Hs.eg.db,
                         keyType="ENTREZID")
kegg_gsea@compareClusterResult <- subset(kegg_gsea@compareClusterResult, p.adjust < 0.05)
kegg_gsea_dat <- as.data.frame(kegg_gsea)
qs::qsave(kegg_gsea, file = "analysis/data/04_MPS_subcluster_analysis/03_function_analysis/02_kegg_gsea.qs")
