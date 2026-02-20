#### Load pacakges ####
library(Seurat)
library(tidyverse)
library(UCell)
library(scTookit)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
keggdb <- qs::qread("~/public_data/kegg_db/kegg_human_pathway.qs")
TNF_geneset <- list(TNF_signaling = keggdb %>% dplyr::filter(term == "TNF signaling pathway") %>% pull(gene))

#### UCell ####
DefaultAssay(seu) <- "RNA"
seu <- AddModuleScore_UCell(obj        = seu, 
                            features   = TNF_geneset, 
                            storeRanks = TRUE,
                            ncores     = 25,
                            force.gc   = TRUE,
                            name       = "")
metadata <- seu@meta.data
qs::qsave(metadata, file = "analysis/data/05_determine_key_mps_subcluster/04_TNF_signaling_UCell/01_TNF_signaling_UCell_metadata.qs")

#### umap ####
seu_list <- SplitObject(seu, split.by = "ERScore_group")
cuttoff <- median(seu$TNF_signaling)
p1 <- ThreshPlot(seu = seu_list[["High_ERScore"]], feature = "TNF_signaling", threshold = cuttoff, legend.position = "right", dot.size = 1, title = "High_ERScore", title.size = 12)
p2 <- ThreshPlot(seu = seu_list[["Low_ERScore"]], feature = "TNF_signaling", threshold = cuttoff, legend.position = "right", dot.size = 1, title = "Low_ERScore", title.size = 12)
p <- patchwork::wrap_plots(list(p1, p2), ncol = 2)
print(p)
ggsave("analysis/figure/05_determine_key_mps_subcluster/04_TNF_signaling_UCell/01_TNF_signaling_score_group_umap.pdf", p, width = 10.68, height = 4.08)

#### Violin Plot ####
seu_list <- SplitObject(seu, split.by = "mps_celltype")
