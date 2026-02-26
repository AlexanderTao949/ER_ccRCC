#### Load pacakges ####
library(monocle)
library(Seurat)
library(tidyverse)
library(data.table)

#### Load data ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)

#### Build CellDataSet ####
metadata <- seu@meta.data
pdata <- metadata[, c("orig.ident", "ERScore_group", "mps_celltype")]
fData <- data.frame(gene_short_name = row.names(seu@assays$RNA@counts), row.names = row.names(seu@assays$RNA@counts))
fd <- new("AnnotatedDataFrame",
          data=fData)
pd <- new("AnnotatedDataFrame",
          data=pdata)
HSMM  <- newCellDataSet(as(as.matrix(seu@assays$RNA@counts), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size(),
                      lowerDetectionLimit = 0.5)


#### Estimate size factors and dispersions ####
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

#### QC ####
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= nrow(metadata) *0.005))
valid_cells <- row.names(subset(pData(HSMM),
                                num_genes_expressed >= 100))
HSMM <- HSMM[,valid_cells]
qs::qsave(HSMM, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/07_mps_trajectory/01_monocle2/HSMM.qs")
qs::qsave(expressed_genes, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/07_mps_trajectory/01_monocle2/expressed_genes.qs")
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
melted_dens_df <- reshape2::melt(t(scale(t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

#### Constructing Single Cell Trajectories ####
HSMM_myo <- HSMM
##### Step 1: choosing genes that define progress ####
cat("Step 1: choosing genes that define progress")
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~mps_celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)

##### Step 2: reducing the dimensionality of the data #####
cat("Step 2: reducing the dimensionality of the data")
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree')

##### Step 3: ordering the cells in pseudotime #####
cat("Step 3: ordering the cells in pseudotime")
HSMM_myo <- orderCells(HSMM_myo)

qs::qsave(HSMM_myo, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/07_mps_trajectory/01_monocle2/01_monocle2_res.qs")


#### Plot ####
mps_colors <- c(
  "Mono_CD14_S100A8" = "#4DBBD5", 
  "Mono_CD16_TCF7L2" = "#00798C", 
  "Maph_Inflam_CCL3" = "#3CB371", 
  "Maph_LA_APOE"     = "#195F43", 
  "cDC2_FCER1A"      = "#A2C865", 
  "cDC1_CCSER1"      = "#6E8B3D"  
)
p1 <- plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime", show_branch_points = FALSE)
p2 <- plot_cell_trajectory(HSMM_myo, color_by = "mps_celltype", show_branch_points = FALSE) + 
  scale_color_manual(values = mps_colors) +
  theme(legend.position = "right")
p <- p1+p2
print(p)

# #### Change State ####
# mycds <- orderCells(mycds,root_state = 2)
# qs::qsave(mycds, file = file.path(homepath, "data/08_MPS_trajectory/0802_monocle/0802_01_monocle_change_state_res.qs"))
# 
# mycds <- qs::qread("data/08_MPS_trajectory/0802_monocle/0802_01_monocle_change_state_res.qs")
# p1 <- plot_cell_trajectory(mycds, color_by = "State")
# p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
# p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
# plot <- p1+p2+p3

### Save to Seurat ####
pdata <- Biobase::pData(HSMM_myo)
pdata <- pdata %>% dplyr::select(Pseudotime, State)
colnames(pdata) <- c("monocle2_ptime", "monocle2_state")
qs::qsave(pdata, file = "analysis/data/07_mps_trajectory/01_monocle2/02_ptime_state.qs")

