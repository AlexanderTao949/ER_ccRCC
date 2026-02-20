#### Load packages ####
# install.packages("moduleColor", update = F, ask = F)
# devtools::install_github("Japrin/sscVis", force = T, upgrade = "never")
# devtools::install_github("Japrin/STARTRAC", force = T, upgrade = "never")
library(scTookit)
library(Seurat)
library(tidyverse)
library(Startrac)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)

#### Load data ####
seu <- qs::qread("analysis/data/04_MPS_subcluster_analysis/02_annotation/03_mps_anno.qs", nthreads = 50)
metadata <- seu@meta.data

#### Run Ro/e ####
Roe <- calTissueDist(metadata,
                     byPatient = F,
                     colname.cluster ="mps_celltype",
                     colname.patient ="orig.ident", 
                     colname.tissue ="ERScore_group", 
                     method ="chisq",  
                     min.rowSum =0) %>% as.data.frame()

Roe <- Roe %>%
  pivot_wider(
    id_cols = Var1,           
    names_from = Var2,        
    values_from = Freq        
  ) %>%
  as.data.frame() %>%        
  column_to_rownames(var = "Var1") %>%  
  dplyr::arrange(desc(High_ERScore)) %>% 
  as.matrix()

qs::qsave(Roe, file = "analysis/data/05_determine_key_mps_subcluster/03_cell_frequency/05_roe/01_roe_res.qs")

col_fun  <-  colorRamp2(c(min(Roe, na.rm  =  TRUE),1,max(Roe, na.rm  =TRUE)),  
                        c("#f6f8e6",  "#f9a33e",  "red"))
p <- Heatmap(Roe,  
             show_heatmap_legend = TRUE, 
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_names_side = 'left',
             column_names_side = "bottom",
             show_column_names = TRUE,
             show_row_names = TRUE,
             col = col_fun,
             width = unit(ncol(Roe_matrix) * 1, "cm"),   
             height = unit(nrow(Roe_matrix) * 0.75, "cm"),    
             row_names_gp = gpar(fontsize = 10),
             column_names_rot = 45,                 
             column_names_gp = gpar(fontsize = 10), 
             heatmap_legend_param = list(
               title = "Ro/e",
               at = c(0, max(Roe_matrix, na.rm = TRUE)),
               labels = c("0", "Max.")
             ),
             cell_fun = function(j, i, x, y, width, height, fill) {
               value <- Roe[i, j] 
               symbol <- if(value == 0) {
                 "−"
               } else if(value > 0 & value < 0.2) {
                 "+/−"
               } else if(value >= 0.2 & value <= 0.8) {
                 "+"
               } else if(value > 0.8 & value <= 1) {
                 "++"
               } else if(value > 1) {
                 "+++"
               }
               
               grid.text(symbol, x, y, gp = gpar(fontsize = 10, col = "black"))
             })

pdf("analysis/figure/05_determine_key_mps_subcluster/03_cell_frequency/05_roe/01_roe_heatmap.pdf", width = 3.65, height = 3.21)
draw(p)
dev.off()
