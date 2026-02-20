#### Load pacakges ####
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(viridis)

#### Load data and process ####
seu <- qs::qread("~/projects/Everolimus_Resistance_ccRCC/analysis/data/02_scRNA_data_process/05_pseudobulk_groups/04_seu_group.qs", nthreads = 50)
metadata <- seu@meta.data
pca.embeddings <- Embeddings(seu, reduction = "pca")
harmony.embeddings <- Embeddings(seu, reduction = "harmony")
umap.embeddings <- Embeddings(seu, reduction = "umap")

sce <- SingleCellExperiment(assays=list(counts=seu[["RNA"]]@counts),
                            colData = metadata,
                            reducedDims=SimpleList(PCA=pca.embeddings,
                                                   HARMONY=harmony.embeddings,
                                                   UMAP=umap.embeddings))
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

organism = "human"
if(organism == "human"){
  
  lr_network_all = 
    readRDS("~/public_data/NicheNet_reference/lr_network_human_allInfo_30112033.rds") %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS("~/public_data/NicheNet_reference/ligand_target_matrix_nsga2r_final.rds")
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS("~/public_data/NicheNet_reference/lr_network_mouse_allInfo_30112033.rds") %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS("~/public_data/NicheNet_reference/ligand_target_matrix_nsga2r_final_mouse.rds")
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

#### multinichenetr ####
sample_id = "orig.ident"
group_id = "ERScore_group"
celltype_id = "major_celltype"
covariates = NA
batches = NA

contrasts_oi = c("'High_ERScore-Low_ERScore','Low_ERScore-High_ERScore'") 
contrast_tbl = tibble(
  contrast = c('High_ERScore-Low_ERScore','Low_ERScore-High_ERScore'), 
  group = c("High_ERScore","Low_ERScore"))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)]

conditions_keep = c("High_ERScore","Low_ERScore")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]

min_cells = 10
min_sample_prop = 0.50
fraction_cutoff = 0.05
empirical_pval = FALSE
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE
top_n_target = 250
n.cores = 50
scenario = "regular"
ligand_activity_down = FALSE

multinichenet_output <- multi_nichenet_analysis(
  sce = sce, 
  celltype_id = celltype_id, 
  sample_id = sample_id, 
  group_id = group_id, 
  batches = batches, 
  covariates = covariates, 
  lr_network = lr_network, 
  ligand_target_matrix = ligand_target_matrix, 
  contrasts_oi = contrasts_oi, 
  contrast_tbl = contrast_tbl, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
)

qs::qsave(multinichenet_output, file = "~/projects/Everolimus_Resistance_ccRCC/analysis/data/03_determine_key_major_celltype/01_multinichenetr_major_celltype/01_multinichenet_output.qs")

#### Get data ####
all_lr <- multinichenet_output[["prioritization_tables"]][["group_prioritization_table_source"]] %>% as.data.frame()
TNF_lr <- all_lr %>% dplyr::filter(ligand == "TNF" & receiver == "Tumor_Cells" & direction_regulation == "up")

plot_data <- TNF_lr %>% dplyr::select(group, sender, receiver, ligand, receptor, fraction_ligand_group, fraction_receptor_group, prioritization_score)
plot_data$sender_receiver <- paste(plot_data$sender, plot_data$receiver, sep = " - ")
plot_data$ligand_receptor <- paste(plot_data$ligand, plot_data$receptor, sep = " -> ")
sender_receiver_levels <- plot_data %>% dplyr::filter(receptor == "TNFRSF1A" & group == "High_ERScore") %>% dplyr::arrange(desc(prioritization_score)) %>% pull(sender_receiver)
plot_data <- plot_data %>%
  dplyr::mutate(
    sender_receiver = factor(sender_receiver, levels = rev(sender_receiver_levels)),
    group = factor(group, levels = c("High_ERScore", "Low_ERScore")),
    ligand_receptor = factor(ligand_receptor, levels = c("TNF -> TNFRSF1A", "TNF -> LTBR", "TNF -> TNFRSF21", "TNF -> TRAF2"))
  )

p <- ggplot(plot_data, aes(x = group, y = sender_receiver)) +
  geom_point(aes(size = fraction_ligand_group,
                 color = prioritization_score),
             alpha = 0.9,
             stroke = 0.5) +
  facet_wrap(~ ligand_receptor, ncol = 4, scales = "fixed") +

  scale_color_viridis_c(
    option = "plasma",
    direction = -1,
    begin = 0.1,
    end = 0.9,
    name = "Prioritization Score"
  ) +

  scale_size_continuous(
    range = c(2, 8),
    name = "Ligand Expression Proportion",
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 10, hjust = 1),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    strip.background = element_rect(fill = "grey95", color = "black", size = 0.8),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p)
ggsave("analysis/figure/03_determine_key_major_celltype/01_multinichenetr_major_celltype/multinichenetr_TNF_dotplot.pdf", p, width = 10.14, height = 3.43)
