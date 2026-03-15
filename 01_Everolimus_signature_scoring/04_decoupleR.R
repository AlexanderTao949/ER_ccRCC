Sys.setenv(http_proxy = "http://iyun70.com:7890")
Sys.setenv(https_proxy = "http://iyun70.com:7890") 
#### Load pacakges ####
# devtools::install_github('saezlab/OmnipathR', ref = '041dace', dependencies = TRUE, upgrade = "never", force = T)
# devtools::install_version("httr2", version = "1.0.2", dependencies = TRUE, upgrade = "never", force = T)
library(decoupleR)
library(tidyverse)
library(tibble)
library(ggplot2)
library(pheatmap)
library(ggrepel)

#### Load data and process ####
net <- get_progeny(organism = 'human', top = 500)
tfnet <- get_collectri(organism = 'human', split_complexes = FALSE)
qs::qsave(net, file = "analysis/resource/decoupleR_PROGENy_pathway_net.qs")
qs::qsave(tfnet, file = "analysis/resource/decoupleR_CollecTRI_TF_net.qs")
exprSet <- qs::qread("RNA_seq_upstream/GSE99875/02_annotated_expression_matrices/vsd.qs")
exprSet <- exprSet %>%
  dplyr::mutate_if(~ any(is.na(.x)), 
                   ~ dplyr::if_else(is.na(.x), 0, .x)) %>% as.matrix() 
metadata <- qs::qread("RNA_seq_upstream/GSE99875/resource/metadata.qs")
de_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/01_DESeq2_res.qs")
de_res$stat <- ifelse(de_res$lfcSE == 0, 0, 
                      de_res$log2FoldChange / de_res$lfcSE)
de_res <- de_res %>% dplyr::select(gene, stat) %>% column_to_rownames("gene") %>% as.matrix()

#### pathway activity ####
sample_acts <- decoupleR::run_mlm(mat = exprSet, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05,2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20)

#### DE pathway activity ####
contrast_acts <- decoupleR::run_mlm(mat  =de_res, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)
# Plot
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("Pathways")

print(p)

pathway <- 'NFkB'

df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(de_res), rownames(df)))

df <- df[inter, ]

df['t_value'] <- de_res[inter, ]

df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = df, 
                     mapping = ggplot2::aes(x = weight, 
                                            y = t_value, 
                                            color = color)) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(pathway)

p

#### TF activity ####
sample_tf_acts <- decoupleR::run_ulm(mat = exprSet, 
                                     net = tfnet, 
                                     .source = 'source', 
                                     .target = 'target',
                                     .mor = 'mor', 
                                     minsize = 5)
n_tfs <- 25

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20)