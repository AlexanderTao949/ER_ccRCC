library(tidygraph)
library(dplyr)
library(igraph)
library(ggkegg)
library(ggfx)
library(cols4all)
library(org.Hs.eg.db)

pathwayplot <- function(pathway_id){
  
  g <- pathway(pathway_id) %>% mutate(Log2FC = assign_deseq2(de, org_db = org.Hs.eg.db),
                                      padj = assign_deseq2(de, column = 'padj', org_db = org.Hs.eg.db),
                                      converted_name = convert_id('hsa'))
  
  p <- ggraph(g, layout = 'manual', x = x, y = y) +
    geom_edge_parallel(width = 0.5,
                       arrow = arrow(length = unit(1, 'mm')),
                       start_cap = square(1, 'cm'),
                       end_cap = square(1.5, 'cm'),
                       aes(color = subtype_name)) +
    geom_node_rect(aes(fill = Log2FC,
                       filter = type == 'gene'),
                   color = 'black') +
    ggfx::with_outer_glow(geom_node_text(aes(label = converted_name,
                                             filter = type != 'group'),
                                         size = 3),
                          colour = 'white', expand = 1) +
    scale_fill_continuous_c4a_div('benedictus', reverse = T) +
    scale_edge_color_manual(values = rev(c4a('pastel', 11))) +
    theme_void() +
    theme(plot.margin = unit(c(0, 1, 0, 0), "cm"))  
  return(p)
}

de <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/01_DESeq2_res.qs")
de <- de %>% column_to_rownames("gene")
gsea_kegg_res <- qs::qread("analysis/data/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/03_gsea_kegg.qs")
gsea_kegg_dat <- as.data.frame(gsea_kegg_res)

p <- pathwayplot("hsa04657")
print(p)
ggsaves("analysis/figure/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/IL-17_pathwayview.pdf", p, width = 17.52, height = 8.18)

p <- pathwayplot("hsa04064")
print(p)
ggsave("analysis/figure/01_Everolimus_signature_scoring/01_Everolimus_treat_enrich_4h/IL-17_pathwayview.pdf", p, width = 17.52, height = 8.18)
