#!/bin/bash
set -eo pipefail

## inputs
f_loom_grn=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/data/06_GRN/01_pySCENIC_process/mc_mat_for_step1.loom

## outputs
grn_output=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/data/06_GRN/01_pySCENIC_process/pyscenic_step1_adj.tsv
ctx_output=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/data/06_GRN/01_pySCENIC_process/pyscenic_step2_reg.tsv

## reference (hg38 or mm10)
f_tfs=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/resource/cisTargetDB/allTFs_hg38.txt
f_motif_path=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/resource/cisTargetDB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
f_db_500bp=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/resource/cisTargetDB/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_db_10kb=/media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/resource/cisTargetDB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


#### 1. build GRN
time arboreto_with_multiprocessing.py \
$f_loom_grn \
$f_tfs \
--method grnboost2 \
--output $grn_output \
--num_workers 25 


#### 2. cisTarget
time pyscenic ctx \
$grn_output \
$f_db_500bp $f_db_10kb \
--annotations_fname $f_motif_path \
--expression_mtx_fname $f_loom_grn \
--output $ctx_output \
--num_workers 25

python /media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/script/06_GRN/01_pySCENIC_process/regulon2gmt.py /media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/data/06_GRN/01_pySCENIC_process/pyscenic_step2_reg.tsv /media/desk16/iyun4605/projects/Everolimus_Resistance_ccRCC/analysis/data/06_GRN/01_pySCENIC_process/pyscenic_TF