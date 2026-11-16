#!/usr/bin/env bash
set -euo pipefail

DATA=/lustre/isaac24/proj/UTK0505/module4_rnaseq_data/data
cd /lustre/isaac24/proj/UTK0505/bnsoh/test4

mkdir -p data/reads data/reference scripts logs
mkdir -p results/0_sample_table results/1_qc_raw results/2_trimmed results/3_star_index results/4_mapping
mkdir -p results/5_counts results/6_deseq2 results/7_heatmap results/8_cold_genes

ln -sf $DATA/raw/*.fastq.gz data/reads/
ln -sf $DATA/raw/metadata.txt data/reads/
ln -sf $DATA/reference/GCA_000001735.2_TAIR10.1_genomic.fna data/reference/
ln -sf $DATA/reference/genomic_modified.gff data/reference/

cat > results/0_sample_table/samples.csv <<TABLE
sample,time,replicate,sra_run
Col_0h_rep1,0h,rep1,SRR18313928
Col_0h_rep2,0h,rep2,SRR18313927
Col_0h_rep3,0h,rep3,SRR18313918
Col_3h_rep1,3h,rep1,SRR18313902
Col_3h_rep2,3h,rep2,SRR18313922
Col_3h_rep3,3h,rep3,SRR18313921
TABLE
