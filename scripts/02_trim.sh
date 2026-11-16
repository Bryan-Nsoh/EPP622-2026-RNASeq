#!/usr/bin/env bash
set -euo pipefail

IN=data/reads
OUT=results/2_trimmed

for s in Col_0h_rep1 Col_0h_rep2 Col_0h_rep3 Col_3h_rep1 Col_3h_rep2 Col_3h_rep3; do
  fastp -w 16 --detect_adapter_for_pe --trim_poly_g \
    -i $IN/${s}_1.fastq.gz -I $IN/${s}_2.fastq.gz \
    -o $OUT/${s}_1.fastq.gz -O $OUT/${s}_2.fastq.gz \
    -h $OUT/$s.fastp.html -j $OUT/$s.fastp.json
done
fastqc -t 12 --outdir $OUT $OUT/*.fastq.gz
multiqc -f $OUT -o $OUT
