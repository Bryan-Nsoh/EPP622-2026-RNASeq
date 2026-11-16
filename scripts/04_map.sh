#!/usr/bin/env bash
set -euo pipefail

IN=results/2_trimmed
INDEX=results/3_star_index
OUT=results/4_mapping

for s in Col_0h_rep1 Col_0h_rep2 Col_0h_rep3 Col_3h_rep1 Col_3h_rep2 Col_3h_rep3; do
  STAR --runThreadN 16 --genomeDir $INDEX \
    --readFilesIn $IN/${s}_1.fastq.gz $IN/${s}_2.fastq.gz --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $OUT/$s.
  mv $OUT/$s.Aligned.sortedByCoord.out.bam $OUT/$s.bam
  samtools index $OUT/$s.bam
  samtools stats $OUT/$s.bam > $OUT/$s.stats.txt
done
multiqc -f $OUT -o $OUT
