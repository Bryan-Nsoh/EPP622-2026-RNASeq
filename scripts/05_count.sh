#!/usr/bin/env bash
set -euo pipefail

GFF=data/reference/genomic_modified.gff
OUT=results/5_counts

featureCounts -p -s 2 -T 16 -t gene -g ID -F GTF -a $GFF -o $OUT/counts.txt results/4_mapping/*.bam
multiqc -f results/4_mapping $OUT -o $OUT
