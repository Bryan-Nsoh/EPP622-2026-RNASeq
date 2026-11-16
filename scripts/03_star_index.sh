#!/usr/bin/env bash
set -euo pipefail

REF=data/reference
OUT=results/3_star_index

STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir $OUT --outFileNamePrefix $OUT/ \
  --genomeFastaFiles $REF/GCA_000001735.2_TAIR10.1_genomic.fna \
  --sjdbGTFfile $REF/genomic_modified.gff --sjdbGTFtagExonParentTranscript Parent \
  --genomeSAindexNbases 12 --sjdbOverhang 149
