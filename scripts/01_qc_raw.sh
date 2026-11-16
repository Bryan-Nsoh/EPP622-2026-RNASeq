#!/usr/bin/env bash
set -euo pipefail

OUT=results/1_qc_raw

fastqc -t 12 --outdir $OUT data/reads/*.fastq.gz
multiqc -f $OUT -o $OUT
