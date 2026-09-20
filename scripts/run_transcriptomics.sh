#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

conda run --no-capture-output -n zcz_env Rscript \
  scripts/01_transcriptomics_deseq2.R \
  data/transcriptomics/frontal_cortex_gene_counts.csv \
  config/samples.csv \
  results/01_transcriptomics
