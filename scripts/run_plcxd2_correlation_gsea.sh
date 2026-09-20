#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

conda run --no-capture-output -n zcz_env Rscript \
  scripts/04_plcxd2_correlation_gsea.R \
  results/01_transcriptomics/tables/vst_expression_filtered.csv \
  results/02_proteomics/tables/protein_expression_filtered.csv \
  config/samples.csv \
  results/04_plcxd2_correlation_gsea \
  data/reference
