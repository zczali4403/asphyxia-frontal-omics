#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

conda run --no-capture-output -n zcz_env Rscript \
  scripts/02_proteomics_limma.R \
  data/proteomics/protein_expression_log2.xlsx \
  config/samples.csv \
  results/02_proteomics
