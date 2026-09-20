#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

conda run --no-capture-output -n zcz_env Rscript \
  scripts/03_transcriptome_proteome_integration.R \
  results/01_transcriptomics/tables/deseq2_all_genes.csv \
  results/02_proteomics/tables/limma_all_proteins.csv \
  results/01_transcriptomics/tables/vst_expression_filtered.csv \
  results/02_proteomics/tables/protein_expression_filtered.csv \
  config/samples.csv \
  results/03_transcriptome_proteome_integration
