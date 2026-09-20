#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

conda run --no-capture-output -n zcz_env Rscript \
  scripts/04_pathway_gsea_integration.R \
  results/01_transcriptomics/tables/deseq2_all_genes.csv \
  results/02_proteomics/tables/limma_all_proteins.csv \
  results/04_transcriptome_proteome_pathways \
  data/reference
