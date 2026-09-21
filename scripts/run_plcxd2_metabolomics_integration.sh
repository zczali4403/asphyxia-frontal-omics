#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

run_analysis() {
  "$@" scripts/06_plcxd2_metabolomics_integration.R \
    results/05_metabolomics/tables/oplsda_vip_all_metabolites.csv \
    results/05_metabolomics/tables/metabolite_expression_filtered.csv \
    results/01_transcriptomics/tables/vst_expression_filtered.csv \
    results/02_proteomics/tables/protein_expression_filtered.csv \
    config/samples.csv \
    results/06_plcxd2_metabolomics_integration
}

if [[ "${CONDA_DEFAULT_ENV:-}" == "zcz_env" ]]; then
  run_analysis Rscript
else
  run_analysis conda run --no-capture-output -n zcz_env Rscript
fi

