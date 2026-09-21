#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

run_analysis() {
  "$@" scripts/05_metabolomics_oplsda.R \
    data/metabolomics/metabolite_expression.xlsx \
    config/samples.csv \
    results/05_metabolomics
}

if [[ "${CONDA_DEFAULT_ENV:-}" == "zcz_env" ]]; then
  run_analysis Rscript
else
  run_analysis conda run --no-capture-output -n zcz_env Rscript
fi
