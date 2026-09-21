# Multi-Omics Analysis of Mechanical Asphyxia Effects on the Mouse Frontal Cortex

This repository contains transcriptomic, proteomic, and metabolomic data for investigating the effects of mechanical asphyxia on the mouse frontal cortex. The study groups are labeled `Control` and `Asphyxia`, and the frontal cortex is designated as region `R1` in the source data. Original sample identifiers retain the legacy `CON` prefix.

## Study Design

| Omics layer | Asphyxia | Control | Tissue |
|---|---:|---:|---|
| Transcriptomics | 3 | 3 | Frontal cortex (R1) |
| Proteomics | 4 | 4 | Frontal cortex (R1) |
| Metabolomics | 6 | 6 | Frontal cortex (R1), plus 3 QC samples |

Samples with the same group and replicate number originate from the same animal. For example, transcriptomic samples `CON_R1_1` and `Asphyxia_R1_1` correspond to subjects `CON-1` and `Asphyxia-1`, respectively.

- Complete three-omics pairing: three animals per group, `CON-1` through `CON-3` and `Asphyxia-1` through `Asphyxia-3`.
- Proteomics-metabolomics pairing: four animals per group, `CON-1` through `CON-4` and `Asphyxia-1` through `Asphyxia-4`.
- Metabolomics-only animals: `CON-5`, `CON-6`, `Asphyxia-5`, and `Asphyxia-6`.

## Repository Structure

```text
.
├── config/
│   └── samples.csv
├── data/
│   ├── transcriptomics/
│   ├── proteomics/
│   ├── metabolomics/
│   ├── metadata/
│   └── README.md
├── scripts/
└── results/
```

## Recommended Input Matrices

- Transcriptomic differential expression: `data/transcriptomics/frontal_cortex_gene_counts.csv`
- Transcriptomic visualization and multi-omics integration: `data/transcriptomics/frontal_cortex_vst_expression.csv`
- Proteomics: the processed expression-matrix worksheet in `data/proteomics/protein_expression_log2.xlsx`
- Metabolomics: the processed expression-matrix worksheet in `data/metabolomics/metabolite_expression.xlsx`

`all_regions_vst_expression.csv` contains transcriptomic samples from all four regions (R1-R4) and can be used for region-specific comparisons. Frontal-cortex multi-omics integration should use `frontal_cortex_vst_expression.csv`.

## Analysis Notes

- Use the raw transcriptomic count matrix with count-based models such as DESeq2. Do not use VST values as DESeq2 input.
- The existing proteomic differential list was selected using nominal `P < 0.05` and fold-change thresholds, but no proteins pass `q < 0.05` after multiple-testing correction.
- The existing metabolomic differential list likewise contains no metabolites passing `q < 0.05`; these findings should be treated as exploratory.
- Large raw FASTQ, mass-spectrometry RAW, and tar archive files are intentionally excluded from this repository.
- Metadata retain both the original region code (`R1`) and the standardized region name (`frontal_cortex`).
- Individual-level three-omics analyses should use the six completely paired animals. Proteomics-metabolomics paired analyses can use eight animals.
- Chinese biological annotations and original worksheet names in vendor-provided data files are intentionally preserved.

See `data/README.md` for data provenance, worksheet details, and file checksums.

## Analysis Workflow

Run all commands from the repository root. Each wrapper uses the `zcz_env`
Conda environment and writes results to the corresponding numbered directory
under `results/`.

| Step | Analysis | Command | Main output directory |
|---:|---|---|---|
| 1 | Transcriptomic differential expression with DESeq2 | `bash scripts/run_transcriptomics.sh` | `results/01_transcriptomics/` |
| 2 | Proteomic differential expression with limma | `bash scripts/run_proteomics.sh` | `results/02_proteomics/` |
| 3 | Transcriptome-proteome integration and nine-quadrant analysis | `bash scripts/run_transcriptome_proteome_integration.sh` | `results/03_transcriptome_proteome_integration/` |
| 4 | Plcxd2-centered transcriptome/proteome correlation GSEA | `bash scripts/run_plcxd2_correlation_gsea.sh` | `results/04_plcxd2_correlation_gsea/` |
| 5 | Metabolomic OPLS-DA, permutation validation, and VIP analysis | `bash scripts/run_metabolomics.sh` | `results/05_metabolomics/` |
| 6 | Plcxd2-metabolomics pathway integration | `bash scripts/run_plcxd2_metabolomics_integration.sh` | `results/06_plcxd2_metabolomics_integration/` |

Step 4 ranks genes and proteins by their unadjusted Pearson correlation with
`Plcxd2`; it does not regress out the Control/Asphyxia group effect. GSEA
multiple-testing results are nevertheless reported with BH-adjusted FDR, and
the shared-pathway summary retains concordant pathways significant in both
omics layers.

Step 5 applies feature-level detection and QC filters before fitting an
OPLS-DA model. VIP candidates require `VIP > 1`, nominal `P < 0.05`, and
`|log2FC| >= 0.5`. OPLS-DA separation and VIP values are exploratory and must
be interpreted together with cross-validation, permutation testing, and
multiple-testing results.

Step 6 uses the KEGG annotations already present in the metabolomics table.
It performs VIP-candidate pathway over-representation analysis, calculates
per-animal pathway activity with ssGSEA, and correlates Plcxd2 RNA and protein
abundance with metabolites and pathway scores in matched animals. Its final
network connects Plcxd2 RNA/protein, selected pathways, and representative
VIP metabolites. Nominal pathway enrichment (`P < 0.05`) is distinguished
from FDR significance in the exported tables and figures.

These Plcxd2-centered correlations and networks represent cross-omics
co-variation, not evidence that Plcxd2 directly regulates a metabolite or
pathway. See `scripts/README.md` for detailed inputs, thresholds, outputs, and
implementation notes.
