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

## Running the Transcriptomic Analysis

The first-stage DESeq2 workflow is available in `scripts/01_transcriptomics_deseq2.R`. Run it with:

```bash
bash scripts/run_transcriptomics.sh
```

See `scripts/README.md` for inputs, outputs, statistical thresholds, and implementation details.

The second-stage proteomic workflow is available in
`scripts/02_proteomics_limma.R` and can be run with:

```bash
bash scripts/run_proteomics.sh
```

After both single-omics workflows are complete, run the transcriptome-proteome
integration with:

```bash
bash scripts/run_transcriptome_proteome_integration.sh
```

Run ranked transcriptome-proteome GSEA and pathway-level integration with:

```bash
bash scripts/run_pathway_gsea.sh
```
