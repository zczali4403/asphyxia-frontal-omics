# Analysis Scripts

## Frontal-Cortex Transcriptomics

Run the complete DESeq2 workflow from the repository root:

```bash
bash scripts/run_transcriptomics.sh
```

The wrapper uses the `zcz_env` Conda environment. The underlying R script can also be called directly:

```bash
conda run --no-capture-output -n zcz_env Rscript \
  scripts/01_transcriptomics_deseq2.R \
  data/transcriptomics/frontal_cortex_gene_counts.csv \
  config/samples.csv \
  results/01_transcriptomics
```

The three positional arguments are the count matrix, sample metadata, and output directory. If omitted, the script uses the paths shown above.

The workflow performs:

- validation of sample names, groups, and integer counts;
- low-expression filtering (`count >= 10` in at least three samples);
- DESeq2 normalization and differential expression for Asphyxia versus Control;
- adaptive log2-fold-change shrinkage with `apeglm`;
- PCA, sample-correlation, library-size, and VST-distribution plots;
- volcano, MA, and differential-gene heatmap plots;
- export of full results, FDR-significant genes, effect-size-filtered genes, normalized counts, VST expression, PCA scores, and session information.

The workflow reports two nested differential-expression sets. The FDR set uses `padj < 0.05`; the effect-size-filtered set additionally requires an absolute shrunken log2 fold change of at least 0.5 (approximately a 1.41-fold change).

Each figure is exported as a publication-ready vector PDF and a 320-dpi PNG.

## Frontal-Cortex Proteomics

Run the limma workflow from the repository root:

```bash
bash scripts/run_proteomics.sh
```

The workflow reads the `数据矩阵` worksheet from
`data/proteomics/protein_expression_log2.xlsx`. It retains proteins with an
original missing fraction no greater than 50% in both groups, performs an
unpaired limma comparison of Asphyxia versus Control, and uses BH correction.
The primary protein set requires `padj < 0.05` and `|log2FC| >= 0.5`.

Because nominal-P findings do not control the false-discovery rate, they are
exported separately as exploratory results. The workflow produces detection,
distribution, PCA, correlation, volcano, MA, and protein heatmap figures as
PDF and 320-dpi PNG files.
