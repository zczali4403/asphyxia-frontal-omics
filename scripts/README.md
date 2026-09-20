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
- analysis of all input genes without explicit low-count prefiltering or
  DESeq2 independent filtering (all-zero or otherwise non-testable genes remain
  in the exported table with unavailable statistics);
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

## Transcriptome-Proteome Integration

After running the transcriptomic and proteomic workflows, run:

```bash
bash scripts/run_transcriptome_proteome_integration.sh
```

The workflow selects the most abundant representative when multiple features
share a gene symbol, matches RNA and protein results by gene symbol, and
produces a nine-quadrant analysis. RNA evidence is
defined by `padj < 0.05` and `|shrunken log2FC| >= 0.5`; protein evidence is
exploratory and uses nominal `P < 0.05` and `|log2FC| >= 0.5` because no
proteins pass FDR correction.

Individual-level RNA-protein correlations use the six animals with paired
transcriptomic and proteomic data (subjects 1-3 in each group). These
correlations are exported as exploratory tables because of the small paired
sample size. The candidate expression figure instead compares Control and
Asphyxia using all frontal-cortex samples available in each omics dataset.

## Ranked GSEA and Pathway Integration

After completing the transcriptomic and proteomic analyses, run:

```bash
bash scripts/run_pathway_gsea.sh
```

RNA genes are ranked by the DESeq2 Wald statistic, and proteins are ranked by
the limma moderated t statistic. Ranked GSEA is performed independently for
GO Biological Process and KEGG pathways without applying a fold-change cutoff.
GO gene sets use propagated `GOALL` Biological Process annotations so that
genes are included in both directly annotated terms and their parent terms.
KEGG mouse gene-pathway mappings are downloaded once from the KEGG REST API
and cached under `data/reference/` for reproducible reruns.

The pathway integration matches pathways tested in both omics layers and
compares their normalized enrichment scores. Figures and prioritized pathway
tables use nominal `P < 0.05` for both omics layers and are explicitly
exploratory. FDR values and separate FDR-significant tables are retained for
multiple-testing-aware interpretation. Integrated pathway figures show only
shared pathways with concordant RNA and protein NES directions. A dedicated
GSEA running-score plot is also produced for the GO Biological Process
`lipid catabolic process`, with the position of `Plcxd2` highlighted.
