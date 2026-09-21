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

## Step 4: Plcxd2-Centered Correlation GSEA

Run the Plcxd2-centered co-response analysis with:

```bash
bash scripts/run_plcxd2_correlation_gsea.sh
```

For each omics layer, Pearson correlations are calculated between `Plcxd2`
and every measured gene or protein across all available frontal-cortex
samples. The calculation intentionally does not adjust for Control/Asphyxia,
so the results describe expression patterns that co-vary with `Plcxd2`,
including shared responses to asphyxia. GO Biological Process, Molecular
Function, Cellular Component,
and KEGG GSEA are performed separately for the transcriptome and proteome,
followed by pathway-level integration. Shared
figures retain pathways with `FDR < 0.05` in both layers and concordant NES
directions, displaying up to eight negative and eight positive pathways.
`Plcxd2` itself is retained in each ranked list with its self-correlation fixed
at `r = 1`. A Plcxd2-centered network additionally displays the five strongest
positive and five strongest negative correlates from each omics layer; edge
color indicates correlation direction and edge width indicates absolute
Pearson correlation. A dedicated enrichment curve shows the significant
transcriptomic GSEA result for GO:0042578 (`phosphoric ester hydrolase
activity`) and marks the position of `Plcxd2` in the ranked list. These results
are exploratory and do not establish regulation by `Plcxd2`.

## Step 5: Metabolomics OPLS-DA and VIP Analysis

Run the frontal-cortex metabolomics workflow with:

```bash
bash scripts/run_metabolomics.sh
```

The workflow uses the original-intensity worksheet (`缺失值数据矩阵`) for
detection-rate and QC-RSD assessment and the processed log2 worksheet
(`数据矩阵`) for statistical modeling. Features must be detected in at least
50% of either biological group, detected in all three QC samples, have QC-RSD
no greater than 30%, and retain non-zero biological variance.

OPLS-DA is fitted with `ropls` using Pareto scaling, one predictive component,
one orthogonal component, six-fold cross-validation, and 200 label
permutations. VIP candidates require `VIP > 1`, nominal `P < 0.05`, and
`|log2FC| >= 0.5`. BH-adjusted FDR results are exported separately because VIP
does not replace multiple-testing correction. The workflow exports QC, PCA,
OPLS-DA score, permutation-validation, VIP, S-plot, volcano, heatmap, and
candidate-abundance figures as PDF and 320-dpi PNG files.

## Step 6: Metabolic Pathways and Plcxd2 Integration

Run the local KEGG pathway and Plcxd2-metabolomics integration workflow with:

```bash
bash scripts/run_plcxd2_metabolomics_integration.sh
```

The workflow parses the KEGG compound and pathway annotations already present
in the metabolomics result table. It performs exploratory over-representation
analysis of OPLS-DA VIP candidates. No online query or identifier substitution
is performed.

Per-animal KEGG pathway activity is calculated with single-sample GSEA
(`GSVA::ssGSEA`) using pathways containing at least three annotated
metabolites. Unadjusted Pearson correlations are then calculated between
Plcxd2 abundance and both individual metabolites and ssGSEA pathway scores,
using six RNA-metabolomics matched animals and eight protein-metabolomics
matched animals. The final three-layer network connects Plcxd2 RNA/protein to
directionally concordant pathways and representative member metabolites.
These small-sample correlations describe co-response and do not establish
direct regulation by Plcxd2.
