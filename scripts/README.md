# Analysis Scripts

Run all scripts from the repository root in an R environment containing the
required packages.

## Step 1: Frontal-Cortex Transcriptomics

Run the complete DESeq2 workflow from the repository root:

```bash
Rscript scripts/01_transcriptomics_deseq2.R
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

### Step 1 Supplement: Separate GO and KEGG Enrichment

After Step 1, run:

```bash
Rscript scripts/01_transcriptomics_go_kegg_enrichment.R
```

This standalone analysis uses the Step 1 effect-size-filtered DEG table
(`FDR < 0.05`, `|shrunken log2FC| >= 0.5`). GO BP/MF/CC annotations come from
`org.Mm.eg.db` and KEGG annotations from the local mouse reference in
`data/reference/`. Each ontology/database uses its own annotated, measured
genes with finite DESeq2 P values as background. Hypergeometric tests use
sets of 10–500 measured genes; BH correction is applied separately within
BP, MF, CC and KEGG, including eligible terms with zero overlap.

Outputs are saved under `results/01_transcriptomics/enrichment/`. GO is shown
as a faceted dot plot (up to eight FDR-significant terms per ontology, gene
ratio on the x-axis, overlap count as point size). KEGG is shown as a horizontal
bar plot (up to fifteen FDR-significant pathways, `-log10(FDR)` on the x-axis
and gene counts alongside bars). Both figures use color intensity for
enrichment FDR and are exported as PDF and 320-dpi PNG. If no terms pass FDR,
the corresponding figure states that no enriched terms were found.

The three optional positional arguments are the Step 1 result directory,
the KEGG reference directory, and the enrichment output directory. Full
tables, FDR-filtered tables, plotted terms, gene mapping, background counts,
parameters, summary and session information are also exported. This script
uses transcriptomic results only and does not run Step 7 or differential
expression again. Its per-ontology/database backgrounds differ from the
combined annotation background used by Step 7, so enrichment statistics can
differ.

## Step 2: Frontal-Cortex Proteomics

Run the limma workflow from the repository root:

```bash
Rscript scripts/02_proteomics_limma.R
```

The workflow reads the `数据矩阵` worksheet from
`data/proteomics/protein_expression_log2.xlsx`. It retains proteins with an
original missing fraction no greater than 50% in both groups, performs an
unpaired limma comparison of Asphyxia versus Control, and uses BH correction.
The primary protein set requires `padj < 0.05` and `|log2FC| >= 0.5`.

Proteins meeting `P value < 0.05` are exported separately as exploratory
results. The workflow produces detection,
distribution, PCA, correlation, volcano, MA, and protein heatmap figures as
PDF and 320-dpi PNG files.

### Step 2 Supplement: Separate GO and KEGG Enrichment

After Step 2, run:

```bash
Rscript scripts/02_proteomics_go_kegg_enrichment.R
```

The standalone workflow uses `limma_exploratory_nominal_p_proteins.csv`,
selected by `P value < 0.05` and `|log2FC| >= 0.5`. Protein `GeneID` values are
validated as mouse Entrez IDs; query and background are deduplicated by gene
ID. The background consists of proteins with finite limma P values and valid
annotations in the respective GO ontology or KEGG database. GO BP/MF/CC
annotations use `org.Mm.eg.db` GOALL; KEGG uses the local mouse reference.
Hypergeometric tests use sets of 10–500 measured genes, with BH correction
within BP, MF, CC and KEGG including eligible zero-overlap terms.

Outputs are saved under `results/02_proteomics/enrichment/`: GO is a faceted
dot plot with up to eight FDR-significant terms per ontology, and KEGG is a
horizontal bar plot with up to fifteen FDR-significant pathways. Color
intensity represents enrichment FDR. Counts and ratios refer to unique
protein-associated genes. Both figures are exported as PDF and 320-dpi PNG;
if no terms meet FDR < 0.05, the corresponding figure states this explicitly.
Full results, FDR-filtered tables, plotted terms, protein mapping, background
counts, parameters, summary and session information are also exported.

The three optional positional arguments are the Step 2 result directory,
the KEGG reference directory, and the enrichment output directory. The
script uses protein results only. Its per-ontology/database backgrounds
differ from the combined annotation background in Step 7, so statistics can
differ. Enrichment FDR does not change the statistical status of the input
protein candidates.

## Step 3: Transcriptome-Proteome Integration

After running the transcriptomic and proteomic workflows, run:

```bash
Rscript scripts/03_transcriptome_proteome_integration.R
```

The workflow selects the most abundant representative when multiple features
share a gene symbol, matches RNA and protein results by gene symbol, and
produces a nine-quadrant analysis. RNA evidence is
defined by `padj < 0.05` and `|shrunken log2FC| >= 0.5`; protein evidence is
exploratory and uses `P value < 0.05` and `|log2FC| >= 0.5` because no
proteins pass FDR correction.

Individual-level RNA-protein correlations use the six animals with paired
transcriptomic and proteomic data (subjects 1-3 in each group). These
correlations are exported as exploratory tables because of the small paired
sample size. The candidate expression figure instead compares Control and
Asphyxia using all frontal-cortex samples available in each omics dataset.

## Step 4: Plcxd2-Centered Correlation GSEA

Run the Plcxd2-centered co-response analysis with:

```bash
Rscript scripts/04_plcxd2_correlation_gsea.R
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
Rscript scripts/05_metabolomics_oplsda.R
```

The workflow uses the original-intensity worksheet (`缺失值数据矩阵`) for
detection-rate and QC-RSD assessment and the processed log2 worksheet
(`数据矩阵`) for statistical modeling. Features must be detected in at least
50% of either biological group, detected in all three QC samples, have QC-RSD
no greater than 30%, and retain non-zero biological variance.

OPLS-DA is fitted with `ropls` using Pareto scaling, one predictive component,
one orthogonal component, six-fold cross-validation, and 200 label
permutations. VIP candidates require `VIP > 1`, uncorrected `P < 0.05`, and
`|log2FC| >= 0.5`. BH-adjusted FDR results are exported separately because VIP
does not replace multiple-testing correction. The workflow exports QC, PCA,
OPLS-DA score, permutation-validation, VIP, S-plot, volcano, heatmap, and
candidate-abundance figures as PDF and 320-dpi PNG files.

## Step 6: Metabolic Pathways and Plcxd2 Integration

Run the local KEGG pathway and Plcxd2-metabolomics integration workflow with:

```bash
Rscript scripts/06_plcxd2_metabolomics_integration.R
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

## Step 7: DEG and Exploratory DEP Pathway Enrichment

After Steps 1 and 2, run:

```bash
Rscript scripts/07_differential_feature_enrichment.R
```

The DEG set requires `FDR < 0.05` and `|shrunken log2FC| >= 0.5`.
Because no proteins pass the differential-expression FDR threshold, the
protein set is exploratory (`P value < 0.05`, `|log2FC| >= 0.5`).
Over-representation tests use separately measured/tested genes or proteins
as the respective background, GO BP/MF/CC annotations from `org.Mm.eg.db`,
and the local mouse KEGG reference in `data/reference/`. Gene sets contain
10–500 measured genes; BH correction is performed within each annotation
class. The workflow exports full and FDR-filtered tables, plus a separate
DEG lollipop figure and an exploratory DEP enrichment bar plot, with bars
colored by GO BP, GO MF, GO CC, or KEGG pathway class. Enrichment FDR values
do not make the input protein selection FDR-significant.
