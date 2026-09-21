#!/usr/bin/env Rscript

# Differential transcriptomic analysis of the R1 frontal cortex.
# Comparison: Asphyxia vs Control.

required_packages <- c(
  "DESeq2", "apeglm", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- if (length(args) >= 1) args[[1]] else
  "data/transcriptomics/frontal_cortex_gene_counts.csv"
metadata_file <- if (length(args) >= 2) args[[2]] else
  "config/samples.csv"
output_dir <- if (length(args) >= 3) args[[3]] else
  "results/01_transcriptomics"

# Primary statistical thresholds.
alpha <- 0.05
lfc_threshold <- 0.5
top_labels <- 12L
top_heatmap <- 50L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
qc_dir <- file.path(figures_dir, "qc")
de_dir <- file.path(figures_dir, "differential_expression")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading count matrix: ", counts_file)
counts_input <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_count_columns <- c("GeneID", "GeneName")
if (!all(required_count_columns %in% colnames(counts_input))) {
  stop("The count matrix must contain GeneID and GeneName columns.", call. = FALSE)
}
if (anyDuplicated(counts_input$GeneID)) {
  stop("GeneID values must be unique.", call. = FALSE)
}

message("Reading sample metadata: ", metadata_file)
metadata_all <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_metadata_columns <- c(
  "omics", "sample_id", "subject_id", "group", "region", "sample_role"
)
if (!all(required_metadata_columns %in% colnames(metadata_all))) {
  stop(
    "The metadata table is missing required columns: ",
    paste(setdiff(required_metadata_columns, colnames(metadata_all)), collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata_all[
  metadata_all$omics == "transcriptomics" &
    metadata_all$sample_role == "biological" &
    metadata_all$region == "frontal_cortex",
  ,
  drop = FALSE
]

if (nrow(metadata) == 0) {
  stop("No frontal-cortex transcriptomic samples were found in metadata.", call. = FALSE)
}
if (anyDuplicated(metadata$sample_id)) {
  stop("Transcriptomic sample_id values must be unique.", call. = FALSE)
}

sample_columns <- metadata$sample_id
missing_samples <- setdiff(sample_columns, colnames(counts_input))
extra_samples <- setdiff(
  setdiff(colnames(counts_input), required_count_columns),
  sample_columns
)
if (length(missing_samples) > 0) {
  stop(
    "Samples in metadata but absent from the count matrix: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}
if (length(extra_samples) > 0) {
  warning(
    "Count columns not selected by the metadata were ignored: ",
    paste(extra_samples, collapse = ", ")
  )
}

# Standardize the display label while retaining compatibility with legacy metadata.
metadata$group[metadata$group == "CON"] <- "Control"
metadata$group <- factor(metadata$group, levels = c("Control", "Asphyxia"))
if (anyNA(metadata$group)) {
  stop("Transcriptomic groups must be Control (or legacy CON) and Asphyxia.", call. = FALSE)
}
metadata <- metadata[order(metadata$group, metadata$subject_id), , drop = FALSE]
sample_columns <- metadata$sample_id
rownames(metadata) <- metadata$sample_id
sample_display <- setNames(
  sub("^CON(?=[-_])", "Control", sample_columns, perl = TRUE),
  sample_columns
)

count_matrix <- as.matrix(counts_input[, sample_columns, drop = FALSE])
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(count_matrix < 0) || any(count_matrix %% 1 != 0)) {
  stop("All expression values must be non-negative integer counts.", call. = FALSE)
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts_input$GeneID

gene_annotation <- counts_input[, c("GeneID", "GeneName"), drop = FALSE]
gene_annotation$GeneName[
  is.na(gene_annotation$GeneName) | gene_annotation$GeneName == ""
] <- gene_annotation$GeneID[
  is.na(gene_annotation$GeneName) | gene_annotation$GeneName == ""
]
rownames(gene_annotation) <- gene_annotation$GeneID

write.csv(
  metadata,
  file.path(tables_dir, "sample_metadata_used.csv"),
  row.names = FALSE,
  quote = TRUE
)

message(
  "Analyzing all ", nrow(count_matrix),
  " input genes without explicit low-count prefiltering."
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ group
)
dds <- DESeq(dds, quiet = TRUE)

raw_results <- results(
  dds,
  contrast = c("group", "Asphyxia", "Control"),
  alpha = alpha,
  independentFiltering = FALSE
)

coefficient_name <- "group_Asphyxia_vs_Control"
if (!coefficient_name %in% resultsNames(dds)) {
  stop(
    "Expected DESeq2 coefficient was not found. Available coefficients: ",
    paste(resultsNames(dds), collapse = ", "),
    call. = FALSE
  )
}

shrunk_results <- lfcShrink(
  dds,
  coef = coefficient_name,
  type = "apeglm",
  quiet = TRUE
)

result_table <- data.frame(
  GeneID = rownames(raw_results),
  GeneName = gene_annotation[rownames(raw_results), "GeneName"],
  baseMean = raw_results$baseMean,
  log2FoldChange_raw = raw_results$log2FoldChange,
  log2FoldChange = shrunk_results$log2FoldChange,
  lfcSE = shrunk_results$lfcSE,
  stat = raw_results$stat,
  pvalue = raw_results$pvalue,
  padj = raw_results$padj,
  stringsAsFactors = FALSE
)

result_table$FDRRegulation <- "Not significant"
result_table$FDRRegulation[
  !is.na(result_table$padj) &
    result_table$padj < alpha &
    result_table$log2FoldChange > 0
] <- "FDR up"
result_table$FDRRegulation[
  !is.na(result_table$padj) &
    result_table$padj < alpha &
    result_table$log2FoldChange < 0
] <- "FDR down"

result_table$ThresholdRegulation <- "Not significant"
result_table$ThresholdRegulation[
  !is.na(result_table$padj) &
    result_table$padj < alpha &
    result_table$log2FoldChange >= lfc_threshold
] <- "Up"
result_table$ThresholdRegulation[
  !is.na(result_table$padj) &
    result_table$padj < alpha &
    result_table$log2FoldChange <= -lfc_threshold
] <- "Down"

result_table <- result_table[
  order(is.na(result_table$padj), result_table$padj, -abs(result_table$log2FoldChange)),
  ,
  drop = FALSE
]
fdr_table <- result_table[result_table$FDRRegulation != "Not significant", , drop = FALSE]
threshold_table <- result_table[
  result_table$ThresholdRegulation != "Not significant",
  ,
  drop = FALSE
]

write.csv(
  result_table,
  file.path(tables_dir, "deseq2_all_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  fdr_table,
  file.path(tables_dir, "deseq2_fdr_significant_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  threshold_table,
  file.path(tables_dir, "deseq2_effect_size_filtered_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)

normalized_counts <- counts(dds, normalized = TRUE)
normalized_output <- data.frame(
  GeneID = rownames(normalized_counts),
  GeneName = gene_annotation[rownames(normalized_counts), "GeneName"],
  normalized_counts,
  check.names = FALSE
)
write.csv(
  normalized_output,
  file.path(tables_dir, "normalized_counts.csv"),
  row.names = FALSE,
  quote = TRUE
)

vst_object <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_object)
vst_output <- data.frame(
  GeneID = rownames(vst_matrix),
  GeneName = gene_annotation[rownames(vst_matrix), "GeneName"],
  vst_matrix,
  check.names = FALSE
)
write.csv(
  vst_output,
  file.path(tables_dir, "vst_expression_filtered.csv"),
  row.names = FALSE,
  quote = TRUE
)

group_colors <- c(Control = "#4C78A8", Asphyxia = "#E45756")
threshold_colors <- c(
  Down = "#4C78A8",
  `Not significant` = "#C9CED6",
  Up = "#E45756"
)

theme_project <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = element_text(color = "#4B5563", margin = margin(b = 10)),
      axis.title = element_text(face = "bold", color = "#273142"),
      axis.text = element_text(color = "#374151"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E8EBEF", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(12, 16, 12, 12)
    )
}

save_ggplot <- function(plot_object, filename, width, height) {
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    paste0(filename, ".png"),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

# Library-size plot.
library_data <- data.frame(
  Sample = sample_columns,
  Reads = colSums(count_matrix[, sample_columns, drop = FALSE]),
  Group = metadata[sample_columns, "group"]
)
library_data$Sample <- factor(
  library_data$Sample,
  levels = sample_columns,
  labels = unname(sample_display[sample_columns])
)

library_plot <- ggplot(library_data, aes(x = Sample, y = Reads / 1e6, fill = Group)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.1f M", Reads / 1e6)),
    vjust = -0.45,
    size = 3.4,
    fontface = "bold",
    color = "#374151"
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Sequencing library sizes",
    subtitle = "R1 frontal cortex transcriptomes before normalization",
    x = NULL,
    y = "Total counts (millions)",
    fill = "Group"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_ggplot(library_plot, file.path(qc_dir, "01_library_sizes"), 8.2, 5.4)

# VST expression distribution.
distribution_data <- data.frame(
  Sample = rep(colnames(vst_matrix), each = nrow(vst_matrix)),
  Expression = as.vector(vst_matrix),
  stringsAsFactors = FALSE
)
distribution_data$Group <- metadata[distribution_data$Sample, "group"]
distribution_data$Sample <- factor(
  distribution_data$Sample,
  levels = sample_columns,
  labels = unname(sample_display[sample_columns])
)

distribution_plot <- ggplot(
  distribution_data,
  aes(x = Sample, y = Expression, fill = Group)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.45,
    alpha = 0.88
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Normalized expression distributions",
    subtitle = "Variance-stabilizing transformation without gene prefiltering",
    x = NULL,
    y = "VST expression",
    fill = "Group"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_ggplot(distribution_plot, file.path(qc_dir, "02_vst_distributions"), 8.2, 5.4)

# PCA.
pca_data <- plotPCA(vst_object, intgroup = "group", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data$Sample <- rownames(pca_data)
pca_data$DisplaySample <- unname(sample_display[pca_data$Sample])

pca_hulls <- do.call(
  rbind,
  lapply(split(pca_data, pca_data$group), function(group_data) {
    group_data[chull(group_data$PC1, group_data$PC2), , drop = FALSE]
  })
)

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = group, label = DisplaySample)
) +
  geom_polygon(
    data = pca_hulls,
    aes(x = PC1, y = PC2, fill = group, group = group),
    inherit.aes = FALSE,
    alpha = 0.07,
    color = NA,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "#D9DDE3", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#D9DDE3", linewidth = 0.35) +
  geom_point(size = 4.1, alpha = 0.95) +
  geom_text_repel(
    size = 3.5,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    seed = 20260920,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(
    title = "Principal component analysis",
    subtitle = "R1 frontal cortex transcriptomes",
    x = paste0("PC1 (", percent_variance[[1]], "%)"),
    y = paste0("PC2 (", percent_variance[[2]], "%)"),
    color = "Group"
  ) +
  theme_project() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(clip = "off")
save_ggplot(pca_plot, file.path(qc_dir, "03_pca"), 7.2, 6.0)

write.csv(
  pca_data[, c("Sample", "group", "PC1", "PC2")],
  file.path(tables_dir, "pca_scores.csv"),
  row.names = FALSE,
  quote = TRUE
)

# Sample correlation heatmap.
sample_cor <- cor(vst_matrix, method = "pearson")
annotation_col <- data.frame(Group = metadata[colnames(sample_cor), "group"])
rownames(annotation_col) <- colnames(sample_cor)
annotation_colors <- list(Group = group_colors)

save_pheatmap <- function(filename, width, height, ...) {
  cairo_pdf(paste0(filename, ".pdf"), width = width, height = height)
  pheatmap(...)
  dev.off()
  png(
    paste0(filename, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 320,
    type = "cairo",
    bg = "white"
  )
  pheatmap(...)
  dev.off()
}

cor_breaks <- seq(min(sample_cor), 1, length.out = 101)
cor_palette <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)
save_pheatmap(
  file.path(qc_dir, "04_sample_correlation"),
  width = 7.4,
  height = 6.6,
  mat = sample_cor,
  color = cor_palette,
  breaks = cor_breaks,
  border_color = "white",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  labels_row = unname(sample_display[rownames(sample_cor)]),
  labels_col = unname(sample_display[colnames(sample_cor)]),
  display_numbers = TRUE,
  number_format = "%.3f",
  fontsize = 10,
  fontsize_number = 8,
  main = "Sample correlation | VST expression"
)

# Volcano plot.
volcano_data <- result_table
volcano_data$minus_log10_padj <- -log10(volcano_data$padj)
finite_y <- volcano_data$minus_log10_padj[is.finite(volcano_data$minus_log10_padj)]
replacement_y <- if (length(finite_y) > 0) max(finite_y) + 0.5 else 1
volcano_data$minus_log10_padj[is.infinite(volcano_data$minus_log10_padj)] <- replacement_y
volcano_data$minus_log10_padj[is.na(volcano_data$minus_log10_padj)] <- 0
volcano_x_limit <- 3
volcano_data$log2FoldChange_plot <- pmax(
  pmin(volcano_data$log2FoldChange, volcano_x_limit),
  -volcano_x_limit
)

label_candidates <- volcano_data[
  volcano_data$ThresholdRegulation != "Not significant",
  ,
  drop = FALSE
]
label_candidates <- head(label_candidates[order(label_candidates$padj), , drop = FALSE], top_labels)

volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FoldChange_plot, y = minus_log10_padj, color = ThresholdRegulation)
) +
  geom_point(alpha = 0.70, size = 1.65) +
  geom_vline(
    xintercept = c(-lfc_threshold, lfc_threshold),
    linetype = "dashed",
    color = "#68707C",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    color = "#68707C",
    linewidth = 0.45
  ) +
  geom_text_repel(
    data = label_candidates,
    aes(label = GeneName),
    size = 3.15,
    box.padding = 0.42,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 20260920,
    show.legend = FALSE
  ) +
  scale_color_manual(values = threshold_colors) +
  labs(
    title = "Differential expression in the frontal cortex",
    subtitle = "asphyxia vs control",
    x = "log2 fold change",
    y = expression(-log[10](adjusted~P)),
    color = "Classification"
  ) +
  theme_project()
save_ggplot(volcano_plot, file.path(de_dir, "01_volcano"), 8.0, 6.4)

# MA plot.
ma_data <- result_table
ma_label_candidates <- ma_data[
  ma_data$ThresholdRegulation != "Not significant",
  ,
  drop = FALSE
]
ma_plot <- ggplot(
  ma_data,
  aes(x = log10(baseMean + 1), y = log2FoldChange, color = ThresholdRegulation)
) +
  geom_point(alpha = 0.67, size = 1.55) +
  geom_hline(yintercept = 0, color = "#68707C", linewidth = 0.45) +
  geom_hline(
    yintercept = c(-lfc_threshold, lfc_threshold),
    linetype = "dashed",
    color = "#A0A7B2",
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = ma_label_candidates,
    aes(label = GeneName),
    size = 3.1,
    box.padding = 0.42,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 20260920,
    show.legend = FALSE
  ) +
  scale_color_manual(values = threshold_colors) +
  labs(
    title = "MA plot",
    subtitle = "Effect size across the expression range",
    x = expression(log[10](baseMean + 1)),
    y = "Shrunken log2 fold change",
    color = "Classification"
  ) +
  theme_project() +
  coord_cartesian(ylim = c(-1.25, 1.25), clip = "on")
save_ggplot(ma_plot, file.path(de_dir, "02_ma_plot"), 8.2, 6.2)

# Differential-expression heatmap. Prefer genes meeting both thresholds.
if (nrow(threshold_table) >= 2) {
  heatmap_candidates <- head(threshold_table, top_heatmap)
  heatmap_title <- paste0(
    "Differentially expressed genes (n = ",
    nrow(heatmap_candidates), ")"
  )
} else {
  heatmap_candidates <- head(fdr_table, min(30L, nrow(fdr_table)))
  heatmap_title <- paste0(
    "Top FDR-significant genes; fewer than two effect-size-filtered hits (n = ",
    nrow(heatmap_candidates), ")"
  )
}

if (nrow(heatmap_candidates) >= 2) {
  heatmap_matrix <- vst_matrix[heatmap_candidates$GeneID, sample_columns, drop = FALSE]
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  heatmap_matrix[heatmap_matrix > 2] <- 2
  heatmap_matrix[heatmap_matrix < -2] <- -2
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  row_labels <- make.unique(
    ifelse(
      is.na(heatmap_candidates$GeneName) | heatmap_candidates$GeneName == "",
      heatmap_candidates$GeneID,
      heatmap_candidates$GeneName
    )
  )
  rownames(heatmap_matrix) <- row_labels

  heatmap_palette <- colorRampPalette(c("#3B6FB6", "#F7F7F7", "#D1495B"))(101)
  save_pheatmap(
    file.path(de_dir, "03_differential_gene_heatmap"),
    width = 7.7,
    height = max(5.6, min(12, 2.7 + 0.24 * nrow(heatmap_matrix))),
    mat = heatmap_matrix,
    color = heatmap_palette,
    breaks = seq(-2, 2, length.out = 102),
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col[sample_columns, , drop = FALSE],
    annotation_colors = annotation_colors,
    labels_col = unname(sample_display[sample_columns]),
    show_colnames = TRUE,
    fontsize = 10,
    fontsize_row = if (nrow(heatmap_matrix) > 35) 7 else 9,
    main = heatmap_title
  )
}

summary_table <- data.frame(
  Metric = c(
    "Input genes",
    "Genes analyzed without prefiltering",
    "Control samples",
    "Asphyxia samples",
    "FDR-significant upregulated genes",
    "FDR-significant downregulated genes",
    "FDR- and effect-size-filtered upregulated genes",
    "FDR- and effect-size-filtered downregulated genes",
    "Adjusted P-value threshold",
    "Absolute shrunken log2FC threshold"
  ),
  Value = c(
    nrow(count_matrix),
    nrow(count_matrix),
    sum(metadata$group == "Control"),
    sum(metadata$group == "Asphyxia"),
    sum(result_table$FDRRegulation == "FDR up"),
    sum(result_table$FDRRegulation == "FDR down"),
    sum(result_table$ThresholdRegulation == "Up"),
    sum(result_table$ThresholdRegulation == "Down"),
    alpha,
    lfc_threshold
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary_table,
  file.path(output_dir, "analysis_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

session_lines <- capture.output(sessionInfo())
session_lines <- sub(
  "^BLAS/LAPACK: .*;[[:space:]]*LAPACK version:",
  "BLAS/LAPACK: <system library>; LAPACK version:",
  session_lines
)
writeLines(session_lines, file.path(output_dir, "session_info.txt"))

message(
  "Analysis complete. FDR-significant genes: ", nrow(fdr_table),
  " (Up = ", sum(fdr_table$FDRRegulation == "FDR up"),
  ", Down = ", sum(fdr_table$FDRRegulation == "FDR down"),
  "). FDR- and effect-size-filtered genes: ", nrow(threshold_table), "."
)
message("Results written to: ", normalizePath(output_dir))
