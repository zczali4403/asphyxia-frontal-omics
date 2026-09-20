#!/usr/bin/env Rscript

# Differential proteomic analysis of the R1 frontal cortex.
# Comparison: Asphyxia vs Control.

required_packages <- c(
  "limma", "readxl", "ggplot2", "ggrepel", "pheatmap", "statmod"
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
  library(limma)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
protein_file <- if (length(args) >= 1) args[[1]] else
  "data/proteomics/protein_expression_log2.xlsx"
metadata_file <- if (length(args) >= 2) args[[2]] else
  "config/samples.csv"
output_dir <- if (length(args) >= 3) args[[3]] else
  "results/02_proteomics"

matrix_sheet <- "数据矩阵"
max_missing_fraction <- 0.50
alpha <- 0.05
lfc_threshold <- 0.5
top_labels <- 12L
top_heatmap <- 30L

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
qc_dir <- file.path(figures_dir, "qc")
de_dir <- file.path(figures_dir, "differential_expression")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading protein matrix: ", protein_file, " [", matrix_sheet, "]")
protein_input <- as.data.frame(
  read_excel(protein_file, sheet = matrix_sheet, .name_repair = "minimal"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Accession", "Gene Name", "Description",
  "Missing_Percent(Asphyxia)", "Missing_Percent(CON)"
)
if (!all(required_columns %in% colnames(protein_input))) {
  stop(
    "The protein matrix is missing required columns: ",
    paste(setdiff(required_columns, colnames(protein_input)), collapse = ", "),
    call. = FALSE
  )
}
if (anyNA(protein_input$Accession) || anyDuplicated(protein_input$Accession)) {
  stop("Accession values must be present and unique.", call. = FALSE)
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
  metadata_all$omics == "proteomics" &
    metadata_all$sample_role == "biological" &
    metadata_all$region == "frontal_cortex",
  ,
  drop = FALSE
]
if (nrow(metadata) == 0 || anyDuplicated(metadata$sample_id)) {
  stop("Valid frontal-cortex proteomic metadata were not found.", call. = FALSE)
}
metadata$group[metadata$group == "CON"] <- "Control"
metadata$group <- factor(metadata$group, levels = c("Control", "Asphyxia"))
if (anyNA(metadata$group)) {
  stop("Proteomic groups must be Control (or legacy CON) and Asphyxia.", call. = FALSE)
}
metadata <- metadata[order(metadata$group, metadata$subject_id), , drop = FALSE]
sample_columns <- metadata$sample_id
rownames(metadata) <- sample_columns

missing_samples <- setdiff(sample_columns, colnames(protein_input))
if (length(missing_samples) > 0) {
  stop(
    "Samples in metadata but absent from the protein matrix: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}

sample_display <- setNames(
  sub("^CON(?=[-_])", "Control", sample_columns, perl = TRUE),
  sample_columns
)

expression_matrix <- as.matrix(protein_input[, sample_columns, drop = FALSE])
suppressWarnings(storage.mode(expression_matrix) <- "numeric")
if (anyNA(expression_matrix) || any(!is.finite(expression_matrix))) {
  stop(
    "The processed protein expression matrix must contain finite log2 values.",
    call. = FALSE
  )
}
rownames(expression_matrix) <- protein_input$Accession

missing_asphyxia <- as.numeric(protein_input[["Missing_Percent(Asphyxia)"]])
missing_control <- as.numeric(protein_input[["Missing_Percent(CON)"]])
if (
  anyNA(missing_asphyxia) || anyNA(missing_control) ||
    any(missing_asphyxia < 0 | missing_asphyxia > 1) ||
    any(missing_control < 0 | missing_control > 1)
) {
  stop("Protein missingness fractions must be numeric values from 0 to 1.", call. = FALSE)
}

# Require detection in at least two of four samples in each group.
keep <- missing_asphyxia <= max_missing_fraction &
  missing_control <= max_missing_fraction
filtered_expression <- expression_matrix[keep, , drop = FALSE]
if (nrow(filtered_expression) < 2) {
  stop("Too few proteins remain after missingness filtering.", call. = FALSE)
}
message(
  "Retained ", nrow(filtered_expression), " of ", nrow(expression_matrix),
  " proteins after missingness filtering."
)

gene_name <- protein_input$`Gene Name`
if ("Gene_Name" %in% colnames(protein_input)) {
  fallback <- is.na(gene_name) | trimws(gene_name) == ""
  gene_name[fallback] <- protein_input$Gene_Name[fallback]
}
fallback <- is.na(gene_name) | trimws(gene_name) == ""
gene_name[fallback] <- protein_input$Accession[fallback]

annotation <- data.frame(
  Accession = protein_input$Accession,
  GeneName = gene_name,
  GeneID = if ("GeneID" %in% colnames(protein_input))
    as.character(protein_input$GeneID) else NA_character_,
  Description = protein_input$Description,
  MissingFractionControl = missing_control,
  MissingFractionAsphyxia = missing_asphyxia,
  stringsAsFactors = FALSE
)
rownames(annotation) <- annotation$Accession

write.csv(
  metadata,
  file.path(tables_dir, "sample_metadata_used.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  annotation[keep, , drop = FALSE],
  file.path(tables_dir, "protein_annotation_filtered.csv"),
  row.names = FALSE,
  quote = TRUE
)

# Limma model: unpaired comparison of four biological replicates per group.
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- levels(metadata$group)
contrast_matrix <- makeContrasts(
  Asphyxia_vs_Control = Asphyxia - Control,
  levels = design
)
fit <- lmFit(filtered_expression, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

limma_table <- topTable(
  fit,
  coef = "Asphyxia_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
limma_table$Accession <- rownames(limma_table)

result_table <- data.frame(
  Accession = limma_table$Accession,
  GeneName = annotation[limma_table$Accession, "GeneName"],
  GeneID = annotation[limma_table$Accession, "GeneID"],
  Description = annotation[limma_table$Accession, "Description"],
  MissingFractionControl = annotation[
    limma_table$Accession, "MissingFractionControl"
  ],
  MissingFractionAsphyxia = annotation[
    limma_table$Accession, "MissingFractionAsphyxia"
  ],
  AveExpr = limma_table$AveExpr,
  log2FoldChange = limma_table$logFC,
  t = limma_table$t,
  pvalue = limma_table$P.Value,
  padj = limma_table$adj.P.Val,
  B = limma_table$B,
  stringsAsFactors = FALSE
)

result_table$FDRRegulation <- "Not significant"
result_table$FDRRegulation[
  result_table$padj < alpha & result_table$log2FoldChange > 0
] <- "FDR up"
result_table$FDRRegulation[
  result_table$padj < alpha & result_table$log2FoldChange < 0
] <- "FDR down"

result_table$ThresholdRegulation <- "Not significant"
result_table$ThresholdRegulation[
  result_table$padj < alpha & result_table$log2FoldChange >= lfc_threshold
] <- "Up"
result_table$ThresholdRegulation[
  result_table$padj < alpha & result_table$log2FoldChange <= -lfc_threshold
] <- "Down"

result_table$ExploratoryRegulation <- "Not significant"
result_table$ExploratoryRegulation[
  result_table$pvalue < alpha & result_table$log2FoldChange >= lfc_threshold
] <- "Up"
result_table$ExploratoryRegulation[
  result_table$pvalue < alpha & result_table$log2FoldChange <= -lfc_threshold
] <- "Down"

fdr_table <- result_table[
  result_table$FDRRegulation != "Not significant", , drop = FALSE
]
threshold_table <- result_table[
  result_table$ThresholdRegulation != "Not significant", , drop = FALSE
]
exploratory_table <- result_table[
  result_table$ExploratoryRegulation != "Not significant", , drop = FALSE
]

write.csv(
  result_table,
  file.path(tables_dir, "limma_all_proteins.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  fdr_table,
  file.path(tables_dir, "limma_fdr_significant_proteins.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  threshold_table,
  file.path(tables_dir, "limma_effect_size_filtered_proteins.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  exploratory_table,
  file.path(tables_dir, "limma_exploratory_nominal_p_proteins.csv"),
  row.names = FALSE,
  quote = TRUE
)

expression_output <- data.frame(
  Accession = rownames(filtered_expression),
  GeneName = annotation[rownames(filtered_expression), "GeneName"],
  filtered_expression,
  check.names = FALSE
)
write.csv(
  expression_output,
  file.path(tables_dir, "protein_expression_filtered.csv"),
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
    paste0(filename, ".pdf"), plot = plot_object,
    width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    paste0(filename, ".png"), plot = plot_object,
    width = width, height = height, units = "in",
    dpi = 320, bg = "white"
  )
}

save_pheatmap <- function(filename, width, height, ...) {
  cairo_pdf(paste0(filename, ".pdf"), width = width, height = height)
  pheatmap(...)
  dev.off()
  png(
    paste0(filename, ".png"), width = width, height = height,
    units = "in", res = 320, type = "cairo", bg = "white"
  )
  pheatmap(...)
  dev.off()
}

# Original detection completeness.
detection_data <- rbind(
  data.frame(
    Group = "Control",
    DetectedSamples = round(4 * (1 - missing_control))
  ),
  data.frame(
    Group = "Asphyxia",
    DetectedSamples = round(4 * (1 - missing_asphyxia))
  )
)
detection_data$Group <- factor(
  detection_data$Group,
  levels = c("Control", "Asphyxia")
)
detection_summary <- aggregate(
  list(Proteins = rep(1L, nrow(detection_data))),
  by = list(
    Group = detection_data$Group,
    DetectedSamples = detection_data$DetectedSamples
  ),
  FUN = sum
)
detection_plot <- ggplot(
  detection_summary,
  aes(x = factor(DetectedSamples), y = Proteins, fill = Group)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Protein detection completeness",
    subtitle = "Number of original observations before missing-value filling",
    x = "Samples with detected protein (out of 4)",
    y = "Number of proteins",
    fill = "Group"
  ) +
  theme_project()
save_ggplot(detection_plot, file.path(qc_dir, "01_detection_completeness"), 8.0, 5.5)

# Expression distributions.
distribution_data <- data.frame(
  Sample = rep(sample_columns, each = nrow(filtered_expression)),
  Expression = as.vector(filtered_expression),
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
  geom_boxplot(width = 0.64, outlier.shape = NA, linewidth = 0.45, alpha = 0.88) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Protein expression distributions",
    subtitle = "Filtered log2 protein expression matrix",
    x = NULL,
    y = "log2 expression",
    fill = "Group"
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_ggplot(distribution_plot, file.path(qc_dir, "02_expression_distributions"), 9.0, 5.6)

# PCA.
pca <- prcomp(t(filtered_expression), center = TRUE, scale. = FALSE)
percent_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_data <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
pca_data$Group <- metadata[pca_data$Sample, "group"]
pca_data$DisplaySample <- unname(sample_display[pca_data$Sample])
pca_hulls <- do.call(
  rbind,
  lapply(split(pca_data, pca_data$Group), function(group_data) {
    group_data[chull(group_data$PC1, group_data$PC2), , drop = FALSE]
  })
)
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Group, label = DisplaySample)
) +
  geom_polygon(
    data = pca_hulls,
    aes(x = PC1, y = PC2, fill = Group, group = Group),
    inherit.aes = FALSE, alpha = 0.07, color = NA, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "#D9DDE3", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#D9DDE3", linewidth = 0.35) +
  geom_point(size = 4.1, alpha = 0.95) +
  geom_text_repel(
    size = 3.4, box.padding = 0.45, point.padding = 0.35,
    min.segment.length = 0, seed = 20260920, show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(
    title = "Principal component analysis",
    subtitle = "R1 frontal cortex proteomes",
    x = paste0("PC1 (", percent_variance[[1]], "%)"),
    y = paste0("PC2 (", percent_variance[[2]], "%)"),
    color = "Group"
  ) +
  theme_project() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(clip = "off")
save_ggplot(pca_plot, file.path(qc_dir, "03_pca"), 7.5, 6.2)
write.csv(
  pca_data[, c("Sample", "Group", "PC1", "PC2")],
  file.path(tables_dir, "pca_scores.csv"),
  row.names = FALSE,
  quote = TRUE
)

# Sample correlation heatmap.
sample_cor <- cor(filtered_expression, method = "pearson")
annotation_col <- data.frame(Group = metadata[colnames(sample_cor), "group"])
rownames(annotation_col) <- colnames(sample_cor)
annotation_colors <- list(Group = group_colors)
cor_breaks <- seq(min(sample_cor), 1, length.out = 101)
cor_palette <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)
save_pheatmap(
  file.path(qc_dir, "04_sample_correlation"),
  width = 8.0,
  height = 7.0,
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
  fontsize = 9,
  fontsize_number = 7,
  main = "Sample correlation | log2 protein expression"
)

make_volcano <- function(data, regulation_column, p_column, subtitle, filename) {
  plot_data <- data
  plot_data$PlotRegulation <- plot_data[[regulation_column]]
  plot_data$minus_log10_p <- -log10(plot_data[[p_column]])
  finite_y <- plot_data$minus_log10_p[is.finite(plot_data$minus_log10_p)]
  replacement_y <- if (length(finite_y) > 0) max(finite_y) + 0.5 else 1
  plot_data$minus_log10_p[is.infinite(plot_data$minus_log10_p)] <- replacement_y
  plot_data$minus_log10_p[is.na(plot_data$minus_log10_p)] <- 0

  labels <- plot_data[
    plot_data$PlotRegulation != "Not significant", , drop = FALSE
  ]
  labels <- head(labels[order(labels[[p_column]]), , drop = FALSE], top_labels)

  plot_object <- ggplot(
    plot_data,
    aes(x = log2FoldChange, y = minus_log10_p, color = PlotRegulation)
  ) +
    geom_point(alpha = 0.70, size = 1.65) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed", color = "#68707C", linewidth = 0.45
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed", color = "#68707C", linewidth = 0.45
    ) +
    geom_text_repel(
      data = labels,
      aes(label = GeneName),
      size = 3.1, box.padding = 0.42, point.padding = 0.25,
      max.overlaps = Inf, min.segment.length = 0,
      seed = 20260920, show.legend = FALSE
    ) +
    scale_color_manual(values = threshold_colors) +
    labs(
      title = "Differential protein expression in the frontal cortex",
      subtitle = subtitle,
      x = "log2 fold change",
      y = if (p_column == "padj")
        expression(-log[10](adjusted~P)) else expression(-log[10](P)),
      color = "Classification"
    ) +
    theme_project()
  save_ggplot(plot_object, file.path(de_dir, filename), 8.0, 6.4)
}

make_volcano(
  result_table,
  "ExploratoryRegulation",
  "pvalue",
  "asphyxia vs control",
  "01_exploratory_volcano"
)

# Exploratory MA plot. FDR-controlled results remain in the exported tables.
ma_labels <- result_table[
  result_table$ExploratoryRegulation != "Not significant", , drop = FALSE
]
ma_labels <- head(ma_labels[order(ma_labels$pvalue), , drop = FALSE], top_labels)
ma_plot <- ggplot(
  result_table,
  aes(x = AveExpr, y = log2FoldChange, color = ExploratoryRegulation)
) +
  geom_point(alpha = 0.67, size = 1.55) +
  geom_hline(yintercept = 0, color = "#68707C", linewidth = 0.45) +
  geom_hline(
    yintercept = c(-lfc_threshold, lfc_threshold),
    linetype = "dashed", color = "#A0A7B2", linewidth = 0.4
  ) +
  geom_text_repel(
    data = ma_labels,
    aes(label = GeneName),
    size = 3.1, box.padding = 0.42, point.padding = 0.25,
    max.overlaps = Inf, min.segment.length = 0,
    seed = 20260920, show.legend = FALSE
  ) +
  scale_color_manual(values = threshold_colors) +
  labs(
    title = "MA plot",
    subtitle = "asphyxia vs control",
    x = "Average log2 expression",
    y = "log2 fold change",
    color = "Classification"
  ) +
  theme_project()
save_ggplot(ma_plot, file.path(de_dir, "02_exploratory_ma_plot"), 8.2, 6.0)

# Heatmap: primary proteins when available, otherwise clearly labeled exploratory hits.
if (nrow(threshold_table) >= 2) {
  heatmap_candidates <- head(threshold_table, top_heatmap)
  heatmap_title <- paste0(
    "Differentially expressed proteins (n = ", nrow(heatmap_candidates), ")"
  )
} else if (nrow(exploratory_table) >= 2) {
  heatmap_candidates <- head(exploratory_table, top_heatmap)
  heatmap_title <- paste0(
    "Differentially expressed proteins (n = ",
    nrow(heatmap_candidates), ")"
  )
} else {
  heatmap_candidates <- head(result_table, min(top_heatmap, nrow(result_table)))
  heatmap_title <- paste0(
    "Top-ranked proteins (n = ", nrow(heatmap_candidates), ")"
  )
}

heatmap_matrix <- filtered_expression[
  heatmap_candidates$Accession, sample_columns, drop = FALSE
]
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix[heatmap_matrix > 2] <- 2
heatmap_matrix[heatmap_matrix < -2] <- -2
heatmap_matrix[is.na(heatmap_matrix)] <- 0
rownames(heatmap_matrix) <- make.unique(heatmap_candidates$GeneName)
heatmap_palette <- colorRampPalette(c("#3B6FB6", "#F7F7F7", "#D1495B"))(101)
save_pheatmap(
  file.path(de_dir, "03_exploratory_protein_heatmap"),
  width = 8.2,
  height = max(5.8, min(10, 3 + 0.24 * nrow(heatmap_matrix))),
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
  fontsize_row = if (nrow(heatmap_matrix) > 25) 7.5 else 9,
  main = heatmap_title
)

summary_table <- data.frame(
  Metric = c(
    "Input proteins",
    "Proteins retained after missingness filtering",
    "Control samples",
    "Asphyxia samples",
    "FDR-significant proteins",
    "FDR- and effect-size-filtered upregulated proteins",
    "FDR- and effect-size-filtered downregulated proteins",
    "Exploratory nominal-P upregulated proteins",
    "Exploratory nominal-P downregulated proteins",
    "Adjusted P-value threshold",
    "Nominal P-value threshold for exploratory results",
    "Absolute log2FC threshold",
    "Maximum missing fraction per group"
  ),
  Value = c(
    nrow(expression_matrix),
    nrow(filtered_expression),
    sum(metadata$group == "Control"),
    sum(metadata$group == "Asphyxia"),
    nrow(fdr_table),
    sum(result_table$ThresholdRegulation == "Up"),
    sum(result_table$ThresholdRegulation == "Down"),
    sum(result_table$ExploratoryRegulation == "Up"),
    sum(result_table$ExploratoryRegulation == "Down"),
    alpha,
    alpha,
    lfc_threshold,
    max_missing_fraction
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
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

message(
  "Analysis complete. FDR-significant proteins: ", nrow(fdr_table),
  "; FDR- and effect-size-filtered proteins: ", nrow(threshold_table),
  "; exploratory nominal-P proteins: ", nrow(exploratory_table), "."
)
message("Results written to: ", normalizePath(output_dir))
