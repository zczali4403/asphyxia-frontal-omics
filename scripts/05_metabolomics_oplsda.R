#!/usr/bin/env Rscript

# Frontal-cortex metabolomics analysis with OPLS-DA, VIP, permutation
# validation, and limma statistics for Asphyxia versus Control.

required_packages <- c(
  "readxl", "ggplot2", "ggrepel", "limma", "pheatmap", "ropls", "scales"
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
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(pheatmap)
  library(ropls)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else
  "data/metabolomics/metabolite_expression.xlsx"
metadata_file <- if (length(args) >= 2) args[[2]] else "config/samples.csv"
output_dir <- if (length(args) >= 3) args[[3]] else "results/05_metabolomics"

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
models_dir <- file.path(output_dir, "models")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

processed_sheet <- "数据矩阵"
original_sheet <- "缺失值数据矩阵"
detection_threshold <- 0.50
qc_rsd_threshold <- 30
vip_threshold <- 1
pvalue_threshold <- 0.05
fdr_threshold <- 0.05
logfc_threshold <- 0.5
permutation_count <- 200L
cross_validation_folds <- 6L
top_vip_count <- 30L
set.seed(20260921)

group_colors <- c(Control = "#4C78A8", Asphyxia = "#E45756", QC = "#59A14F")
direction_colors <- c(
  Down = "#4C78A8", `Not selected` = "#C8CED6", Up = "#E45756"
)

theme_project <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold", size = base_size + 3, hjust = 0.5, color = "#111827"
      ),
      plot.subtitle = element_text(
        color = "#4B5563", hjust = 0.5, margin = margin(b = 10)
      ),
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

message("Reading metabolomics matrices and metadata")
processed <- read_excel(input_file, sheet = processed_sheet)
original <- read_excel(input_file, sheet = original_sheet)
metadata <- read.csv(
  metadata_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_annotation <- c(
  "ID", "Metabolites", "Metabolites_cn", "KEGG", "HMDB",
  "Super Class", "Class", "Sub Class", "level", "Score"
)
if (!all(required_annotation %in% colnames(processed))) {
  stop("The metabolomics worksheet is missing required annotation columns.")
}
metabolomics_metadata <- metadata[
  metadata$omics == "metabolomics" &
    metadata$region == "frontal_cortex",
  ,
  drop = FALSE
]
biological_metadata <- metabolomics_metadata[
  metabolomics_metadata$sample_role == "biological",
  ,
  drop = FALSE
]
qc_metadata <- metabolomics_metadata[
  metabolomics_metadata$sample_role == "quality_control",
  ,
  drop = FALSE
]
biological_metadata$group <- factor(
  biological_metadata$group,
  levels = c("Control", "Asphyxia")
)
biological_samples <- biological_metadata$sample_id
qc_samples <- qc_metadata$sample_id
all_samples <- c(biological_samples, qc_samples)

if (!all(all_samples %in% colnames(processed)) ||
    !all(all_samples %in% colnames(original))) {
  stop("Metabolomics sample names do not match the metadata.")
}
if (!identical(as.character(processed$ID), as.character(original$ID))) {
  stop("Feature IDs differ between the processed and original worksheets.")
}

processed_matrix <- as.matrix(processed[, all_samples, drop = FALSE])
original_matrix <- as.matrix(original[, all_samples, drop = FALSE])
storage.mode(processed_matrix) <- "numeric"
storage.mode(original_matrix) <- "numeric"
rownames(processed_matrix) <- processed$ID
rownames(original_matrix) <- original$ID
if (anyNA(processed_matrix)) {
  stop("The processed metabolomics matrix contains missing values.")
}

control_samples <- biological_metadata$sample_id[
  biological_metadata$group == "Control"
]
asphyxia_samples <- biological_metadata$sample_id[
  biological_metadata$group == "Asphyxia"
]
control_detection <- rowMeans(original_matrix[, control_samples, drop = FALSE] > 0)
asphyxia_detection <- rowMeans(original_matrix[, asphyxia_samples, drop = FALSE] > 0)
qc_detection <- rowMeans(original_matrix[, qc_samples, drop = FALSE] > 0)
qc_mean <- rowMeans(original_matrix[, qc_samples, drop = FALSE])
qc_sd <- apply(original_matrix[, qc_samples, drop = FALSE], 1, sd)
qc_rsd <- 100 * qc_sd / qc_mean
processed_variance <- apply(
  processed_matrix[, biological_samples, drop = FALSE], 1, var
)

keep_detection <- control_detection >= detection_threshold |
  asphyxia_detection >= detection_threshold
keep_qc <- qc_detection == 1 & is.finite(qc_rsd) &
  qc_rsd <= qc_rsd_threshold
keep_variance <- is.finite(processed_variance) & processed_variance > 0
keep_feature <- keep_detection & keep_qc & keep_variance

qc_table <- data.frame(
  ID = processed$ID,
  Metabolites = processed$Metabolites,
  Metabolites_cn = processed$Metabolites_cn,
  ControlDetectionRate = control_detection,
  AsphyxiaDetectionRate = asphyxia_detection,
  QCDetectionRate = qc_detection,
  QCMean = qc_mean,
  QCRSD = qc_rsd,
  ProcessedVariance = processed_variance,
  PassDetection = keep_detection,
  PassQCRSD = keep_qc,
  PassVariance = keep_variance,
  Retained = keep_feature,
  stringsAsFactors = FALSE
)
write.csv(
  qc_table, file.path(tables_dir, "metabolite_qc_metrics.csv"),
  row.names = FALSE, quote = TRUE
)

processed_filtered <- processed[keep_feature, , drop = FALSE]
expression_filtered <- processed_matrix[
  keep_feature, biological_samples, drop = FALSE
]
feature_keys <- make.unique(as.character(processed_filtered$ID))
rownames(expression_filtered) <- feature_keys
annotation_columns <- setdiff(colnames(processed_filtered), all_samples)
write.csv(
  processed_filtered[, c(annotation_columns, all_samples), drop = FALSE],
  file.path(tables_dir, "metabolite_expression_filtered.csv"),
  row.names = FALSE, quote = TRUE
)

qc_plot_data <- qc_table[is.finite(qc_table$QCRSD), , drop = FALSE]
qc_rsd_plot <- ggplot(qc_plot_data, aes(x = QCRSD)) +
  geom_histogram(
    bins = 35, fill = "#59A14F", color = "white", linewidth = 0.35,
    alpha = 0.88
  ) +
  geom_vline(
    xintercept = qc_rsd_threshold,
    color = "#E45756", linetype = "dashed", linewidth = 0.8
  ) +
  annotate(
    "text", x = qc_rsd_threshold, y = Inf,
    label = paste0(qc_rsd_threshold, "% cutoff"),
    hjust = 1.08, vjust = 1.5,
    color = "#E45756", fontface = "bold", size = 3.6
  ) +
  labs(
    title = "Metabolite QC reproducibility",
    x = "QC relative standard deviation (%)",
    y = "Number of metabolite features"
  ) +
  theme_project()
save_ggplot(
  qc_rsd_plot, file.path(figures_dir, "01_qc_rsd_distribution"), 8.2, 5.8
)

pareto_scale <- function(sample_by_feature_matrix) {
  centered <- sweep(
    sample_by_feature_matrix, 2, colMeans(sample_by_feature_matrix), "-"
  )
  feature_sd <- apply(sample_by_feature_matrix, 2, sd)
  sweep(centered, 2, sqrt(feature_sd), "/")
}

all_expression <- t(processed_matrix[keep_feature, all_samples, drop = FALSE])
pca_scaled <- pareto_scale(all_expression)
pca_model <- prcomp(pca_scaled, center = FALSE, scale. = FALSE)
pca_variance <- 100 * pca_model$sdev^2 / sum(pca_model$sdev^2)
pca_scores <- data.frame(
  Sample = rownames(pca_model$x),
  PC1 = pca_model$x[, 1],
  PC2 = pca_model$x[, 2],
  stringsAsFactors = FALSE
)
pca_scores <- merge(
  pca_scores,
  metabolomics_metadata[, c("sample_id", "group", "sample_role")],
  by.x = "Sample", by.y = "sample_id", sort = FALSE
)
pca_scores$DisplayGroup <- ifelse(
  pca_scores$sample_role == "quality_control", "QC", pca_scores$group
)
pca_scores$DisplayGroup <- factor(
  pca_scores$DisplayGroup, levels = c("Control", "Asphyxia", "QC")
)
write.csv(
  pca_scores, file.path(tables_dir, "pca_scores.csv"),
  row.names = FALSE, quote = TRUE
)

pca_plot <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = DisplayGroup, shape = DisplayGroup)
) +
  stat_ellipse(
    data = pca_scores[pca_scores$DisplayGroup != "QC", , drop = FALSE],
    aes(group = DisplayGroup, color = DisplayGroup),
    type = "norm", level = 0.80, linewidth = 0.75,
    show.legend = FALSE
  ) +
  geom_point(size = 4.2, alpha = 0.95) +
  geom_text_repel(
    aes(label = Sample), size = 3.2, show.legend = FALSE,
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
    seed = 20260921
  ) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(Control = 16, Asphyxia = 16, QC = 17)) +
  labs(
    title = "Metabolomics PCA",
    x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
    y = sprintf("PC2 (%.1f%%)", pca_variance[2]),
    color = "Group", shape = "Group"
  ) +
  theme_project()
save_ggplot(pca_plot, file.path(figures_dir, "02_pca"), 8.4, 6.6)

message("Running limma differential metabolite analysis")
group <- biological_metadata$group
design <- model.matrix(~group)
fit <- eBayes(lmFit(expression_filtered, design))
limma_table <- topTable(
  fit, coef = "groupAsphyxia", number = Inf, sort.by = "none"
)
limma_table <- limma_table[feature_keys, , drop = FALSE]

message("Fitting OPLS-DA model with permutation validation")
x_model <- t(expression_filtered)
colnames(x_model) <- feature_keys
rownames(x_model) <- biological_samples
set.seed(20260921)
opls_model <- opls(
  x_model,
  group,
  predI = 1,
  orthoI = 1,
  algoC = "nipals",
  crossvalI = cross_validation_folds,
  permI = permutation_count,
  scaleC = "pareto",
  fig.pdfC = "none",
  info.txtC = "none"
)
saveRDS(opls_model, file.path(models_dir, "metabolomics_oplsda_model.rds"))

vip <- getVipVn(opls_model)[feature_keys]
predictive_loading <- as.numeric(getLoadingMN(opls_model)[feature_keys, 1])
predictive_score <- as.numeric(getScoreMN(opls_model)[, 1])
orthogonal_score <- as.numeric(opls_model@orthoScoreMN[, 1])
s_correlation <- vapply(
  seq_len(ncol(x_model)),
  function(index) cor(x_model[, index], predictive_score),
  numeric(1)
)

result_table <- processed_filtered[, annotation_columns, drop = FALSE]
result_table$FeatureKey <- feature_keys
result_table$log2FC <- limma_table$logFC
result_table$AverageExpression <- limma_table$AveExpr
result_table$t <- limma_table$t
result_table$Pvalue <- limma_table$P.Value
result_table$padj <- limma_table$adj.P.Val
result_table$B <- limma_table$B
result_table$VIP <- as.numeric(vip)
result_table$PredictiveLoading <- predictive_loading
result_table$SCorrelation <- s_correlation
result_table$FDRSignificant <- result_table$padj < fdr_threshold &
  abs(result_table$log2FC) >= logfc_threshold
result_table$VIPCandidate <- result_table$VIP > vip_threshold &
  result_table$Pvalue < pvalue_threshold &
  abs(result_table$log2FC) >= logfc_threshold
result_table$Direction <- ifelse(
  result_table$VIPCandidate & result_table$log2FC > 0,
  "Up",
  ifelse(
    result_table$VIPCandidate & result_table$log2FC < 0,
    "Down", "Not selected"
  )
)
result_table <- result_table[
  order(!result_table$VIPCandidate, -result_table$VIP, result_table$Pvalue),
  ,
  drop = FALSE
]
write.csv(
  result_table, file.path(tables_dir, "oplsda_vip_all_metabolites.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  result_table[result_table$VIPCandidate, , drop = FALSE],
  file.path(tables_dir, "oplsda_vip_candidate_metabolites.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  result_table[result_table$FDRSignificant, , drop = FALSE],
  file.path(tables_dir, "fdr_significant_metabolites.csv"),
  row.names = FALSE, quote = TRUE
)

predicted_group <- fitted(opls_model)
score_table <- data.frame(
  Sample = biological_samples,
  Group = group,
  PredictiveScore = predictive_score,
  OrthogonalScore = orthogonal_score,
  PredictedGroup = as.character(predicted_group[biological_samples]),
  stringsAsFactors = FALSE
)
write.csv(
  score_table, file.path(tables_dir, "oplsda_scores.csv"),
  row.names = FALSE, quote = TRUE
)

score_plot <- ggplot(
  score_table,
  aes(x = PredictiveScore, y = OrthogonalScore, color = Group)
) +
  stat_ellipse(
    aes(fill = Group), geom = "polygon", type = "norm", level = 0.80,
    alpha = 0.12, linewidth = 0.7, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "#A7AFBA", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#A7AFBA", linewidth = 0.4) +
  geom_point(size = 4.6, alpha = 0.96) +
  geom_text_repel(
    aes(label = Sample), size = 3.3, show.legend = FALSE,
    max.overlaps = Inf, box.padding = 0.45, point.padding = 0.3,
    seed = 20260921
  ) +
  scale_color_manual(values = group_colors[c("Control", "Asphyxia")]) +
  scale_fill_manual(values = group_colors[c("Control", "Asphyxia")]) +
  labs(
    title = "OPLS-DA score plot",
    x = "Predictive component score",
    y = "Orthogonal component score",
    color = "Group"
  ) +
  theme_project()
save_ggplot(
  score_plot, file.path(figures_dir, "03_oplsda_score_plot"), 8.4, 6.6
)

model_summary <- getSummaryDF(opls_model)
permutation_matrix <- as.data.frame(opls_model@suppLs$permMN)
permutation_matrix$Permutation <- seq_len(nrow(permutation_matrix)) - 1L
permutation_matrix$Observed <- permutation_matrix$sim == 1
write.csv(
  permutation_matrix,
  file.path(tables_dir, "oplsda_permutation_results.csv"),
  row.names = FALSE, quote = TRUE
)

permutation_long <- rbind(
  data.frame(
    Similarity = permutation_matrix$sim,
    Metric = "R2Y",
    Value = permutation_matrix[["R2Y(cum)"]],
    Observed = permutation_matrix$Observed
  ),
  data.frame(
    Similarity = permutation_matrix$sim,
    Metric = "Q2",
    Value = permutation_matrix[["Q2(cum)"]],
    Observed = permutation_matrix$Observed
  )
)
permutation_long$Metric <- factor(
  permutation_long$Metric, levels = c("R2Y", "Q2")
)
validation_plot <- ggplot(
  permutation_long[!permutation_long$Observed, , drop = FALSE],
  aes(x = Similarity, y = Value, color = Metric)
) +
  geom_hline(yintercept = 0, color = "#A7AFBA", linewidth = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2.0, alpha = 0.60) +
  geom_point(
    data = permutation_long[permutation_long$Observed, , drop = FALSE],
    aes(fill = Metric), shape = 23, size = 4.4,
    color = "#273142", stroke = 0.8
  ) +
  annotate(
    "label", x = 0.03, y = Inf,
    label = sprintf(
      "R2Y = %.2f  |  Q2 = %.2f\np(R2Y) = %.3g  |  p(Q2) = %.3g",
      model_summary$`R2Y(cum)`, model_summary$`Q2(cum)`,
      model_summary$pR2Y, model_summary$pQ2
    ),
    hjust = 0, vjust = 1.3, size = 3.5,
    color = "#273142", fill = "white", label.size = 0.25
  ) +
  scale_color_manual(values = c(R2Y = "#E45756", Q2 = "#4C78A8")) +
  scale_fill_manual(values = c(R2Y = "#E45756", Q2 = "#4C78A8")) +
  labs(
    title = "OPLS-DA permutation validation",
    x = "Class-label correlation",
    y = "Model metric",
    color = NULL,
    fill = NULL
  ) +
  theme_project()
save_ggplot(
  validation_plot,
  file.path(figures_dir, "04_oplsda_permutation_validation"),
  8.4,
  6.2
)

make_feature_label <- function(name, id) {
  output <- ifelse(is.na(name) | trimws(name) == "", id, name)
  make.unique(output)
}
wrap_plot_label <- function(x, width = 34) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}
short_plot_label <- function(x, width = 34) {
  vapply(
    x,
    function(value) {
      if (nchar(value) <= width) value else
        paste0(substr(value, 1, width - 3), "...")
    },
    character(1)
  )
}
result_table$FeatureLabel <- make_feature_label(
  result_table$Metabolites, result_table$ID
)
top_vip <- head(result_table[order(-result_table$VIP), , drop = FALSE], top_vip_count)
write.csv(
  top_vip, file.path(tables_dir, "top_vip_metabolites.csv"),
  row.names = FALSE, quote = TRUE
)
top_vip$PlotLabel <- wrap_plot_label(top_vip$FeatureLabel)
top_vip$PlotLabel <- factor(
  top_vip$PlotLabel,
  levels = rev(top_vip$PlotLabel)
)
top_vip$PlotDirection <- ifelse(top_vip$log2FC > 0, "Up", "Down")
vip_plot <- ggplot(top_vip, aes(x = VIP, y = PlotLabel, fill = PlotDirection)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_vline(
    xintercept = vip_threshold, color = "#6B7280",
    linetype = "dashed", linewidth = 0.65
  ) +
  scale_fill_manual(values = direction_colors[c("Down", "Up")]) +
  labs(
    title = "Top OPLS-DA VIP metabolites",
    x = "Variable importance in projection (VIP)",
    y = NULL,
    fill = "Asphyxia vs Control"
  ) +
  theme_project() +
  theme(panel.grid.major.y = element_blank())
save_ggplot(
  vip_plot, file.path(figures_dir, "05_top_vip_metabolites"),
  10.0, max(9.2, 3.2 + 0.32 * nrow(top_vip))
)

s_plot_data <- result_table
s_plot_data$Label <- ""
s_positive <- which(s_plot_data$log2FC > 0)
s_negative <- which(s_plot_data$log2FC < 0)
s_label_index <- c(
  head(s_positive[order(-s_plot_data$VIP[s_positive])], 5),
  head(s_negative[order(-s_plot_data$VIP[s_negative])], 5)
)
s_plot_data$Label[s_label_index] <- short_plot_label(
  s_plot_data$FeatureLabel[s_label_index], 32
)
s_plot <- ggplot(
  s_plot_data,
  aes(x = PredictiveLoading, y = SCorrelation, color = log2FC)
) +
  geom_hline(yintercept = 0, color = "#A7AFBA", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#A7AFBA", linewidth = 0.4) +
  geom_point(alpha = 0.72, size = 2.0) +
  geom_text_repel(
    data = s_plot_data[s_plot_data$Label != "", , drop = FALSE],
    aes(label = Label), size = 2.8, max.overlaps = Inf,
    box.padding = 0.55, point.padding = 0.25, force = 2,
    seed = 20260921,
    show.legend = FALSE
  ) +
  scale_color_gradient2(
    low = "#4C78A8", mid = "#D6DAE0", high = "#E45756", midpoint = 0
  ) +
  labs(
    title = "OPLS-DA S-plot",
    x = "Predictive covariance loading",
    y = "Predictive correlation loading",
    color = "log2 fold change"
  ) +
  theme_project()
save_ggplot(s_plot, file.path(figures_dir, "06_oplsda_s_plot"), 8.4, 6.6)

volcano_data <- result_table
volcano_data$MinusLog10P <- -log10(pmax(volcano_data$Pvalue, 1e-300))
volcano_data$Direction <- factor(
  volcano_data$Direction, levels = c("Down", "Not selected", "Up")
)
volcano_candidates <- volcano_data[
  volcano_data$VIPCandidate, , drop = FALSE
]
volcano_up <- volcano_candidates[volcano_candidates$log2FC > 0, , drop = FALSE]
volcano_down <- volcano_candidates[volcano_candidates$log2FC < 0, , drop = FALSE]
volcano_labels <- rbind(
  head(volcano_up[order(-volcano_up$VIP), , drop = FALSE], 5),
  head(volcano_down[order(-volcano_down$VIP), , drop = FALSE], 5)
)
volcano_labels$PlotLabel <- short_plot_label(volcano_labels$FeatureLabel, 32)
volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FC, y = MinusLog10P, color = Direction)
) +
  geom_hline(
    yintercept = -log10(pvalue_threshold),
    color = "#6B7280", linetype = "dashed", linewidth = 0.6
  ) +
  geom_vline(
    xintercept = c(-logfc_threshold, logfc_threshold),
    color = "#6B7280", linetype = "dashed", linewidth = 0.6
  ) +
  geom_point(alpha = 0.74, size = 2.0) +
  geom_text_repel(
    data = volcano_labels,
    aes(label = PlotLabel), size = 2.9, max.overlaps = Inf,
    box.padding = 0.55, point.padding = 0.25, force = 2,
    seed = 20260921, show.legend = FALSE
  ) +
  scale_color_manual(values = direction_colors, drop = FALSE) +
  labs(
    title = "Differential metabolite abundance",
    subtitle = "Asphyxia vs Control",
    x = "log2 fold change",
    y = expression(-log[10](P)),
    color = "VIP candidate"
  ) +
  theme_project()
save_ggplot(
  volcano_plot, file.path(figures_dir, "07_differential_metabolite_volcano"),
  8.6, 6.6
)

candidate_table <- result_table[result_table$VIPCandidate, , drop = FALSE]
heatmap_features <- head(
  as.character(candidate_table$FeatureKey[order(-candidate_table$VIP)]),
  top_vip_count
)
if (length(heatmap_features) > 1) {
  heatmap_matrix <- expression_filtered[
    heatmap_features, biological_samples, drop = FALSE
  ]
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  heatmap_matrix[heatmap_matrix > 2] <- 2
  heatmap_matrix[heatmap_matrix < -2] <- -2
  label_map <- setNames(
    short_plot_label(result_table$FeatureLabel, 46), result_table$FeatureKey
  )
  rownames(heatmap_matrix) <- make.unique(label_map[rownames(heatmap_matrix)])
  column_annotation <- data.frame(Group = group)
  rownames(column_annotation) <- biological_samples
  annotation_colors <- list(
    Group = group_colors[c("Control", "Asphyxia")]
  )
  heatmap_colors <- colorRampPalette(c("#3B6FB6", "#F4F5F7", "#D95555"))(101)

  draw_heatmap <- function(filename, device_type) {
    if (device_type == "pdf") {
      cairo_pdf(filename, width = 9.4, height = 9.2, bg = "white")
    } else {
      png(
        filename, width = 9.4, height = 9.2, units = "in",
        res = 320, type = "cairo", bg = "white"
      )
    }
    pheatmap(
      heatmap_matrix,
      color = heatmap_colors,
      breaks = seq(-2, 2, length.out = 102),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = column_annotation,
      annotation_colors = annotation_colors,
      show_colnames = TRUE,
      border_color = NA,
      fontsize = 10,
      fontsize_row = 8.5,
      fontsize_col = 9,
      angle_col = 45,
      main = "Top OPLS-DA VIP candidate metabolites"
    )
    dev.off()
  }
  draw_heatmap(
    file.path(figures_dir, "08_vip_candidate_heatmap.pdf"), "pdf"
  )
  draw_heatmap(
    file.path(figures_dir, "08_vip_candidate_heatmap.png"), "png"
  )
}

top_box_features <- head(candidate_table$FeatureKey[order(-candidate_table$VIP)], 6)
if (length(top_box_features) > 0) {
  boxplot_label_map <- setNames(
    wrap_plot_label(
      short_plot_label(
        unname(setNames(
          result_table$FeatureLabel, result_table$FeatureKey
        )[top_box_features]),
        width = 34
      ),
      width = 17
    ),
    top_box_features
  )
  boxplot_data <- do.call(
    rbind,
    lapply(top_box_features, function(feature) {
      data.frame(
        FeatureKey = feature,
        Metabolite = unname(boxplot_label_map[feature]),
        Sample = biological_samples,
        Group = group,
        Abundance = as.numeric(expression_filtered[feature, biological_samples]),
        stringsAsFactors = FALSE
      )
    })
  )
  boxplot_data$Group <- factor(
    boxplot_data$Group, levels = c("Control", "Asphyxia")
  )
  boxplot_data$Metabolite <- factor(
    boxplot_data$Metabolite,
    levels = unname(boxplot_label_map)
  )
  boxplot_figure <- ggplot(
    boxplot_data,
    aes(x = Group, y = Abundance, color = Group, fill = Group)
  ) +
    geom_boxplot(width = 0.58, alpha = 0.16, outlier.shape = NA, linewidth = 0.75) +
    geom_jitter(width = 0.10, size = 2.8, alpha = 0.90) +
    facet_wrap(
      ~Metabolite, scales = "free_y", ncol = 3
    ) +
    scale_color_manual(values = group_colors[c("Control", "Asphyxia")]) +
    scale_fill_manual(values = group_colors[c("Control", "Asphyxia")]) +
    labs(
      title = "Top VIP metabolite abundance",
      x = NULL,
      y = "Processed log2 abundance"
    ) +
    theme_project() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", color = "#273142", size = 9)
    )
  save_ggplot(
    boxplot_figure,
    file.path(figures_dir, "09_top_vip_metabolite_abundance"),
    10.2,
    8.4
  )
}

summary_table <- data.frame(
  Metric = c(
    "Biological samples",
    "Control samples",
    "Asphyxia samples",
    "QC samples",
    "Input metabolite features",
    "Retained metabolite features",
    "Median QC RSD (%)",
    "OPLS-DA predictive components",
    "OPLS-DA orthogonal components",
    "R2X cumulative",
    "R2Y cumulative",
    "Q2 cumulative",
    "Permutation p-value R2Y",
    "Permutation p-value Q2",
    "Metabolites with VIP > 1",
    "Metabolites with nominal P < 0.05",
    "FDR- and effect-size-filtered metabolites",
    "VIP candidate metabolites"
  ),
  Value = c(
    nrow(biological_metadata),
    sum(group == "Control"),
    sum(group == "Asphyxia"),
    nrow(qc_metadata),
    nrow(processed),
    sum(keep_feature),
    median(qc_rsd[is.finite(qc_rsd)]),
    model_summary$pre,
    model_summary$ort,
    model_summary$`R2X(cum)`,
    model_summary$`R2Y(cum)`,
    model_summary$`Q2(cum)`,
    model_summary$pR2Y,
    model_summary$pQ2,
    sum(result_table$VIP > vip_threshold),
    sum(result_table$Pvalue < pvalue_threshold),
    sum(result_table$FDRSignificant),
    sum(result_table$VIPCandidate)
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary_table, file.path(output_dir, "analysis_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
parameters <- data.frame(
  Parameter = c(
    "Comparison", "Processed worksheet", "Original worksheet",
    "Biological detection requirement", "QC RSD threshold", "Scaling",
    "OPLS-DA algorithm", "Predictive components", "Orthogonal components",
    "Cross-validation folds", "Permutations", "VIP threshold",
    "Nominal P threshold", "FDR threshold", "Absolute log2FC threshold",
    "VIP candidate definition"
  ),
  Value = c(
    "Asphyxia vs Control", processed_sheet, original_sheet,
    paste0(">= ", detection_threshold * 100, "% in either group"),
    paste0("<= ", qc_rsd_threshold, "%"), "Pareto",
    "NIPALS implemented by ropls", 1, 1,
    cross_validation_folds, permutation_count, vip_threshold,
    pvalue_threshold, fdr_threshold, logfc_threshold,
    "VIP > 1; nominal P < 0.05; absolute log2FC >= 0.5"
  ),
  stringsAsFactors = FALSE
)
write.table(
  parameters, file.path(output_dir, "analysis_parameters.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "session_info.txt"))

message(
  "Metabolomics OPLS-DA complete. Retained features: ", sum(keep_feature),
  "; VIP candidates: ", sum(result_table$VIPCandidate),
  "; FDR significant: ", sum(result_table$FDRSignificant), "."
)
message("Results written to: ", normalizePath(output_dir))
