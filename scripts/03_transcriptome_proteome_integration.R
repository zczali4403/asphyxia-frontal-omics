#!/usr/bin/env Rscript

# Transcriptome-proteome integration for the R1 frontal cortex.
# Comparison: Asphyxia vs Control.

required_packages <- c("ggplot2", "ggrepel", "ggh4x")
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
  library(ggplot2)
  library(ggrepel)
  library(ggh4x)
})

args <- commandArgs(trailingOnly = TRUE)
transcript_results_file <- if (length(args) >= 1) args[[1]] else
  "results/01_transcriptomics/tables/deseq2_all_genes.csv"
protein_results_file <- if (length(args) >= 2) args[[2]] else
  "results/02_proteomics/tables/limma_all_proteins.csv"
transcript_expression_file <- if (length(args) >= 3) args[[3]] else
  "results/01_transcriptomics/tables/vst_expression_filtered.csv"
protein_expression_file <- if (length(args) >= 4) args[[4]] else
  "results/02_proteomics/tables/protein_expression_filtered.csv"
metadata_file <- if (length(args) >= 5) args[[5]] else
  "config/samples.csv"
output_dir <- if (length(args) >= 6) args[[6]] else
  "results/03_transcriptome_proteome_integration"

rna_alpha <- 0.05
protein_exploratory_alpha <- 0.05
lfc_threshold <- 0.5
top_labels <- 12L

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading transcriptomic differential results: ", transcript_results_file)
rna <- read.csv(
  transcript_results_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
message("Reading proteomic differential results: ", protein_results_file)
protein <- read.csv(
  protein_results_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_rna <- c(
  "GeneID", "GeneName", "baseMean", "log2FoldChange", "pvalue", "padj"
)
required_protein <- c(
  "Accession", "GeneName", "AveExpr", "log2FoldChange", "pvalue", "padj"
)
if (!all(required_rna %in% colnames(rna))) {
  stop("The transcriptomic result table is missing required columns.", call. = FALSE)
}
if (!all(required_protein %in% colnames(protein))) {
  stop("The proteomic result table is missing required columns.", call. = FALSE)
}

valid_rna <- !is.na(rna$GeneName) & trimws(rna$GeneName) != ""
valid_protein <- !is.na(protein$GeneName) & trimws(protein$GeneName) != ""
rna <- rna[valid_rna, , drop = FALSE]
protein <- protein[valid_protein, , drop = FALSE]

# Choose the most abundant feature for the few duplicated gene symbols.
rna <- rna[order(rna$GeneName, -rna$baseMean), , drop = FALSE]
rna$RNAFeatureCount <- ave(rna$GeneName, rna$GeneName, FUN = length)
rna_representative <- rna[!duplicated(rna$GeneName), , drop = FALSE]
protein <- protein[order(protein$GeneName, -protein$AveExpr), , drop = FALSE]
protein$ProteinFeatureCount <- ave(
  protein$GeneName, protein$GeneName, FUN = length
)
protein_representative <- protein[
  !duplicated(protein$GeneName), , drop = FALSE
]

matched <- merge(
  rna_representative,
  protein_representative,
  by = "GeneName",
  suffixes = c("_RNA", "_Protein"),
  all = FALSE,
  sort = FALSE
)
if (nrow(matched) < 2) {
  stop("Too few gene-protein pairs were matched by GeneName.", call. = FALSE)
}

integration <- data.frame(
  GeneName = matched$GeneName,
  GeneID = matched$GeneID_RNA,
  Accession = matched$Accession,
  RNAFeatureCount = matched$RNAFeatureCount,
  ProteinFeatureCount = matched$ProteinFeatureCount,
  RNABaseMean = matched$baseMean,
  RNA_log2FC = matched$log2FoldChange_RNA,
  RNA_pvalue = matched$pvalue_RNA,
  RNA_padj = matched$padj_RNA,
  ProteinAveExpr = matched$AveExpr,
  Protein_log2FC = matched$log2FoldChange_Protein,
  Protein_pvalue = matched$pvalue_Protein,
  Protein_padj = matched$padj_Protein,
  stringsAsFactors = FALSE
)

integration$RNAEffect <- "Stable"
integration$RNAEffect[integration$RNA_log2FC >= lfc_threshold] <- "Up"
integration$RNAEffect[integration$RNA_log2FC <= -lfc_threshold] <- "Down"
integration$ProteinEffect <- "Stable"
integration$ProteinEffect[
  integration$Protein_log2FC >= lfc_threshold
] <- "Up"
integration$ProteinEffect[
  integration$Protein_log2FC <= -lfc_threshold
] <- "Down"

integration$RNAEvidence <-
  !is.na(integration$RNA_padj) &
  integration$RNA_padj < rna_alpha &
  abs(integration$RNA_log2FC) >= lfc_threshold
integration$ProteinEvidence <-
  !is.na(integration$Protein_pvalue) &
  integration$Protein_pvalue < protein_exploratory_alpha &
  abs(integration$Protein_log2FC) >= lfc_threshold

integration$Evidence <- "Neither"
integration$Evidence[integration$RNAEvidence] <- "RNA only"
integration$Evidence[integration$ProteinEvidence] <- "Protein only"
integration$Evidence[
  integration$RNAEvidence & integration$ProteinEvidence
] <- "Both"
integration$Evidence <- factor(
  integration$Evidence,
  levels = c("Neither", "RNA only", "Protein only", "Both")
)

integration$NineQuadrant <- "Stable in both"
integration$NineQuadrant[
  integration$RNAEffect == "Up" & integration$ProteinEffect == "Up"
] <- "Concordant up"
integration$NineQuadrant[
  integration$RNAEffect == "Down" & integration$ProteinEffect == "Down"
] <- "Concordant down"
integration$NineQuadrant[
  integration$RNAEffect != "Stable" &
    integration$ProteinEffect != "Stable" &
    integration$RNAEffect != integration$ProteinEffect
] <- "Discordant"
integration$NineQuadrant[
  integration$RNAEffect != "Stable" & integration$ProteinEffect == "Stable"
] <- "RNA-specific"
integration$NineQuadrant[
  integration$RNAEffect == "Stable" & integration$ProteinEffect != "Stable"
] <- "Protein-specific"
integration$NineQuadrant <- factor(
  integration$NineQuadrant,
  levels = c(
    "Stable in both", "RNA-specific", "Protein-specific",
    "Discordant", "Concordant down", "Concordant up"
  )
)

integration$NineCell <- paste(
  paste0("RNA ", tolower(integration$RNAEffect)),
  paste0("Protein ", tolower(integration$ProteinEffect)),
  sep = " / "
)

# Combined evidence score is used only to choose labels, not for inference.
safe_rna_p <- ifelse(is.na(integration$RNA_padj), 1, integration$RNA_padj)
safe_protein_p <- ifelse(
  is.na(integration$Protein_pvalue), 1, integration$Protein_pvalue
)
integration$LabelScore <- -log10(pmax(safe_rna_p, .Machine$double.xmin)) +
  -log10(pmax(safe_protein_p, .Machine$double.xmin))

integration <- integration[
  order(
    integration$Evidence,
    integration$LabelScore,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]
write.csv(
  integration,
  file.path(tables_dir, "matched_transcript_protein_results.csv"),
  row.names = FALSE,
  quote = TRUE
)

joint_candidates <- integration[
  integration$RNAEvidence | integration$ProteinEvidence,
  ,
  drop = FALSE
]
write.csv(
  joint_candidates,
  file.path(tables_dir, "cross_omics_candidate_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)
both_candidates <- integration[
  integration$RNAEvidence & integration$ProteinEvidence,
  ,
  drop = FALSE
]
write.csv(
  both_candidates,
  file.path(tables_dir, "rna_protein_evidence_overlap.csv"),
  row.names = FALSE,
  quote = TRUE
)

nine_levels <- c(
  "RNA down / Protein up", "RNA stable / Protein up", "RNA up / Protein up",
  "RNA down / Protein stable", "RNA stable / Protein stable",
  "RNA up / Protein stable", "RNA down / Protein down",
  "RNA stable / Protein down", "RNA up / Protein down"
)
nine_counts <- as.data.frame(
  table(factor(integration$NineCell, levels = nine_levels)),
  stringsAsFactors = FALSE
)
colnames(nine_counts) <- c("NineCell", "Genes")
write.csv(
  nine_counts,
  file.path(tables_dir, "nine_quadrant_counts.csv"),
  row.names = FALSE,
  quote = TRUE
)

group_colors <- c(Control = "#4C78A8", Asphyxia = "#E45756")
nine_cell_colors <- c(
  `RNA down / Protein up` = "#009E73",
  `RNA stable / Protein up` = "#E69F00",
  `RNA up / Protein up` = "#D62728",
  `RNA down / Protein stable` = "#0072B2",
  `RNA stable / Protein stable` = "#B8BEC8",
  `RNA up / Protein stable` = "#CC79A7",
  `RNA down / Protein down` = "#332288",
  `RNA stable / Protein down` = "#56B4E9",
  `RNA up / Protein down` = "#882255"
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

# Global fold-change correlation is retained as a summary statistic.
fc_cor <- suppressWarnings(cor.test(
  integration$RNA_log2FC,
  integration$Protein_log2FC,
  method = "spearman",
  exact = FALSE
))

# Nine-quadrant plot based on effect-size thresholds; evidence is stored separately.
nine_plot_data <- integration
nine_plot_data$NineCell <- factor(
  nine_plot_data$NineCell,
  levels = nine_levels
)
nine_background <- data.frame(
  xmin = rep(c(-Inf, -lfc_threshold, lfc_threshold), times = 3),
  xmax = rep(c(-lfc_threshold, lfc_threshold, Inf), times = 3),
  ymin = rep(c(lfc_threshold, -lfc_threshold, -Inf), each = 3),
  ymax = rep(c(Inf, lfc_threshold, -lfc_threshold), each = 3),
  NineCell = factor(nine_levels, levels = nine_levels)
)
legend_key_data <- data.frame(
  RNA_log2FC = rep(Inf, length(nine_levels)),
  Protein_log2FC = rep(Inf, length(nine_levels)),
  NineCell = factor(nine_levels, levels = nine_levels)
)
nine_plot_data <- nine_plot_data[
  order(nine_plot_data$NineQuadrant == "Stable in both", decreasing = TRUE),
  ,
  drop = FALSE
]
nine_labels <- nine_plot_data[
  nine_plot_data$Evidence != "Neither", , drop = FALSE
]
nine_labels <- head(
  nine_labels[order(-nine_labels$LabelScore), , drop = FALSE],
  top_labels
)
if ("Plcxd2" %in% nine_plot_data$GeneName) {
  nine_labels <- unique(rbind(
    nine_labels,
    nine_plot_data[nine_plot_data$GeneName == "Plcxd2", , drop = FALSE]
  ))
}

nine_plot <- ggplot(
  nine_plot_data,
  aes(
    x = RNA_log2FC,
    y = Protein_log2FC,
    color = NineCell
  )
) +
  geom_rect(
    data = nine_background,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = NineCell),
    inherit.aes = FALSE, alpha = 0.07, color = NA, show.legend = FALSE
  ) +
  geom_vline(
    xintercept = c(-lfc_threshold, lfc_threshold),
    linetype = "dashed", color = "#68707C", linewidth = 0.45
  ) +
  geom_hline(
    yintercept = c(-lfc_threshold, lfc_threshold),
    linetype = "dashed", color = "#68707C", linewidth = 0.45
  ) +
  geom_vline(xintercept = 0, color = "#B7BDC6", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "#B7BDC6", linewidth = 0.35) +
  geom_point(alpha = 0.72, size = 1.65) +
  geom_point(
    data = legend_key_data,
    aes(x = RNA_log2FC, y = Protein_log2FC, color = NineCell),
    inherit.aes = FALSE, alpha = 0, size = 0, show.legend = TRUE
  ) +
  geom_text_repel(
    data = nine_labels,
    aes(label = GeneName),
    size = 3.0, box.padding = 0.42, point.padding = 0.25,
    max.overlaps = Inf, min.segment.length = 0,
    seed = 20260920, show.legend = FALSE
  ) +
  scale_color_manual(
    values = nine_cell_colors,
    breaks = nine_levels,
    drop = FALSE
  ) +
  scale_fill_manual(values = nine_cell_colors, guide = "none") +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 3))
  ) +
  labs(
    title = "Transcript-protein nine-quadrant analysis",
    subtitle = "asphyxia vs control",
    x = "RNA log2 fold change",
    y = "Protein log2 fold change",
    color = "Expression pattern"
  ) +
  theme_project()
save_ggplot(nine_plot, file.path(figures_dir, "01_nine_quadrant"), 8.6, 7.0)

# Compact count plot for the nine cells, excluding neither axis information.
nine_counts$NineCell <- factor(nine_counts$NineCell, levels = nine_levels)
nine_counts$GenesPlot <- log10(nine_counts$Genes + 1)
count_plot <- ggplot(nine_counts, aes(x = NineCell, y = GenesPlot)) +
  geom_col(width = 0.72, fill = "#607D9B") +
  geom_text(
    aes(label = Genes), hjust = -0.15, size = 3.5,
    fontface = "bold", color = "#374151"
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = log10(c(0, 10, 100, 1000, 7000) + 1),
    labels = c("0", "10", "100", "1,000", "7,000"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Nine-quadrant gene counts",
    subtitle = "|log2 fold change| = 0.5",
    x = NULL,
    y = "Matched genes (log10 scale)"
  ) +
  theme_project() +
  theme(legend.position = "none")
save_ggplot(count_plot, file.path(figures_dir, "02_nine_quadrant_counts"), 8.4, 6.2)

# Paired individual-level RNA-protein correlations for subjects 1-3 per group.
message("Reading paired expression matrices")
rna_expression <- read.csv(
  transcript_expression_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
protein_expression <- read.csv(
  protein_expression_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

rna_metadata <- metadata[
  metadata$omics == "transcriptomics" &
    metadata$sample_role == "biological" &
    metadata$region == "frontal_cortex",
  c("sample_id", "subject_id", "group"),
  drop = FALSE
]
protein_metadata <- metadata[
  metadata$omics == "proteomics" &
    metadata$sample_role == "biological" &
    metadata$region == "frontal_cortex",
  c("sample_id", "subject_id", "group"),
  drop = FALSE
]
paired_samples <- merge(
  rna_metadata,
  protein_metadata,
  by = c("subject_id", "group"),
  suffixes = c("_RNA", "_Protein")
)
paired_samples$group[paired_samples$group == "CON"] <- "Control"
paired_samples$group <- factor(
  paired_samples$group,
  levels = c("Control", "Asphyxia")
)
paired_samples <- paired_samples[
  order(paired_samples$group, paired_samples$subject_id),
  ,
  drop = FALSE
]
if (nrow(paired_samples) < 4) {
  stop("Too few cross-omics paired animals were found.", call. = FALSE)
}
write.csv(
  paired_samples,
  file.path(tables_dir, "paired_samples_used.csv"),
  row.names = FALSE,
  quote = TRUE
)

rna_index <- match(integration$GeneID, rna_expression$GeneID)
protein_index <- match(integration$Accession, protein_expression$Accession)
valid_expression_pairs <- !is.na(rna_index) & !is.na(protein_index)
paired_integration <- integration[valid_expression_pairs, , drop = FALSE]
rna_index <- rna_index[valid_expression_pairs]
protein_index <- protein_index[valid_expression_pairs]

rna_paired_matrix <- as.matrix(
  rna_expression[
    rna_index,
    paired_samples$sample_id_RNA,
    drop = FALSE
  ]
)
protein_paired_matrix <- as.matrix(
  protein_expression[
    protein_index,
    paired_samples$sample_id_Protein,
    drop = FALSE
  ]
)
storage.mode(rna_paired_matrix) <- "numeric"
storage.mode(protein_paired_matrix) <- "numeric"
rownames(rna_paired_matrix) <- paired_integration$GeneName
rownames(protein_paired_matrix) <- paired_integration$GeneName

correlation_rows <- lapply(seq_len(nrow(paired_integration)), function(i) {
  x <- rna_paired_matrix[i, ]
  y <- protein_paired_matrix[i, ]
  pearson <- suppressWarnings(cor.test(x, y, method = "pearson"))
  spearman <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    GeneName = paired_integration$GeneName[[i]],
    GeneID = paired_integration$GeneID[[i]],
    Accession = paired_integration$Accession[[i]],
    PearsonR = unname(pearson$estimate),
    PearsonP = pearson$p.value,
    SpearmanRho = unname(spearman$estimate),
    SpearmanP = spearman$p.value,
    stringsAsFactors = FALSE
  )
})
paired_correlations <- do.call(rbind, correlation_rows)
paired_correlations$PearsonPadj <- p.adjust(
  paired_correlations$PearsonP,
  method = "BH"
)
paired_correlations$SpearmanPadj <- p.adjust(
  paired_correlations$SpearmanP,
  method = "BH"
)
paired_correlations <- paired_correlations[
  order(paired_correlations$SpearmanP, -abs(paired_correlations$SpearmanRho)),
  ,
  drop = FALSE
]
write.csv(
  paired_correlations,
  file.path(tables_dir, "paired_rna_protein_correlations.csv"),
  row.names = FALSE,
  quote = TRUE
)

# Group-level expression plot for the strongest shared evidence candidate.
candidate_gene <- if ("Plcxd2" %in% paired_integration$GeneName) {
  "Plcxd2"
} else if (nrow(both_candidates) > 0) {
  both_candidates$GeneName[[1]]
} else {
  paired_correlations$GeneName[[1]]
}
candidate_index <- match(candidate_gene, rownames(rna_paired_matrix))
candidate_data <- data.frame(
  Subject = paired_samples$subject_id,
  Group = paired_samples$group,
  RNA = as.numeric(scale(rna_paired_matrix[candidate_index, ])),
  Protein = as.numeric(scale(protein_paired_matrix[candidate_index, ])),
  stringsAsFactors = FALSE
)
candidate_data$DisplaySubject <- sub(
  "^CON(?=-)", "Control", candidate_data$Subject, perl = TRUE
)
candidate_cor <- suppressWarnings(cor.test(
  candidate_data$RNA,
  candidate_data$Protein,
  method = "spearman",
  exact = FALSE
))

candidate_result <- paired_integration[candidate_index, , drop = FALSE]
rna_candidate_index <- match(candidate_result$GeneID, rna_expression$GeneID)
protein_candidate_index <- match(
  candidate_result$Accession,
  protein_expression$Accession
)
rna_candidate_values <- as.numeric(
  rna_expression[rna_candidate_index, rna_metadata$sample_id]
)
protein_candidate_values <- as.numeric(
  protein_expression[protein_candidate_index, protein_metadata$sample_id]
)
candidate_group_data <- rbind(
  data.frame(
    SampleID = rna_metadata$sample_id,
    Group = rna_metadata$group,
    Omics = "Transcriptomics",
    Expression = rna_candidate_values,
    Zscore = as.numeric(scale(rna_candidate_values)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    SampleID = protein_metadata$sample_id,
    Group = protein_metadata$group,
    Omics = "Proteomics",
    Expression = protein_candidate_values,
    Zscore = as.numeric(scale(protein_candidate_values)),
    stringsAsFactors = FALSE
  )
)
candidate_group_data$Group <- factor(
  candidate_group_data$Group,
  levels = c("Control", "Asphyxia")
)
candidate_group_data$Omics <- factor(
  candidate_group_data$Omics,
  levels = c("Transcriptomics", "Proteomics")
)
write.csv(
  candidate_group_data,
  file.path(tables_dir, "candidate_expression_across_omics.csv"),
  row.names = FALSE,
  quote = TRUE
)

candidate_plot <- ggplot(
  candidate_group_data,
  aes(x = Group, y = Zscore, color = Group, fill = Group)
) +
  geom_boxplot(
    width = 0.52, alpha = 0.14, linewidth = 0.65,
    outlier.shape = NA, show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0, seed = 20260920),
    size = 3.2, alpha = 0.95, show.legend = FALSE
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 3.6,
    color = "#273142", fill = "white", stroke = 0.8,
    show.legend = FALSE
  ) +
  ggh4x::facet_wrap2(
    ~Omics,
    nrow = 1,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = c("#E9E1F4", "#DDF1E5"),
        color = NA
      ),
      text_x = ggh4x::elem_list_text(
        color = c("#6B4FA3", "#2E7D5B"),
        face = "bold",
        size = 12
      )
    )
  ) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  labs(
    title = paste0(candidate_gene, " expression across omics"),
    subtitle = "asphyxia vs control",
    x = NULL,
    y = "Standardized expression (z-score)"
  ) +
  theme_project() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )
save_ggplot(
  candidate_plot,
  file.path(figures_dir, "03_candidate_expression_across_omics"),
  8.2,
  5.8
)

summary_table <- data.frame(
  Metric = c(
    "Matched transcript-protein genes",
    "RNA evidence genes among matched genes",
    "Exploratory protein evidence genes among matched genes",
    "Genes with evidence in both omics",
    "Concordant-up effect genes",
    "Concordant-down effect genes",
    "Discordant effect genes",
    "RNA-specific effect genes",
    "Protein-specific effect genes",
    "Stable-in-both genes",
    "Fold-change Spearman rho",
    "Fold-change Spearman P value",
    "Paired animals",
    paste0(candidate_gene, " paired-expression Spearman rho"),
    paste0(candidate_gene, " paired-expression Spearman P value"),
    "Absolute log2FC threshold",
    "RNA adjusted P-value threshold",
    "Protein nominal P-value threshold"
  ),
  Value = c(
    nrow(integration),
    sum(integration$RNAEvidence),
    sum(integration$ProteinEvidence),
    sum(integration$RNAEvidence & integration$ProteinEvidence),
    sum(integration$NineQuadrant == "Concordant up"),
    sum(integration$NineQuadrant == "Concordant down"),
    sum(integration$NineQuadrant == "Discordant"),
    sum(integration$NineQuadrant == "RNA-specific"),
    sum(integration$NineQuadrant == "Protein-specific"),
    sum(integration$NineQuadrant == "Stable in both"),
    unname(fc_cor$estimate),
    fc_cor$p.value,
    nrow(paired_samples),
    unname(candidate_cor$estimate),
    candidate_cor$p.value,
    lfc_threshold,
    rna_alpha,
    protein_exploratory_alpha
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
  "Integration complete. Matched genes: ", nrow(integration),
  "; RNA-protein evidence overlap: ",
  sum(integration$RNAEvidence & integration$ProteinEvidence),
  "; paired animals: ", nrow(paired_samples), "."
)
message("Results written to: ", normalizePath(output_dir))
