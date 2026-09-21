#!/usr/bin/env Rscript

# Local KEGG pathway analysis and Plcxd2-centered metabolomics integration.
# Correlations are unadjusted Pearson correlations across matched animals and
# therefore describe co-response, not direct regulation by Plcxd2.

required_packages <- c(
  "ggplot2", "ggrepel", "ggraph", "igraph", "GSVA", "scales",
  "BiocParallel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})
serial_param <- BiocParallel::SerialParam(progressbar = FALSE)
BiocParallel::register(serial_param, default = TRUE)

args <- commandArgs(trailingOnly = TRUE)
metabolite_results_file <- if (length(args) >= 1) args[[1]] else
  "results/05_metabolomics/tables/oplsda_vip_all_metabolites.csv"
metabolite_expression_file <- if (length(args) >= 2) args[[2]] else
  "results/05_metabolomics/tables/metabolite_expression_filtered.csv"
rna_file <- if (length(args) >= 3) args[[3]] else
  "results/01_transcriptomics/tables/vst_expression_filtered.csv"
protein_file <- if (length(args) >= 4) args[[4]] else
  "results/02_proteomics/tables/protein_expression_filtered.csv"
metadata_file <- if (length(args) >= 5) args[[5]] else "config/samples.csv"
output_dir <- if (length(args) >= 6) args[[6]] else
  "results/06_plcxd2_metabolomics_integration"

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
legacy_network_outputs <- c(
  file.path(
    figures_dir,
    paste0("04_plcxd2_associated_metabolite_network", c(".png", ".pdf"))
  ),
  file.path(
    tables_dir,
    c(
      "plcxd2_metabolite_network_nodes.csv",
      "plcxd2_metabolite_network_edges.csv"
    )
  )
)
unlink(legacy_network_outputs, force = TRUE)
legacy_pathway_outputs <- c(
  file.path(
    figures_dir,
    c(
      "01_ranked_metabolic_pathway_enrichment.png",
      "01_ranked_metabolic_pathway_enrichment.pdf",
      "02_vip_candidate_pathway_enrichment.png",
      "02_vip_candidate_pathway_enrichment.pdf"
    )
  ),
  file.path(tables_dir, "ranked_metabolite_pathway_enrichment.csv")
)
unlink(legacy_pathway_outputs, force = TRUE)

target_symbol <- "Plcxd2"
min_pathway_size <- 3L
max_pathway_size <- 500L
set.seed(20260921)

group_colors <- c(Control = "#4C78A8", Asphyxia = "#E45756")
direction_colors <- c(Positive = "#D95F59", Negative = "#4C78A8")
omics_colors <- c(Transcriptomics = "#6B4FA3", Proteomics = "#2E8B68")

theme_project <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold", size = base_size + 3, hjust = 0.5,
        color = "#111827"
      ),
      plot.subtitle = element_text(
        color = "#4B5563", hjust = 0.5, margin = margin(b = 10)
      ),
      axis.title = element_text(face = "bold", color = "#273142"),
      axis.text = element_text(color = "#445064"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E3E8EF", linewidth = 0.45),
      plot.margin = margin(12, 16, 12, 16)
    )
}

save_ggplot <- function(plot_object, path_without_extension, width, height) {
  ggsave(
    paste0(path_without_extension, ".png"), plot_object,
    width = width, height = height, dpi = 320, bg = "white"
  )
  ggsave(
    paste0(path_without_extension, ".pdf"), plot_object,
    width = width, height = height, device = cairo_pdf, bg = "white"
  )
}

short_label <- function(x, width = 24) {
  x <- ifelse(is.na(x) | trimws(x) == "", "Unannotated metabolite", x)
  vapply(
    x,
    function(value) {
      if (nchar(value) <= width) value else paste0(substr(value, 1, width - 3), "...")
    },
    character(1)
  )
}

wrap_label <- function(x, width = 20) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

message("Reading metabolomics, expression, and metadata tables")
metabolite_results <- read.csv(
  metabolite_results_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
metabolite_expression <- read.csv(
  metabolite_expression_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
rna <- read.csv(
  rna_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
protein <- read.csv(
  protein_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
metadata <- read.csv(
  metadata_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_metabolite_columns <- c(
  "ID", "FeatureKey", "Metabolites", "KEGG", "ID Annotation", "Annotation",
  "Class", "Sub Class", "t", "Pvalue", "padj", "VIP", "VIPCandidate",
  "log2FC"
)
if (!all(required_metabolite_columns %in% colnames(metabolite_results))) {
  stop("The metabolomics result table is missing required columns.", call. = FALSE)
}
if (!all(c("GeneID", "GeneName") %in% colnames(rna))) {
  stop("The RNA table must contain GeneID and GeneName.", call. = FALSE)
}
if (!all(c("Accession", "GeneName") %in% colnames(protein))) {
  stop("The protein table must contain Accession and GeneName.", call. = FALSE)
}

metabolite_results$FeatureLabel <- ifelse(
  is.na(metabolite_results$Metabolites) |
    trimws(metabolite_results$Metabolites) == "",
  metabolite_results$FeatureKey,
  metabolite_results$Metabolites
)

parse_pathway_membership <- function(result_table) {
  output <- vector("list", nrow(result_table))
  for (i in seq_len(nrow(result_table))) {
    id_text <- result_table[["ID Annotation"]][i]
    if (is.na(id_text) || trimws(id_text) == "") next
    ids <- trimws(strsplit(id_text, ",", fixed = TRUE)[[1]])
    ids <- ids[grepl("^mmu[0-9]{5}$", ids)]
    if (length(ids) == 0) next
    name_text <- result_table$Annotation[i]
    names_i <- if (is.na(name_text) || trimws(name_text) == "") {
      rep(NA_character_, length(ids))
    } else {
      trimws(strsplit(name_text, "|", fixed = TRUE)[[1]])
    }
    if (length(names_i) < length(ids)) {
      names_i <- c(names_i, rep(NA_character_, length(ids) - length(names_i)))
    }
    names_i <- names_i[seq_along(ids)]
    names_i[is.na(names_i) | names_i == ""] <- ids[is.na(names_i) | names_i == ""]
    output[[i]] <- data.frame(
      FeatureKey = result_table$FeatureKey[i],
      FeatureLabel = result_table$FeatureLabel[i],
      KEGGCompound = result_table$KEGG[i],
      PathwayID = ids,
      PathwayName = names_i,
      stringsAsFactors = FALSE
    )
  }
  output <- output[!vapply(output, is.null, logical(1))]
  if (length(output) == 0) {
    return(data.frame())
  }
  unique(do.call(rbind, output))
}

message("Parsing local KEGG metabolite-pathway annotations")
membership <- parse_pathway_membership(metabolite_results)
if (nrow(membership) == 0) {
  stop("No local KEGG pathway annotations were found.", call. = FALSE)
}
write.csv(
  membership,
  file.path(tables_dir, "metabolite_kegg_pathway_membership.csv"),
  row.names = FALSE, quote = TRUE
)

background_features <- unique(membership$FeatureKey)
candidate_features <- intersect(
  metabolite_results$FeatureKey[metabolite_results$VIPCandidate %in% TRUE],
  background_features
)
pathway_ids <- unique(membership$PathwayID)
ora_rows <- lapply(pathway_ids, function(pathway_id) {
  pathway_features <- unique(membership$FeatureKey[membership$PathwayID == pathway_id])
  pathway_name <- membership$PathwayName[match(
    pathway_id, membership$PathwayID
  )]
  universe_size <- length(background_features)
  candidate_size <- length(candidate_features)
  pathway_size <- length(pathway_features)
  overlap <- intersect(pathway_features, candidate_features)
  overlap_size <- length(overlap)
  if (pathway_size < min_pathway_size || candidate_size == 0) return(NULL)
  pvalue <- phyper(
    overlap_size - 1, pathway_size, universe_size - pathway_size,
    candidate_size, lower.tail = FALSE
  )
  fold_enrichment <- if (overlap_size == 0) 0 else
    (overlap_size / candidate_size) / (pathway_size / universe_size)
  data.frame(
    PathwayID = pathway_id,
    PathwayName = pathway_name,
    BackgroundSize = universe_size,
    CandidateSize = candidate_size,
    PathwaySize = pathway_size,
    OverlapSize = overlap_size,
    FoldEnrichment = fold_enrichment,
    Pvalue = pvalue,
    OverlapFeatures = paste(overlap, collapse = "/"),
    stringsAsFactors = FALSE
  )
})
ora_rows <- ora_rows[!vapply(ora_rows, is.null, logical(1))]
ora <- do.call(rbind, ora_rows)
ora$padj <- p.adjust(ora$Pvalue, method = "BH")
ora <- ora[order(ora$Pvalue, -ora$FoldEnrichment), , drop = FALSE]
write.csv(
  ora, file.path(tables_dir, "vip_candidate_pathway_ora.csv"),
  row.names = FALSE, quote = TRUE
)

pathway_list <- split(membership$FeatureKey, membership$PathwayID)
pathway_list <- lapply(pathway_list, unique)
pathway_list <- pathway_list[
  vapply(pathway_list, length, integer(1)) >= min_pathway_size
]
pathway_name_map <- setNames(membership$PathwayName, membership$PathwayID)

ora_plot_data <- head(ora[ora$OverlapSize > 0, , drop = FALSE], 15)
if (nrow(ora_plot_data) > 0) {
  ora_plot_data$PlotLabel <- factor(
    paste0("[KEGG] ", ora_plot_data$PathwayName),
    levels = rev(paste0("[KEGG] ", ora_plot_data$PathwayName))
  )
  ora_plot_data$MinusLog10P <- -log10(pmax(ora_plot_data$Pvalue, 1e-300))
  ora_plot <- ggplot(
    ora_plot_data,
    aes(
      x = FoldEnrichment, y = PlotLabel, size = OverlapSize,
      color = MinusLog10P
    )
  ) +
    geom_point(alpha = 0.92) +
    scale_color_gradient(low = "#7DB9A6", high = "#2E6F5E") +
    scale_size_continuous(range = c(3.5, 10)) +
    labs(
      title = "VIP candidate pathway enrichment",
      subtitle = "Asphyxia vs Control",
      x = "Fold enrichment", y = NULL,
      color = expression(-log[10](P)), size = "Overlap"
    ) +
    theme_project(11) +
    theme(axis.text.y = element_text(size = 9.3))
  save_ggplot(
    ora_plot,
    file.path(figures_dir, "01_vip_candidate_pathway_enrichment"),
    10.8, 7.6
  )
}

get_omics_metadata <- function(omics_name) {
  metadata[
    metadata$omics == omics_name &
      metadata$sample_role == "biological" &
      metadata$region == "frontal_cortex",
    ,
    drop = FALSE
  ]
}

calculate_metabolite_correlations <- function(
  target_table, target_id_column, target_symbol_column,
  target_metadata, metabolite_metadata, omics_label
) {
  shared_subjects <- intersect(
    target_metadata$subject_id, metabolite_metadata$subject_id
  )
  target_metadata <- target_metadata[
    match(shared_subjects, target_metadata$subject_id), , drop = FALSE
  ]
  metabolite_metadata <- metabolite_metadata[
    match(shared_subjects, metabolite_metadata$subject_id), , drop = FALSE
  ]
  target_samples <- target_metadata$sample_id
  metabolite_samples <- metabolite_metadata$sample_id
  if (!all(target_samples %in% colnames(target_table)) ||
      !all(metabolite_samples %in% colnames(metabolite_expression))) {
    stop("Matched sample columns are missing for ", omics_label, call. = FALSE)
  }
  target_rows <- which(
    tolower(target_table[[target_symbol_column]]) == tolower(target_symbol)
  )
  if (length(target_rows) == 0) {
    stop(target_symbol, " was not found in ", omics_label, call. = FALSE)
  }
  target_matrix <- as.matrix(
    target_table[target_rows, target_samples, drop = FALSE]
  )
  storage.mode(target_matrix) <- "numeric"
  target_row <- target_rows[which.max(rowMeans(target_matrix))]
  target_values <- as.numeric(target_table[target_row, target_samples])
  metabolite_matrix <- as.matrix(
    metabolite_expression[, metabolite_samples, drop = FALSE]
  )
  storage.mode(metabolite_matrix) <- "numeric"
  target_centered <- target_values - mean(target_values)
  metabolite_centered <- metabolite_matrix - rowMeans(metabolite_matrix)
  target_ss <- sum(target_centered^2)
  feature_ss <- rowSums(metabolite_centered^2)
  correlation <- as.numeric(metabolite_centered %*% target_centered) /
    sqrt(feature_ss * target_ss)
  correlation[!is.finite(correlation)] <- NA_real_
  correlation <- pmax(pmin(correlation, 1), -1)
  n_samples <- length(shared_subjects)
  correlation_t <- correlation * sqrt(
    (n_samples - 2) / pmax(1 - correlation^2, .Machine$double.eps)
  )
  pvalue <- 2 * pt(-abs(correlation_t), df = n_samples - 2)
  annotation_match <- match(metabolite_expression$ID, metabolite_results$ID)
  output <- data.frame(
    FeatureKey = metabolite_results$FeatureKey[annotation_match],
    Metabolite = metabolite_results$FeatureLabel[annotation_match],
    KEGG = metabolite_results$KEGG[annotation_match],
    Class = metabolite_results$Class[annotation_match],
    SubClass = metabolite_results[["Sub Class"]][annotation_match],
    VIP = metabolite_results$VIP[annotation_match],
    VIPCandidate = metabolite_results$VIPCandidate[annotation_match],
    log2FC = metabolite_results$log2FC[annotation_match],
    PearsonR = correlation,
    CorrelationT = correlation_t,
    Pvalue = pvalue,
    padj = p.adjust(pvalue, method = "BH"),
    N = n_samples,
    Omics = omics_label,
    stringsAsFactors = FALSE
  )
  output <- output[is.finite(output$PearsonR), , drop = FALSE]
  output <- output[order(-abs(output$PearsonR)), , drop = FALSE]
  list(
    table = output,
    subjects = data.frame(
      Subject = shared_subjects,
      TargetSample = target_samples,
      MetaboliteSample = metabolite_samples,
      stringsAsFactors = FALSE
    ),
    target_feature = target_table[[target_id_column]][target_row]
  )
}

rna_metadata <- get_omics_metadata("transcriptomics")
protein_metadata <- get_omics_metadata("proteomics")
metabolite_metadata <- get_omics_metadata("metabolomics")

message("Calculating matched-animal Plcxd2-metabolite correlations")
rna_cor <- calculate_metabolite_correlations(
  rna, "GeneID", "GeneName", rna_metadata, metabolite_metadata,
  "Transcriptomics"
)
protein_cor <- calculate_metabolite_correlations(
  protein, "Accession", "GeneName", protein_metadata, metabolite_metadata,
  "Proteomics"
)
write.csv(
  rna_cor$table,
  file.path(tables_dir, "rna_plcxd2_metabolite_correlations.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  protein_cor$table,
  file.path(tables_dir, "protein_plcxd2_metabolite_correlations.csv"),
  row.names = FALSE, quote = TRUE
)
rna_pairs <- rna_cor$subjects
rna_pairs$Omics <- "Transcriptomics"
protein_pairs <- protein_cor$subjects
protein_pairs$Omics <- "Proteomics"
write.csv(
  rbind(rna_pairs, protein_pairs),
  file.path(tables_dir, "matched_samples_used.csv"),
  row.names = FALSE, quote = TRUE
)

shared <- merge(
  rna_cor$table[, c(
    "FeatureKey", "Metabolite", "KEGG", "Class", "SubClass", "VIP",
    "VIPCandidate", "log2FC", "PearsonR", "Pvalue", "padj", "N"
  )],
  protein_cor$table[, c(
    "FeatureKey", "PearsonR", "Pvalue", "padj", "N"
  )],
  by = "FeatureKey", suffixes = c("_RNA", "_Protein")
)
shared$DirectionPattern <- paste0(
  ifelse(shared$PearsonR_RNA >= 0, "RNA positive", "RNA negative"),
  " / ",
  ifelse(shared$PearsonR_Protein >= 0, "Protein positive", "Protein negative")
)
shared$Concordant <- sign(shared$PearsonR_RNA) == sign(shared$PearsonR_Protein)
shared$MeanAbsR <- rowMeans(abs(shared[, c("PearsonR_RNA", "PearsonR_Protein")]))
shared$MinAbsR <- pmin(abs(shared$PearsonR_RNA), abs(shared$PearsonR_Protein))
shared <- shared[order(!shared$Concordant, -shared$MinAbsR), , drop = FALSE]
write.csv(
  shared,
  file.path(tables_dir, "shared_plcxd2_metabolite_correlations.csv"),
  row.names = FALSE, quote = TRUE
)

lipid_pattern <- paste(
  c("lipid", "phosph", "sphing", "fatty", "acyl", "carnitine", "sterol"),
  collapse = "|"
)
lipid_related <- shared[
  grepl(
    lipid_pattern,
    paste(shared$Metabolite, shared$Class, shared$SubClass),
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]
write.csv(
  lipid_related,
  file.path(tables_dir, "plcxd2_lipid_related_metabolites.csv"),
  row.names = FALSE, quote = TRUE
)

message("Calculating per-animal KEGG pathway activity with ssGSEA")
metabolite_samples_all <- metabolite_metadata$sample_id
if (!all(metabolite_samples_all %in% colnames(metabolite_expression))) {
  stop("Metabolomics sample columns required for ssGSEA are missing.", call. = FALSE)
}
ssgsea_matrix <- as.matrix(
  metabolite_expression[, metabolite_samples_all, drop = FALSE]
)
storage.mode(ssgsea_matrix) <- "numeric"
rownames(ssgsea_matrix) <- metabolite_expression$ID
ssgsea_pathways <- lapply(
  pathway_list,
  function(features) intersect(unique(features), rownames(ssgsea_matrix))
)
ssgsea_pathways <- ssgsea_pathways[
  vapply(ssgsea_pathways, length, integer(1)) >= min_pathway_size &
    vapply(ssgsea_pathways, length, integer(1)) <= max_pathway_size
]
ssgsea_parameter <- GSVA::ssgseaParam(
  exprData = ssgsea_matrix,
  geneSets = ssgsea_pathways,
  minSize = min_pathway_size,
  maxSize = max_pathway_size,
  alpha = 0.25,
  normalize = TRUE
)
pathway_score_matrix <- GSVA::gsva(
  ssgsea_parameter,
  verbose = FALSE,
  BPPARAM = serial_param
)
pathway_score_matrix <- as.matrix(pathway_score_matrix)
pathway_score_table <- data.frame(
  PathwayID = rownames(pathway_score_matrix),
  PathwayName = unname(pathway_name_map[rownames(pathway_score_matrix)]),
  SetSize = vapply(
    ssgsea_pathways[rownames(pathway_score_matrix)], length, integer(1)
  ),
  pathway_score_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(
  pathway_score_table,
  file.path(tables_dir, "kegg_pathway_ssgsea_scores.csv"),
  row.names = FALSE, quote = TRUE
)

calculate_pathway_correlations <- function(
  target_table, target_id_column, target_symbol_column,
  target_metadata, metabolite_metadata, omics_label
) {
  shared_subjects <- intersect(
    target_metadata$subject_id, metabolite_metadata$subject_id
  )
  target_metadata <- target_metadata[
    match(shared_subjects, target_metadata$subject_id), , drop = FALSE
  ]
  metabolite_metadata <- metabolite_metadata[
    match(shared_subjects, metabolite_metadata$subject_id), , drop = FALSE
  ]
  target_samples <- target_metadata$sample_id
  score_samples <- metabolite_metadata$sample_id
  target_rows <- which(
    tolower(target_table[[target_symbol_column]]) == tolower(target_symbol)
  )
  if (length(target_rows) == 0) {
    stop(target_symbol, " was not found in ", omics_label, call. = FALSE)
  }
  target_values_all <- as.matrix(
    target_table[target_rows, target_samples, drop = FALSE]
  )
  storage.mode(target_values_all) <- "numeric"
  target_row <- target_rows[which.max(rowMeans(target_values_all))]
  target_values <- as.numeric(target_table[target_row, target_samples])
  score_matrix <- pathway_score_matrix[, score_samples, drop = FALSE]
  target_centered <- target_values - mean(target_values)
  score_centered <- score_matrix - rowMeans(score_matrix)
  target_ss <- sum(target_centered^2)
  score_ss <- rowSums(score_centered^2)
  correlation <- as.numeric(score_centered %*% target_centered) /
    sqrt(score_ss * target_ss)
  correlation[!is.finite(correlation)] <- NA_real_
  correlation <- pmax(pmin(correlation, 1), -1)
  n_samples <- length(shared_subjects)
  correlation_t <- correlation * sqrt(
    (n_samples - 2) / pmax(1 - correlation^2, .Machine$double.eps)
  )
  pvalue <- 2 * pt(-abs(correlation_t), df = n_samples - 2)
  output <- data.frame(
    PathwayID = rownames(score_matrix),
    PathwayName = unname(pathway_name_map[rownames(score_matrix)]),
    SetSize = vapply(
      ssgsea_pathways[rownames(score_matrix)], length, integer(1)
    ),
    PearsonR = correlation,
    CorrelationT = correlation_t,
    Pvalue = pvalue,
    padj = p.adjust(pvalue, method = "BH"),
    N = n_samples,
    Omics = omics_label,
    stringsAsFactors = FALSE
  )
  output <- output[is.finite(output$PearsonR), , drop = FALSE]
  output <- output[order(-abs(output$PearsonR)), , drop = FALSE]
  list(
    table = output,
    target_feature = target_table[[target_id_column]][target_row]
  )
}

rna_pathway_cor <- calculate_pathway_correlations(
  rna, "GeneID", "GeneName", rna_metadata, metabolite_metadata,
  "Transcriptomics"
)
protein_pathway_cor <- calculate_pathway_correlations(
  protein, "Accession", "GeneName", protein_metadata, metabolite_metadata,
  "Proteomics"
)
write.csv(
  rna_pathway_cor$table,
  file.path(tables_dir, "rna_plcxd2_pathway_correlations.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  protein_pathway_cor$table,
  file.path(tables_dir, "protein_plcxd2_pathway_correlations.csv"),
  row.names = FALSE, quote = TRUE
)

shared_pathways <- merge(
  rna_pathway_cor$table[, c(
    "PathwayID", "PathwayName", "SetSize", "PearsonR", "Pvalue", "padj", "N"
  )],
  protein_pathway_cor$table[, c(
    "PathwayID", "PearsonR", "Pvalue", "padj", "N"
  )],
  by = "PathwayID", suffixes = c("_RNA", "_Protein")
)
shared_pathways$Concordant <-
  sign(shared_pathways$PearsonR_RNA) == sign(shared_pathways$PearsonR_Protein)
shared_pathways$MeanAbsR <- rowMeans(abs(
  shared_pathways[, c("PearsonR_RNA", "PearsonR_Protein")]
))
shared_pathways$MinAbsR <- pmin(
  abs(shared_pathways$PearsonR_RNA),
  abs(shared_pathways$PearsonR_Protein)
)
ora_match <- match(shared_pathways$PathwayID, ora$PathwayID)
shared_pathways$VIPORAFoldEnrichment <- ora$FoldEnrichment[ora_match]
shared_pathways$VIPORAOverlap <- ora$OverlapSize[ora_match]
shared_pathways$VIPORAPvalue <- ora$Pvalue[ora_match]
shared_pathways$VIPORAFDR <- ora$padj[ora_match]
shared_pathways <- shared_pathways[
  order(!shared_pathways$Concordant, -shared_pathways$MinAbsR),
  ,
  drop = FALSE
]
write.csv(
  shared_pathways,
  file.path(tables_dir, "shared_plcxd2_pathway_correlations.csv"),
  row.names = FALSE, quote = TRUE
)

selected_pathway_ids <- head(
  ora$PathwayID[ora$OverlapSize > 0], 6L
)
selected_pathways <- shared_pathways[
  match(selected_pathway_ids, shared_pathways$PathwayID), , drop = FALSE
]
selected_pathways <- selected_pathways[
  !is.na(selected_pathways$PathwayID), , drop = FALSE
]
selected_pathways <- head(selected_pathways, 6L)

pathway_scatter <- shared_pathways
pathway_scatter$Evidence <- ifelse(
  is.finite(pathway_scatter$VIPORAPvalue) &
    pathway_scatter$VIPORAPvalue < 0.05,
  "VIP-enriched (P < 0.05)", "Other"
)
pathway_scatter$Evidence <- factor(
  pathway_scatter$Evidence,
  levels = c("Other", "VIP-enriched (P < 0.05)")
)
pathway_correlation_plot <- ggplot(
  pathway_scatter,
  aes(x = PearsonR_RNA, y = PearsonR_Protein, color = Evidence)
) +
  geom_hline(yintercept = 0, color = "#A5ADB8", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "#A5ADB8", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "#7B8491", linetype = "dashed") +
  geom_point(aes(size = SetSize), alpha = 0.78) +
  geom_text_repel(
    data = pathway_scatter[
      pathway_scatter$PathwayID %in% selected_pathways$PathwayID,
      ,
      drop = FALSE
    ],
    aes(label = short_label(PathwayName, 30)),
    size = 3.0, box.padding = 0.5, point.padding = 0.25,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      Other = "#C8CED6", `VIP-enriched (P < 0.05)` = "#8A5FBF"
    )
  ) +
  scale_size_continuous(range = c(2.5, 7.5)) +
  coord_equal(xlim = c(-1.03, 1.03), ylim = c(-1.03, 1.03)) +
  labs(
    title = "Plcxd2-pathway activity correlations",
    x = "RNA Plcxd2 correlation", y = "Protein Plcxd2 correlation",
    color = "Evidence", size = "Metabolites"
  ) +
  theme_project()
save_ggplot(
  pathway_correlation_plot,
  file.path(figures_dir, "04_plcxd2_pathway_activity_correlations"),
  9.4, 8.2
)

network_member_candidates <- merge(
  membership[
    membership$PathwayID %in% selected_pathways$PathwayID,
    ,
    drop = FALSE
  ],
  shared[, c(
    "FeatureKey", "Metabolite", "VIPCandidate", "MeanAbsR", "MinAbsR"
  )],
  by = "FeatureKey"
)
network_member_rows <- if (nrow(network_member_candidates) > 0) {
  do.call(
    rbind,
    lapply(
      split(network_member_candidates, network_member_candidates$PathwayID),
      function(x) {
        x <- x[order(!x$VIPCandidate, -x$MeanAbsR), , drop = FALSE]
        head(x, 2L)
      }
    )
  )
} else {
  network_member_candidates
}
selected_keys <- unique(network_member_rows$FeatureKey)
selected_network <- shared[
  match(selected_keys, shared$FeatureKey), , drop = FALSE
]
scatter_data <- shared
scatter_data$Evidence <- ifelse(
  scatter_data$VIPCandidate %in% TRUE, "VIP candidate", "Other"
)
scatter_data$Evidence <- factor(
  scatter_data$Evidence,
  levels = c("Other", "VIP candidate")
)
scatter_plot <- ggplot(
  scatter_data,
  aes(x = PearsonR_RNA, y = PearsonR_Protein, color = Evidence)
) +
  geom_hline(yintercept = 0, color = "#A5ADB8", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "#A5ADB8", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "#7B8491", linetype = "dashed") +
  geom_point(alpha = 0.72, size = 2.4) +
  geom_text_repel(
    data = scatter_data[scatter_data$FeatureKey %in% selected_keys, , drop = FALSE],
    aes(label = short_label(Metabolite, 22)),
    size = 3.0, box.padding = 0.45, point.padding = 0.25,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      Other = "#C8CED6", `VIP candidate` = "#F2B134"
    )
  ) +
  coord_equal(xlim = c(-1.03, 1.03), ylim = c(-1.03, 1.03)) +
  labs(
    title = "Plcxd2-metabolite correlation concordance",
    x = "RNA Plcxd2 correlation", y = "Protein Plcxd2 correlation",
    color = "Evidence"
  ) +
  theme_project()
save_ggplot(
  scatter_plot,
  file.path(figures_dir, "03_plcxd2_metabolite_correlation_concordance"),
  9.2, 8.0
)

if (nrow(selected_pathways) > 0 && nrow(network_member_rows) > 0) {
  selected_pathways$Node <- paste0("Pathway_", selected_pathways$PathwayID)
  selected_pathways$Label <- make.unique(wrap_label(
    short_label(selected_pathways$PathwayName, 42), 20
  ))
  selected_pathways$y <- seq(
    from = (nrow(selected_pathways) - 1) / 2,
    to = -(nrow(selected_pathways) - 1) / 2,
    length.out = nrow(selected_pathways)
  )
  pathway_nodes <- data.frame(
    Node = selected_pathways$Node,
    Label = selected_pathways$Label,
    NodeType = "Pathway",
    NodeClass = ifelse(
      selected_pathways$PearsonR_RNA >= 0,
      "Positive pathway", "Negative pathway"
    ),
    x = 0,
    y = selected_pathways$y,
    stringsAsFactors = FALSE
  )

  metabolite_keys <- unique(network_member_rows$FeatureKey)
  metabolite_info <- network_member_rows[
    match(metabolite_keys, network_member_rows$FeatureKey), , drop = FALSE
  ]
  metabolite_nodes <- data.frame(
    Node = metabolite_info$FeatureKey,
    Label = make.unique(wrap_label(
      short_label(metabolite_info$Metabolite, 30), 15
    )),
    NodeType = "Metabolite",
    NodeClass = ifelse(
      metabolite_info$VIPCandidate %in% TRUE,
      "VIP metabolite", "Other metabolite"
    ),
    x = 3.2,
    y = seq(
      from = (length(metabolite_keys) - 1) / 2,
      to = -(length(metabolite_keys) - 1) / 2,
      length.out = length(metabolite_keys)
    ),
    stringsAsFactors = FALSE
  )
  target_nodes <- data.frame(
    Node = c("Plcxd2_RNA", "Plcxd2_Protein"),
    Label = c("Plcxd2\nRNA", "Plcxd2\nprotein"),
    NodeType = c("Transcriptomics", "Proteomics"),
    NodeClass = c("Plcxd2 RNA", "Plcxd2 protein"),
    x = -3.2,
    y = c(1.25, -1.25),
    stringsAsFactors = FALSE
  )
  network_nodes <- rbind(target_nodes, pathway_nodes, metabolite_nodes)

  rna_pathway_edges <- data.frame(
    From = "Plcxd2_RNA", To = selected_pathways$Node,
    EdgeType = "Plcxd2-pathway", Omics = "Transcriptomics",
    Correlation = selected_pathways$PearsonR_RNA,
    Direction = ifelse(
      selected_pathways$PearsonR_RNA >= 0, "Positive", "Negative"
    ),
    x = -3.2, y = 1.25, xend = 0, yend = selected_pathways$y,
    stringsAsFactors = FALSE
  )
  protein_pathway_edges <- data.frame(
    From = "Plcxd2_Protein", To = selected_pathways$Node,
    EdgeType = "Plcxd2-pathway", Omics = "Proteomics",
    Correlation = selected_pathways$PearsonR_Protein,
    Direction = ifelse(
      selected_pathways$PearsonR_Protein >= 0, "Positive", "Negative"
    ),
    x = -3.2, y = -1.25, xend = 0, yend = selected_pathways$y,
    stringsAsFactors = FALSE
  )
  membership_edges <- unique(data.frame(
    From = paste0("Pathway_", network_member_rows$PathwayID),
    To = network_member_rows$FeatureKey,
    EdgeType = "Pathway-metabolite",
    Omics = "Membership",
    Correlation = NA_real_,
    Direction = "Membership",
    stringsAsFactors = FALSE
  ))
  membership_edges$x <- pathway_nodes$x[
    match(membership_edges$From, pathway_nodes$Node)
  ]
  membership_edges$y <- pathway_nodes$y[
    match(membership_edges$From, pathway_nodes$Node)
  ]
  membership_edges$xend <- metabolite_nodes$x[
    match(membership_edges$To, metabolite_nodes$Node)
  ]
  membership_edges$yend <- metabolite_nodes$y[
    match(membership_edges$To, metabolite_nodes$Node)
  ]
  target_pathway_edges <- rbind(rna_pathway_edges, protein_pathway_edges)
  network_edges <- rbind(target_pathway_edges, membership_edges)
  write.csv(
    network_nodes,
    file.path(tables_dir, "plcxd2_pathway_metabolite_network_nodes.csv"),
    row.names = FALSE, quote = TRUE
  )
  write.csv(
    network_edges,
    file.path(tables_dir, "plcxd2_pathway_metabolite_network_edges.csv"),
    row.names = FALSE, quote = TRUE
  )

  network_graph <- igraph::graph_from_data_frame(
    d = network_edges,
    directed = FALSE,
    vertices = network_nodes
  )
  set.seed(20260921)
  network_layout <- ggraph::create_layout(network_graph, layout = "fr")
  target_layout_rows <- match(
    c("Plcxd2_RNA", "Plcxd2_Protein"), network_layout$name
  )
  layout_span_x <- diff(range(network_layout$x, finite = TRUE))
  target_center_x <- mean(network_layout$x[target_layout_rows])
  target_center_y <- mean(network_layout$y[target_layout_rows])
  target_gap <- 0.14 * layout_span_x
  network_layout$x[target_layout_rows] <- target_center_x + c(
    -target_gap / 2, target_gap / 2
  )
  network_layout$y[target_layout_rows] <- target_center_y

  network_plot <- ggraph::ggraph(network_layout) +
    ggraph::geom_edge_link(
      aes(filter = EdgeType == "Pathway-metabolite"),
      edge_colour = "#C5CBD4", edge_width = 0.75,
      edge_alpha = 0.70, lineend = "round", show.legend = FALSE
    ) +
    ggraph::geom_edge_link(
      aes(
        filter = EdgeType == "Plcxd2-pathway",
        edge_colour = Direction,
        edge_linetype = Omics,
        edge_width = abs(Correlation)
      ),
      edge_alpha = 0.82, lineend = "round"
    ) +
    ggraph::geom_node_point(
      aes(fill = NodeClass, size = NodeType),
      shape = 21, stroke = 1.15, color = "white"
    ) +
    ggraph::geom_node_text(
      aes(
        filter = NodeType %in% c("Transcriptomics", "Proteomics"),
        label = Label
      ),
      color = "white", fontface = "bold", size = 3.7,
      lineheight = 0.90
    ) +
    ggraph::geom_node_label(
      aes(filter = NodeType == "Pathway", label = Label),
      repel = TRUE, size = 4.0, color = "#273142",
      fill = scales::alpha("white", 0.92), label.size = 0,
      label.padding = grid::unit(0.14, "lines"),
      box.padding = 0.55, point.padding = 0.45,
      max.overlaps = Inf, segment.color = "#AAB2BE",
      segment.size = 0.35
    ) +
    ggraph::geom_node_label(
      aes(filter = NodeType == "Metabolite", label = Label),
      repel = TRUE, size = 3.5, color = "#273142",
      fill = scales::alpha("white", 0.90), label.size = 0,
      label.padding = grid::unit(0.12, "lines"),
      box.padding = 0.45, point.padding = 0.38,
      max.overlaps = Inf, segment.color = "#B8C0CC",
      segment.size = 0.30
    ) +
    ggraph::scale_edge_colour_manual(values = direction_colors) +
    ggraph::scale_edge_linetype_manual(
      values = c(Transcriptomics = "solid", Proteomics = "longdash")
    ) +
    ggraph::scale_edge_width_continuous(range = c(0.7, 2.6), limits = c(0, 1)) +
    scale_size_manual(
      values = c(
        Transcriptomics = 22.0, Proteomics = 22.0,
        Pathway = 11.5, Metabolite = 8.8
      ),
      guide = "none"
    ) +
    scale_fill_manual(
      values = c(
        `Plcxd2 RNA` = omics_colors[["Transcriptomics"]],
        `Plcxd2 protein` = omics_colors[["Proteomics"]],
        `Positive pathway` = "#D95F59", `Negative pathway` = "#4C78A8",
        `VIP metabolite` = "#F2B134", `Other metabolite` = "#59A14F"
      )
    ) +
    labs(
      title = "Plcxd2-pathway-metabolite network",
      edge_colour = "Correlation", edge_linetype = "Omics",
      edge_width = "|Pearson r|", fill = "Node"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    scale_y_continuous(expand = expansion(mult = 0.14)) +
    coord_cartesian(clip = "off") +
    ggraph::theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(
        face = "bold", size = 18, color = "#111827", hjust = 0.5
      ),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(22, 28, 18, 28)
    )
  save_ggplot(
    network_plot,
    file.path(figures_dir, "05_plcxd2_pathway_metabolite_network"),
    14.0, 10.6
  )
}

network_metabolite_total <- if (exists("metabolite_nodes")) {
  nrow(metabolite_nodes)
} else {
  0L
}
summary_table <- data.frame(
  Metric = c(
    "Metabolite features", "KEGG compound-annotated features",
    "KEGG pathway-annotated features", "VIP candidates",
    "Pathway-annotated VIP candidates", "Tested KEGG pathways",
    "VIP ORA pathways at FDR < 0.05",
    "RNA-metabolomics matched animals",
    "Protein-metabolomics matched animals",
    "Concordant Plcxd2-metabolite correlations",
    "ssGSEA-scored KEGG pathways",
    "Concordant Plcxd2-pathway correlations",
    "Network pathways", "Network metabolites"
  ),
  Value = c(
    nrow(metabolite_results), sum(!is.na(metabolite_results$KEGG)),
    length(background_features), sum(metabolite_results$VIPCandidate %in% TRUE),
    length(candidate_features), nrow(ora),
    sum(ora$padj < 0.05, na.rm = TRUE),
    nrow(rna_cor$subjects), nrow(protein_cor$subjects),
    sum(shared$Concordant), nrow(pathway_score_matrix),
    sum(shared_pathways$Concordant), nrow(selected_pathways),
    network_metabolite_total
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary_table, file.path(output_dir, "analysis_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

parameters <- data.frame(
  Parameter = c(
    "Target", "Correlation", "Group adjustment", "RNA matched animals",
    "Protein matched animals", "Pathway source", "ORA candidate definition",
    "Pathway activity method", "ssGSEA normalization",
    "Minimum pathway size", "Maximum pathway size"
  ),
  Value = c(
    target_symbol, "Pearson", "None", nrow(rna_cor$subjects),
    nrow(protein_cor$subjects), "Local KEGG annotations from metabolomics table",
    "OPLS-DA VIP candidate", "GSVA::ssGSEA (alpha = 0.25)", "Enabled",
    min_pathway_size, max_pathway_size
  ),
  stringsAsFactors = FALSE
)
write.table(
  parameters, file.path(output_dir, "analysis_parameters.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

message(
  "Plcxd2-metabolomics integration complete. Pathway-annotated features: ",
  length(background_features), "; RNA pairs: ", nrow(rna_cor$subjects),
  "; protein pairs: ", nrow(protein_cor$subjects), "."
)
