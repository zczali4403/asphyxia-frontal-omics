#!/usr/bin/env Rscript

# Plcxd2-centered co-response analysis and ranked GSEA.
# Pearson correlations are calculated across all samples in each omics layer
# without adjustment for Control/Asphyxia group.

required_packages <- c(
  "ggplot2", "clusterProfiler", "AnnotationDbi", "org.Mm.eg.db",
  "GO.db", "KEGGREST", "BiocParallel"
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
  library(ggplot2)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
})
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

args <- commandArgs(trailingOnly = TRUE)
rna_file <- if (length(args) >= 1) args[[1]] else
  "results/01_transcriptomics/tables/vst_expression_filtered.csv"
protein_file <- if (length(args) >= 2) args[[2]] else
  "results/02_proteomics/tables/protein_expression_filtered.csv"
metadata_file <- if (length(args) >= 3) args[[3]] else "config/samples.csv"
output_dir <- if (length(args) >= 4) args[[4]] else
  "results/04_plcxd2_correlation_gsea"
reference_dir <- if (length(args) >= 5) args[[5]] else "data/reference"

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

target_symbol <- "Plcxd2"
nominal_threshold <- 0.05
fdr_threshold <- 0.05
min_gene_set_size <- 10L
max_gene_set_size <- 500L
plot_pathways_per_direction <- 8L
max_integrated_plot_pathways <- 30L
set.seed(20260920)

message("Reading expression matrices and metadata")
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

required_metadata <- c("omics", "sample_id", "group", "region", "sample_role")
if (!all(required_metadata %in% colnames(metadata))) {
  stop("The metadata table is missing required columns.", call. = FALSE)
}
if (!all(c("GeneID", "GeneName") %in% colnames(rna))) {
  stop("The RNA expression table must contain GeneID and GeneName.", call. = FALSE)
}
if (!all(c("Accession", "GeneName") %in% colnames(protein))) {
  stop(
    "The protein expression table must contain Accession and GeneName.",
    call. = FALSE
  )
}

rna_metadata <- metadata[
  metadata$omics == "transcriptomics" &
    metadata$sample_role == "biological" &
    metadata$region == "frontal_cortex",
  ,
  drop = FALSE
]
protein_metadata <- metadata[
  metadata$omics == "proteomics" &
    metadata$sample_role == "biological" &
    metadata$region == "frontal_cortex",
  ,
  drop = FALSE
]
rna_samples <- rna_metadata$sample_id
protein_samples <- protein_metadata$sample_id
if (!all(rna_samples %in% colnames(rna))) {
  stop("RNA samples in metadata are missing from the expression table.")
}
if (!all(protein_samples %in% colnames(protein))) {
  stop("Protein samples in metadata are missing from the expression table.")
}

calculate_target_correlations <- function(
  expression_table,
  id_column,
  symbol_column,
  sample_columns,
  target
) {
  matrix_data <- as.matrix(expression_table[, sample_columns, drop = FALSE])
  storage.mode(matrix_data) <- "numeric"
  if (anyNA(matrix_data)) {
    stop("The correlation input matrix contains missing values.", call. = FALSE)
  }
  target_rows <- which(tolower(expression_table[[symbol_column]]) == tolower(target))
  if (length(target_rows) == 0) {
    stop(target, " was not found in the expression matrix.", call. = FALSE)
  }
  if (length(target_rows) > 1) {
    target_rows <- target_rows[
      which.max(rowMeans(matrix_data[target_rows, , drop = FALSE]))
    ]
  }
  target_values <- matrix_data[target_rows, ]
  target_centered <- target_values - mean(target_values)
  target_ss <- sum(target_centered^2)
  if (target_ss <= 0) {
    stop(target, " has no expression variance.", call. = FALSE)
  }

  centered <- matrix_data - rowMeans(matrix_data)
  feature_ss <- rowSums(centered^2)
  correlation <- as.numeric(centered %*% target_centered) /
    sqrt(feature_ss * target_ss)
  correlation[!is.finite(correlation)] <- NA_real_
  correlation <- pmax(pmin(correlation, 1), -1)
  n_samples <- length(sample_columns)
  correlation_t <- correlation * sqrt(
    (n_samples - 2) / pmax(1 - correlation^2, .Machine$double.eps)
  )
  correlation_p <- 2 * pt(-abs(correlation_t), df = n_samples - 2)

  output <- data.frame(
    FeatureID = expression_table[[id_column]],
    GeneName = expression_table[[symbol_column]],
    PearsonR = correlation,
    CorrelationT = correlation_t,
    Pvalue = correlation_p,
    MeanExpression = rowMeans(matrix_data),
    stringsAsFactors = FALSE
  )
  output <- output[-target_rows, , drop = FALSE]
  output <- output[is.finite(output$PearsonR), , drop = FALSE]
  output <- output[order(-output$PearsonR), , drop = FALSE]
  list(
    table = output,
    target_values = target_values,
    target_feature_id = expression_table[[id_column]][target_rows],
    sample_count = n_samples
  )
}

rna_cor <- calculate_target_correlations(
  rna, "GeneID", "GeneName", rna_samples, target_symbol
)
protein_cor <- calculate_target_correlations(
  protein, "Accession", "GeneName", protein_samples, target_symbol
)

rna_ensembl <- sub("\\..*$", "", rna_cor$table$FeatureID)
rna_entrez_map <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(rna_ensembl),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rna_cor$table$EntrezID <- unname(
  rna_entrez_map[match(rna_ensembl, names(rna_entrez_map))]
)

protein_symbols <- protein_cor$table$GeneName
protein_entrez_map <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(protein_symbols[!is.na(protein_symbols) & protein_symbols != ""]),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
protein_cor$table$EntrezID <- unname(
  protein_entrez_map[match(protein_symbols, names(protein_entrez_map))]
)

write.csv(
  rna_cor$table,
  file.path(tables_dir, "rna_plcxd2_correlations.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  protein_cor$table,
  file.path(tables_dir, "protein_plcxd2_correlations.csv"),
  row.names = FALSE,
  quote = TRUE
)

collapse_correlation_rank <- function(correlation_table) {
  ranked <- correlation_table[
    !is.na(correlation_table$EntrezID) &
      correlation_table$EntrezID != "" &
      is.finite(correlation_table$PearsonR),
    c("EntrezID", "PearsonR", "MeanExpression"),
    drop = FALSE
  ]
  ranked <- ranked[
    order(ranked$EntrezID, -abs(ranked$PearsonR), -ranked$MeanExpression),
    ,
    drop = FALSE
  ]
  ranked <- ranked[!duplicated(ranked$EntrezID), , drop = FALSE]
  values <- ranked$PearsonR
  names(values) <- ranked$EntrezID
  sort(values, decreasing = TRUE)
}

rna_rank <- collapse_correlation_rank(rna_cor$table)
protein_rank <- collapse_correlation_rank(protein_cor$table)
write.csv(
  data.frame(EntrezID = names(rna_rank), PearsonR = unname(rna_rank)),
  file.path(tables_dir, "rna_correlation_rank.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  data.frame(EntrezID = names(protein_rank), PearsonR = unname(protein_rank)),
  file.path(tables_dir, "protein_correlation_rank.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Preparing GO BP, MF, CC, and KEGG gene sets")
annotated_entrez <- unique(c(names(rna_rank), names(protein_rank)))
go_annotations <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db,
  keys = annotated_entrez,
  columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "ENTREZID"
))
go_annotations <- unique(go_annotations[
  !is.na(go_annotations$GOALL) &
    go_annotations$ONTOLOGYALL %in% c("BP", "MF", "CC"),
  c("GOALL", "ENTREZID", "ONTOLOGYALL"),
  drop = FALSE
])
colnames(go_annotations) <- c("Term", "Gene", "Ontology")
go_names <- suppressMessages(AnnotationDbi::select(
  GO.db,
  keys = unique(go_annotations$Term),
  columns = "TERM",
  keytype = "GOID"
))
go_names <- unique(go_names[!is.na(go_names$TERM), , drop = FALSE])
colnames(go_names) <- c("Term", "Name")

kegg_gene_file <- file.path(reference_dir, "kegg_mouse_gene_pathways.tsv")
kegg_name_file <- file.path(reference_dir, "kegg_mouse_pathway_names.tsv")
if (!file.exists(kegg_gene_file) || !file.exists(kegg_name_file)) {
  message("Downloading mouse KEGG mappings for the local cache")
  kegg_links <- KEGGREST::keggLink("pathway", "mmu")
  kegg_pathways <- KEGGREST::keggList("pathway", "mmu")
  kegg_gene_map <- unique(data.frame(
    Term = sub("^path:", "", unname(kegg_links)),
    Gene = sub("^mmu:", "", names(kegg_links)),
    stringsAsFactors = FALSE
  ))
  kegg_name_map <- data.frame(
    Term = names(kegg_pathways),
    Name = sub(
      " - Mus musculus \\(house mouse\\)$", "", unname(kegg_pathways)
    ),
    stringsAsFactors = FALSE
  )
  write.table(
    kegg_gene_map, kegg_gene_file, sep = "\t",
    row.names = FALSE, quote = FALSE
  )
  write.table(
    kegg_name_map, kegg_name_file, sep = "\t",
    row.names = FALSE, quote = FALSE
  )
} else {
  kegg_gene_map <- read.delim(
    kegg_gene_file, check.names = FALSE, stringsAsFactors = FALSE
  )
  kegg_name_map <- read.delim(
    kegg_name_file, check.names = FALSE, stringsAsFactors = FALSE
  )
}

run_gsea <- function(rank_vector, term_to_gene, term_to_name, database, omics) {
  set.seed(20260920)
  result <- suppressWarnings(clusterProfiler::GSEA(
    geneList = rank_vector,
    exponent = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = term_to_gene,
    TERM2NAME = term_to_name,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  ))
  output <- as.data.frame(result)
  if (nrow(output) == 0) {
    return(output)
  }
  output$Database <- database
  output$Omics <- omics
  output$Direction <- ifelse(output$NES > 0, "Positive", "Negative")
  output$PathwayKey <- paste(database, output$ID, sep = ":")
  output
}

run_all_go_gsea <- function(rank_vector, omics) {
  do.call(
    rbind,
    lapply(c("BP", "MF", "CC"), function(ontology) {
      ontology_map <- go_annotations[
        go_annotations$Ontology == ontology,
        c("Term", "Gene"),
        drop = FALSE
      ]
      run_gsea(
        rank_vector,
        ontology_map,
        go_names,
        paste("GO", ontology),
        omics
      )
    })
  )
}

message("Running Plcxd2 correlation-ranked GSEA")
rna_gsea <- rbind(
  run_all_go_gsea(rna_rank, "Transcriptomics"),
  run_gsea(rna_rank, kegg_gene_map, kegg_name_map, "KEGG", "Transcriptomics")
)
protein_gsea <- rbind(
  run_all_go_gsea(protein_rank, "Proteomics"),
  run_gsea(
    protein_rank, kegg_gene_map, kegg_name_map, "KEGG", "Proteomics"
  )
)
rna_gsea <- rna_gsea[order(rna_gsea$p.adjust, rna_gsea$pvalue), , drop = FALSE]
protein_gsea <- protein_gsea[
  order(protein_gsea$p.adjust, protein_gsea$pvalue), , drop = FALSE
]

write.csv(
  rna_gsea, file.path(tables_dir, "rna_correlation_gsea_all.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  protein_gsea, file.path(tables_dir, "protein_correlation_gsea_all.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  rna_gsea[rna_gsea$pvalue < nominal_threshold, , drop = FALSE],
  file.path(tables_dir, "rna_correlation_gsea_pvalue.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  protein_gsea[protein_gsea$pvalue < nominal_threshold, , drop = FALSE],
  file.path(tables_dir, "protein_correlation_gsea_pvalue.csv"),
  row.names = FALSE, quote = TRUE
)

rna_integration <- rna_gsea[, c(
  "PathwayKey", "Database", "ID", "Description", "setSize", "NES",
  "pvalue", "p.adjust", "core_enrichment"
)]
protein_integration <- protein_gsea[, c(
  "PathwayKey", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment"
)]
colnames(rna_integration)[5:9] <- c(
  "RNA_setSize", "RNA_NES", "RNA_pvalue", "RNA_padj", "RNA_leadingEdge"
)
colnames(protein_integration)[2:6] <- c(
  "Protein_setSize", "Protein_NES", "Protein_pvalue", "Protein_padj",
  "Protein_leadingEdge"
)
integrated <- merge(
  rna_integration, protein_integration,
  by = "PathwayKey", all = FALSE, sort = FALSE
)
integrated$DirectionPattern <- ifelse(
  integrated$RNA_NES > 0 & integrated$Protein_NES > 0,
  "Both positive",
  ifelse(
    integrated$RNA_NES < 0 & integrated$Protein_NES < 0,
    "Both negative",
    "Opposite"
  )
)
integrated$SharedPvalue <- integrated$RNA_pvalue < nominal_threshold &
  integrated$Protein_pvalue < nominal_threshold
integrated$SharedFDR <- integrated$RNA_padj < fdr_threshold &
  integrated$Protein_padj < fdr_threshold
integrated$JointPvalue <- pmax(
  integrated$RNA_pvalue, integrated$Protein_pvalue
)
integrated <- integrated[
  order(!integrated$SharedPvalue, integrated$JointPvalue),
  ,
  drop = FALSE
]
shared_concordant <- integrated[
  integrated$SharedPvalue & integrated$DirectionPattern != "Opposite",
  ,
  drop = FALSE
]
write.csv(
  integrated,
  file.path(tables_dir, "integrated_correlation_gsea.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  shared_concordant,
  file.path(tables_dir, "shared_concordant_correlation_pathways.csv"),
  row.names = FALSE,
  quote = TRUE
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

wrap_label <- function(x, width = 46) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

select_overview <- function(gsea_table) {
  candidates <- gsea_table[gsea_table$pvalue < nominal_threshold, , drop = FALSE]
  if (nrow(candidates) == 0) {
    candidates <- gsea_table
  }
  positive <- candidates[candidates$NES > 0, , drop = FALSE]
  negative <- candidates[candidates$NES < 0, , drop = FALSE]
  positive <- head(
    positive[order(positive$pvalue, -positive$NES), , drop = FALSE],
    plot_pathways_per_direction
  )
  negative <- head(
    negative[order(negative$pvalue, negative$NES), , drop = FALSE],
    plot_pathways_per_direction
  )
  unique(rbind(negative, positive))
}

make_overview_plot <- function(gsea_table, omics_label) {
  plot_data <- select_overview(gsea_table)
  plot_data$PathwayLabel <- paste0(
    "[", plot_data$Database, "] ", wrap_label(plot_data$Description)
  )
  plot_data$PathwayLabel <- factor(
    plot_data$PathwayLabel,
    levels = plot_data$PathwayLabel[order(plot_data$NES)]
  )
  plot_data$PScore <- -log10(pmax(plot_data$pvalue, 1e-10))
  ggplot(plot_data, aes(x = NES, y = PathwayLabel)) +
    geom_vline(xintercept = 0, color = "#9AA2AE", linewidth = 0.45) +
    geom_segment(
      aes(x = 0, xend = NES, yend = PathwayLabel, color = NES),
      linewidth = 0.8, alpha = 0.55
    ) +
    geom_point(aes(color = NES, size = PScore), alpha = 0.95) +
    scale_color_gradient2(
      low = "#3B6FB6", mid = "#D5DAE1", high = "#D95555", midpoint = 0
    ) +
    scale_size_continuous(range = c(2.6, 6.0)) +
    labs(
      title = paste0("Plcxd2-associated ", tolower(omics_label), " GSEA"),
      subtitle = "P value < 0.05",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = "NES",
      size = "-log10(P value)"
    ) +
    theme_project() +
    theme(panel.grid.major.y = element_blank())
}

rna_plot <- make_overview_plot(rna_gsea, "Transcriptomic")
protein_plot <- make_overview_plot(protein_gsea, "Proteomic")
save_ggplot(
  rna_plot, file.path(figures_dir, "01_rna_correlation_gsea"), 9.2, 7.2
)
save_ggplot(
  protein_plot,
  file.path(figures_dir, "02_protein_correlation_gsea"),
  9.2,
  7.2
)

if (nrow(shared_concordant) > 0) {
  integrated_negative <- shared_concordant[
    shared_concordant$DirectionPattern == "Both negative",
    ,
    drop = FALSE
  ]
  integrated_positive <- shared_concordant[
    shared_concordant$DirectionPattern == "Both positive",
    ,
    drop = FALSE
  ]
  pathways_per_direction <- max_integrated_plot_pathways %/% 2L
  integrated_selected <- rbind(
    head(
      integrated_negative[order(integrated_negative$JointPvalue), , drop = FALSE],
      pathways_per_direction
    ),
    head(
      integrated_positive[order(integrated_positive$JointPvalue), , drop = FALSE],
      pathways_per_direction
    )
  )
  integrated_selected$PathwayLabel <- paste0(
    "[", integrated_selected$Database, "] ",
    wrap_label(integrated_selected$Description)
  )
  integrated_selected$PathwayLabel <- factor(
    integrated_selected$PathwayLabel,
    levels = rev(integrated_selected$PathwayLabel)
  )
  integrated_long <- rbind(
    data.frame(
      PathwayLabel = integrated_selected$PathwayLabel,
      Omics = "Transcriptomics",
      NES = integrated_selected$RNA_NES,
      stringsAsFactors = FALSE
    ),
    data.frame(
      PathwayLabel = integrated_selected$PathwayLabel,
      Omics = "Proteomics",
      NES = integrated_selected$Protein_NES,
      stringsAsFactors = FALSE
    )
  )
  integrated_long$Omics <- factor(
    integrated_long$Omics,
    levels = c("Transcriptomics", "Proteomics")
  )
  integrated_long$PathwayLabel <- factor(
    integrated_long$PathwayLabel,
    levels = levels(integrated_selected$PathwayLabel)
  )

  dotplot <- ggplot() +
    geom_vline(xintercept = 0, color = "#9AA2AE", linewidth = 0.45) +
    geom_segment(
      data = integrated_selected,
      aes(
        x = RNA_NES, xend = Protein_NES,
        y = PathwayLabel, yend = PathwayLabel
      ),
      color = "#C3C8D0", linewidth = 0.8
    ) +
    geom_point(
      data = integrated_long,
      aes(x = NES, y = PathwayLabel, color = Omics),
      size = 3.4
    ) +
    scale_color_manual(
      values = c(Transcriptomics = "#6B4FA3", Proteomics = "#2E8B68")
    ) +
    labs(
      title = "Shared Plcxd2-associated pathways",
      subtitle = "concordant pathways | P value < 0.05",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = "Omics"
    ) +
    theme_project() +
    theme(panel.grid.major.y = element_blank())
  save_ggplot(
    dotplot,
    file.path(figures_dir, "03_shared_concordant_pathways"),
    10.0,
    max(6.4, 2.8 + 0.34 * nrow(integrated_selected))
  )

  heatmap <- ggplot(
    integrated_long,
    aes(x = Omics, y = PathwayLabel, fill = NES)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_gradient2(
      low = "#3B6FB6", mid = "#F4F5F7", high = "#D95555", midpoint = 0
    ) +
    labs(
      title = "Shared Plcxd2-associated pathways",
      subtitle = "concordant pathways | P value < 0.05",
      x = NULL,
      y = NULL,
      fill = "NES"
    ) +
    theme_project() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold")
    )
  save_ggplot(
    heatmap,
    file.path(figures_dir, "04_shared_concordant_pathway_heatmap"),
    8.6,
    max(6.4, 2.8 + 0.34 * nrow(integrated_selected))
  )
}

summary_table <- data.frame(
  Metric = c(
    "RNA samples",
    "Protein samples",
    "RNA features with finite correlation",
    "Protein features with finite correlation",
    "Ranked RNA genes",
    "Ranked protein genes",
    "RNA pathways with P value < 0.05",
    "Protein pathways with P value < 0.05",
    "Shared pathways with P value < 0.05",
    "Shared concordant pathways",
    "Shared opposite-direction pathways",
    "Shared concordant FDR pathways"
  ),
  Value = c(
    rna_cor$sample_count,
    protein_cor$sample_count,
    nrow(rna_cor$table),
    nrow(protein_cor$table),
    length(rna_rank),
    length(protein_rank),
    sum(rna_gsea$pvalue < nominal_threshold),
    sum(protein_gsea$pvalue < nominal_threshold),
    sum(integrated$SharedPvalue),
    nrow(shared_concordant),
    sum(integrated$SharedPvalue & integrated$DirectionPattern == "Opposite"),
    sum(shared_concordant$SharedFDR)
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
write.table(
  data.frame(
    Parameter = c(
      "Target gene", "Correlation method", "Group adjustment",
      "RNA samples", "Protein samples", "GSEA ranking metric"
    ),
    Value = c(
      target_symbol, "Pearson", "None", rna_cor$sample_count,
      protein_cor$sample_count, "Pearson correlation coefficient"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "analysis_parameters.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

message(
  "Plcxd2 correlation GSEA complete. Shared P-value pathways: ",
  sum(integrated$SharedPvalue),
  "; concordant: ", nrow(shared_concordant),
  "; opposite: ",
  sum(integrated$SharedPvalue & integrated$DirectionPattern == "Opposite"),
  "."
)
message("Results written to: ", normalizePath(output_dir))
