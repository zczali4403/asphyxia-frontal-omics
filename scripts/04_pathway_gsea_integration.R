#!/usr/bin/env Rscript

# Ranked GSEA and pathway-level transcriptome-proteome integration.
# Comparison: Asphyxia vs Control in the frontal cortex.

required_packages <- c(
  "ggplot2", "clusterProfiler", "fgsea", "AnnotationDbi",
  "org.Mm.eg.db", "GO.db", "KEGGREST", "BiocParallel", "patchwork"
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
  library(patchwork)
})
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

args <- commandArgs(trailingOnly = TRUE)
rna_file <- if (length(args) >= 1) args[[1]] else
  "results/01_transcriptomics/tables/deseq2_all_genes.csv"
protein_file <- if (length(args) >= 2) args[[2]] else
  "results/02_proteomics/tables/limma_all_proteins.csv"
output_dir <- if (length(args) >= 3) args[[3]] else
  "results/04_transcriptome_proteome_pathways"
reference_dir <- if (length(args) >= 4) args[[4]] else "data/reference"

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

fdr_threshold <- 0.05
nominal_threshold <- 0.05
min_gene_set_size <- 10L
max_gene_set_size <- 500L
plot_pathways_per_direction <- 8L
set.seed(20260920)

required_rna <- c("GeneID", "GeneName", "baseMean", "stat")
required_protein <- c("Accession", "GeneName", "GeneID", "AveExpr", "t")

message("Reading ranked differential-expression results")
rna <- read.csv(
  rna_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
protein <- read.csv(
  protein_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
if (!all(required_rna %in% colnames(rna))) {
  stop("The RNA result table is missing required columns.", call. = FALSE)
}
if (!all(required_protein %in% colnames(protein))) {
  stop("The protein result table is missing required columns.", call. = FALSE)
}

collapse_rank <- function(ids, statistics, abundance) {
  ranked <- data.frame(
    EntrezID = as.character(ids),
    Statistic = as.numeric(statistics),
    Abundance = as.numeric(abundance),
    stringsAsFactors = FALSE
  )
  ranked <- ranked[
    !is.na(ranked$EntrezID) & ranked$EntrezID != "" &
      is.finite(ranked$Statistic),
    ,
    drop = FALSE
  ]
  ranked$Abundance[!is.finite(ranked$Abundance)] <- -Inf
  ranked <- ranked[order(-ranked$Abundance), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$EntrezID), , drop = FALSE]
  values <- ranked$Statistic
  names(values) <- ranked$EntrezID
  sort(values, decreasing = TRUE)
}

rna_ensembl <- sub("\\..*$", "", rna$GeneID)
rna_entrez <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = unique(rna_ensembl),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rna_entrez <- unname(rna_entrez[match(rna_ensembl, names(rna_entrez))])

protein_entrez <- as.character(protein$GeneID)
invalid_protein_id <- is.na(protein_entrez) |
  !grepl("^[0-9]+$", protein_entrez)
if (any(invalid_protein_id)) {
  symbol_entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unique(protein$GeneName[invalid_protein_id]),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  protein_entrez[invalid_protein_id] <- unname(
    symbol_entrez[
      match(protein$GeneName[invalid_protein_id], names(symbol_entrez))
    ]
  )
}

rna_rank <- collapse_rank(rna_entrez, rna$stat, rna$baseMean)
protein_rank <- collapse_rank(protein_entrez, protein$t, protein$AveExpr)

write.csv(
  data.frame(EntrezID = names(rna_rank), Statistic = unname(rna_rank)),
  file.path(tables_dir, "rna_ranked_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  data.frame(EntrezID = names(protein_rank), Statistic = unname(protein_rank)),
  file.path(tables_dir, "protein_ranked_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Preparing local GO Biological Process gene sets")
annotated_entrez <- unique(c(names(rna_rank), names(protein_rank)))
go_annotations <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db,
  keys = annotated_entrez,
  columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "ENTREZID"
))
go_annotations <- unique(go_annotations[
  !is.na(go_annotations$GOALL) & go_annotations$ONTOLOGYALL == "BP",
  c("GOALL", "ENTREZID"),
  drop = FALSE
])
colnames(go_annotations) <- c("Term", "Gene")
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
kegg_info_file <- file.path(reference_dir, "kegg_mouse_cache_info.tsv")
if (!file.exists(kegg_gene_file) || !file.exists(kegg_name_file)) {
  message("Downloading mouse KEGG pathway mappings for the local cache")
  kegg_links <- tryCatch(
    KEGGREST::keggLink("pathway", "mmu"),
    error = function(e) {
      stop(
        "KEGG mappings are not cached and could not be downloaded: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  kegg_pathways <- KEGGREST::keggList("pathway", "mmu")
  kegg_gene_map <- unique(data.frame(
    Term = sub("^path:", "", unname(kegg_links)),
    Gene = sub("^mmu:", "", names(kegg_links)),
    stringsAsFactors = FALSE
  ))
  kegg_name_map <- data.frame(
    Term = names(kegg_pathways),
    Name = sub(
      " - Mus musculus \\(house mouse\\)$",
      "",
      unname(kegg_pathways)
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
  write.table(
    data.frame(
      Source = "KEGG REST API",
      Species = "Mus musculus",
      RetrievalDate = as.character(Sys.Date()),
      stringsAsFactors = FALSE
    ),
    kegg_info_file, sep = "\t", row.names = FALSE, quote = FALSE
  )
} else {
  message("Using cached mouse KEGG pathway mappings")
}
kegg_gene_map <- read.delim(
  kegg_gene_file, check.names = FALSE, stringsAsFactors = FALSE
)
kegg_name_map <- read.delim(
  kegg_name_file, check.names = FALSE, stringsAsFactors = FALSE
)

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
  output$Direction <- ifelse(output$NES > 0, "Up", "Down")
  output$Evidence <- ifelse(
    output$pvalue < nominal_threshold,
    "Nominal P < 0.05",
    "Not significant"
  )
  output$FDRSignificant <- output$p.adjust < fdr_threshold
  output$PathwayKey <- paste(database, output$ID, sep = ":")
  output
}

message("Running ranked GSEA for RNA and protein")
rna_go <- run_gsea(rna_rank, go_annotations, go_names, "GO BP", "Transcriptomics")
protein_go <- run_gsea(
  protein_rank, go_annotations, go_names, "GO BP", "Proteomics"
)
rna_kegg <- run_gsea(
  rna_rank, kegg_gene_map, kegg_name_map, "KEGG", "Transcriptomics"
)
protein_kegg <- run_gsea(
  protein_rank, kegg_gene_map, kegg_name_map, "KEGG", "Proteomics"
)

rna_gsea <- rbind(rna_go, rna_kegg)
protein_gsea <- rbind(protein_go, protein_kegg)
rna_gsea <- rna_gsea[order(rna_gsea$p.adjust, rna_gsea$pvalue), , drop = FALSE]
protein_gsea <- protein_gsea[
  order(protein_gsea$p.adjust, protein_gsea$pvalue), , drop = FALSE
]

write.csv(
  rna_gsea,
  file.path(tables_dir, "rna_gsea_all_pathways.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  protein_gsea,
  file.path(tables_dir, "protein_gsea_all_pathways.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  rna_gsea[rna_gsea$p.adjust < fdr_threshold, , drop = FALSE],
  file.path(tables_dir, "rna_gsea_fdr_significant.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  protein_gsea[protein_gsea$p.adjust < fdr_threshold, , drop = FALSE],
  file.path(tables_dir, "protein_gsea_fdr_significant.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  rna_gsea[rna_gsea$pvalue < nominal_threshold, , drop = FALSE],
  file.path(tables_dir, "rna_gsea_nominal_significant.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  protein_gsea[protein_gsea$pvalue < nominal_threshold, , drop = FALSE],
  file.path(tables_dir, "protein_gsea_nominal_significant.csv"),
  row.names = FALSE,
  quote = TRUE
)

rna_integrated <- rna_gsea[, c(
  "PathwayKey", "Database", "ID", "Description", "setSize", "NES",
  "pvalue", "p.adjust", "core_enrichment"
)]
protein_integrated <- protein_gsea[, c(
  "PathwayKey", "Database", "ID", "Description", "setSize", "NES",
  "pvalue", "p.adjust", "core_enrichment"
)]
colnames(rna_integrated)[5:9] <- c(
  "RNA_setSize", "RNA_NES", "RNA_pvalue", "RNA_padj", "RNA_leadingEdge"
)
colnames(protein_integrated)[5:9] <- c(
  "Protein_setSize", "Protein_NES", "Protein_pvalue", "Protein_padj",
  "Protein_leadingEdge"
)
protein_integrated$Description <- NULL
protein_integrated$Database <- NULL
protein_integrated$ID <- NULL
integrated <- merge(
  rna_integrated,
  protein_integrated,
  by = "PathwayKey",
  all = FALSE,
  sort = FALSE
)
integrated$DirectionPattern <- ifelse(
  integrated$RNA_NES > 0 & integrated$Protein_NES > 0,
  "Both up",
  ifelse(
    integrated$RNA_NES < 0 & integrated$Protein_NES < 0,
    "Both down",
    "Opposite"
  )
)
integrated$RNANominal <- integrated$RNA_pvalue < nominal_threshold
integrated$ProteinNominal <- integrated$Protein_pvalue < nominal_threshold
integrated$RNAFDR <- integrated$RNA_padj < fdr_threshold
integrated$ProteinFDR <- integrated$Protein_padj < fdr_threshold
integrated$EvidencePattern <- ifelse(
  integrated$RNANominal & integrated$ProteinNominal,
  "Both nominal P",
  ifelse(
    integrated$RNANominal,
    "RNA nominal P only",
    ifelse(
      integrated$ProteinNominal,
      "Protein nominal P only",
      "Neither nominal P"
    )
  )
)
integrated$JointNominalP <- pmax(
  integrated$RNA_pvalue,
  integrated$Protein_pvalue
)
integrated$EvidenceRank <- match(
  integrated$EvidencePattern,
  c(
    "Both nominal P", "RNA nominal P only",
    "Protein nominal P only", "Neither nominal P"
  )
)
integrated <- integrated[
  order(integrated$EvidenceRank, integrated$JointNominalP),
  ,
  drop = FALSE
]
write.csv(
  integrated,
  file.path(tables_dir, "integrated_pathway_results.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  integrated[integrated$EvidencePattern != "Neither nominal P", , drop = FALSE],
  file.path(tables_dir, "prioritized_integrated_pathways.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  integrated[integrated$EvidencePattern == "Both nominal P", , drop = FALSE],
  file.path(tables_dir, "shared_nominal_pathways.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  integrated[
    integrated$EvidencePattern == "Both nominal P" &
      integrated$DirectionPattern != "Opposite",
    ,
    drop = FALSE
  ],
  file.path(tables_dir, "concordant_shared_nominal_pathways.csv"),
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

wrap_label <- function(x, width = 44) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

select_overview_pathways <- function(gsea_table, extra_pathway_ids = character()) {
  candidates <- gsea_table[
    gsea_table$pvalue < nominal_threshold, , drop = FALSE
  ]
  evidence_label <- "exploratory nominal P < 0.05"
  if (nrow(candidates) == 0) {
    candidates <- gsea_table
    evidence_label <- "top-ranked pathways"
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
  selected <- unique(rbind(negative, positive))
  if (length(extra_pathway_ids) > 0) {
    extra <- gsea_table[
      gsea_table$ID %in% extra_pathway_ids,
      ,
      drop = FALSE
    ]
    selected <- unique(rbind(selected, extra))
  }
  list(data = selected, label = "P value < 0.05")
}

make_overview_plot <- function(
  gsea_table,
  omics_label,
  extra_pathway_ids = character(),
  extra_pathway_suffix = ""
) {
  selected <- select_overview_pathways(gsea_table, extra_pathway_ids)
  plot_data <- selected$data
  plot_data$PlotLabel <- paste0(
    "[", plot_data$Database, "] ", wrap_label(plot_data$Description)
  )
  plot_data$PlotLabel[
    plot_data$ID %in% extra_pathway_ids
  ] <- paste0(
    plot_data$PlotLabel[plot_data$ID %in% extra_pathway_ids],
    extra_pathway_suffix
  )
  plot_data$PlotLabel <- factor(
    plot_data$PlotLabel,
    levels = plot_data$PlotLabel[order(plot_data$NES)]
  )
  plot_data$PScore <- -log10(pmax(plot_data$pvalue, 1e-10))
  ggplot(plot_data, aes(x = NES, y = PlotLabel)) +
    geom_vline(xintercept = 0, color = "#9AA2AE", linewidth = 0.45) +
    geom_segment(
      aes(x = 0, xend = NES, yend = PlotLabel, color = NES),
      linewidth = 0.8, alpha = 0.55
    ) +
    geom_point(aes(color = NES, size = PScore), alpha = 0.95) +
    scale_color_gradient2(
      low = "#3B6FB6", mid = "#D5DAE1", high = "#D95555", midpoint = 0
    ) +
    scale_size_continuous(range = c(2.6, 6.0)) +
    labs(
      title = paste0(omics_label, " ranked GSEA"),
      subtitle = selected$label,
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = "NES",
      size = "-log10(P value)"
    ) +
    theme_project() +
    theme(panel.grid.major.y = element_blank())
}

rna_plot <- make_overview_plot(
  rna_gsea,
  "Transcriptomic",
  extra_pathway_ids = "GO:0016042",
  extra_pathway_suffix = " (Plcxd2)"
)
protein_plot <- make_overview_plot(protein_gsea, "Proteomic")
save_ggplot(
  rna_plot, file.path(figures_dir, "01_transcriptomic_gsea"), 9.2, 7.2
)
save_ggplot(
  protein_plot, file.path(figures_dir, "02_proteomic_gsea"), 9.2, 7.2
)

make_gsea_curve <- function(
  rank_vector,
  pathway_genes,
  pathway_result,
  highlighted_gene_id,
  highlighted_gene_symbol
) {
  pathway_genes <- intersect(as.character(pathway_genes), names(rank_vector))
  hit <- names(rank_vector) %in% pathway_genes
  if (!any(hit) || all(hit)) {
    stop("The requested pathway cannot be plotted from the ranked gene list.")
  }

  hit_weights <- abs(rank_vector)
  running_score <- cumsum(ifelse(
    hit,
    hit_weights / sum(hit_weights[hit]),
    -1 / sum(!hit)
  ))
  curve_data <- data.frame(
    Rank = seq_along(rank_vector),
    RunningES = running_score,
    Statistic = unname(rank_vector),
    Hit = hit,
    stringsAsFactors = FALSE
  )
  peak_rank <- if (pathway_result$enrichmentScore[[1]] >= 0) {
    which.max(running_score)
  } else {
    which.min(running_score)
  }
  highlighted_rank <- match(highlighted_gene_id, names(rank_vector))

  curve_plot <- ggplot(curve_data, aes(x = Rank, y = RunningES)) +
    geom_hline(yintercept = 0, color = "#9AA2AE", linewidth = 0.4) +
    geom_area(
      aes(y = RunningES), fill = "#5B8CC0", alpha = 0.16
    ) +
    geom_line(color = "#326B9E", linewidth = 1.05) +
    geom_vline(
      xintercept = peak_rank, linetype = "dashed",
      color = "#6B7280", linewidth = 0.55
    ) +
    labs(
      title = "Lipid catabolic process",
      subtitle = paste0(
        "[GO BP] | P value = ",
        formatC(pathway_result$pvalue[[1]], format = "e", digits = 2),
        " | NES = ", formatC(pathway_result$NES[[1]], digits = 2, format = "f")
      ),
      x = NULL,
      y = "Running enrichment score"
    ) +
    theme_project(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 16, 0, 12)
    )

  hit_data <- curve_data[curve_data$Hit, , drop = FALSE]
  hit_plot <- ggplot(hit_data, aes(x = Rank)) +
    geom_linerange(
      aes(ymin = 0, ymax = 1), color = "#364152", linewidth = 0.28,
      alpha = 0.72
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL) +
    theme_project(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 16, 0, 12)
    )
  if (!is.na(highlighted_rank)) {
    hit_plot <- hit_plot +
      geom_linerange(
        data = data.frame(Rank = highlighted_rank),
        aes(x = Rank, ymin = 0, ymax = 1),
        inherit.aes = FALSE,
        color = "#D95555", linewidth = 1.35
      ) +
      annotate(
        "text", x = highlighted_rank + 0.012 * length(rank_vector), y = 0.78,
        label = highlighted_gene_symbol, color = "#D95555", fontface = "bold",
        hjust = 0, size = 3.8
      )
  }

  rank_plot <- ggplot(curve_data, aes(x = Rank)) +
    geom_area(aes(y = pmax(Statistic, 0)), fill = "#E45A5A", alpha = 0.82) +
    geom_area(aes(y = pmin(Statistic, 0)), fill = "#4F7FAF", alpha = 0.82) +
    geom_hline(yintercept = 0, color = "#9AA2AE", linewidth = 0.35) +
    labs(x = "Rank in ordered transcriptomic gene list", y = "Wald statistic") +
    theme_project(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.margin = margin(0, 16, 12, 12)
    )

  combined <- curve_plot / hit_plot / rank_plot +
    plot_layout(heights = c(3.2, 0.65, 1.3))
  list(plot = combined, curve_data = curve_data, highlighted_rank = highlighted_rank)
}

lipid_pathway_id <- "GO:0016042"
lipid_result <- rna_gsea[rna_gsea$ID == lipid_pathway_id, , drop = FALSE]
lipid_genes <- go_annotations$Gene[go_annotations$Term == lipid_pathway_id]
if (nrow(lipid_result) == 1) {
  lipid_gsea <- make_gsea_curve(
    rank_vector = rna_rank,
    pathway_genes = lipid_genes,
    pathway_result = lipid_result,
    highlighted_gene_id = "433022",
    highlighted_gene_symbol = "Plcxd2"
  )
  save_ggplot(
    lipid_gsea$plot,
    file.path(figures_dir, "05_lipid_catabolic_process_gsea"),
    9.2,
    6.4
  )
  write.csv(
    lipid_gsea$curve_data,
    file.path(tables_dir, "lipid_catabolic_process_gsea_curve.csv"),
    row.names = FALSE,
    quote = TRUE
  )
}

integrated_candidates <- integrated[
  integrated$EvidencePattern == "Both nominal P" &
    integrated$DirectionPattern != "Opposite",
  ,
  drop = FALSE
]
if (nrow(integrated_candidates) == 0) {
  integrated_candidates <- integrated
}
integrated_selected <- integrated_candidates[
  order(integrated_candidates$JointNominalP),
  ,
  drop = FALSE
]
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
    Pvalue = integrated_selected$RNA_pvalue,
    Padj = integrated_selected$RNA_padj,
    stringsAsFactors = FALSE
  ),
  data.frame(
    PathwayLabel = integrated_selected$PathwayLabel,
    Omics = "Proteomics",
    NES = integrated_selected$Protein_NES,
    Pvalue = integrated_selected$Protein_pvalue,
    Padj = integrated_selected$Protein_padj,
    stringsAsFactors = FALSE
  )
)
integrated_long$Significance <- ifelse(
  integrated_long$Pvalue < nominal_threshold,
  "Nominal P < 0.05",
  "Not significant"
)
integrated_long$Omics <- factor(
  integrated_long$Omics,
  levels = c("Transcriptomics", "Proteomics")
)
integrated_long$PathwayLabel <- factor(
  integrated_long$PathwayLabel,
  levels = levels(integrated_selected$PathwayLabel)
)

integrated_dotplot <- ggplot() +
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
    title = "Cross-omics pathway enrichment",
    subtitle = "concordant pathways | P value < 0.05",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    color = "Omics"
  ) +
  theme_project() +
  theme(panel.grid.major.y = element_blank())
save_ggplot(
  integrated_dotplot,
  file.path(figures_dir, "03_cross_omics_pathway_enrichment"),
  10.0,
  max(8.0, 2.8 + 0.34 * nrow(integrated_selected))
)

integrated_heatmap <- ggplot(
  integrated_long,
  aes(x = Omics, y = PathwayLabel, fill = NES)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "#F4F5F7", high = "#D95555", midpoint = 0
  ) +
  labs(
    title = "Integrated pathway activity",
    subtitle = "concordant pathways | P value < 0.05",
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  theme_project() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right"
  )
save_ggplot(
  integrated_heatmap,
  file.path(figures_dir, "04_integrated_pathway_heatmap"),
  8.6,
  max(8.0, 2.8 + 0.34 * nrow(integrated_selected))
)

summary_table <- data.frame(
  Metric = c(
    "Ranked RNA genes",
    "Ranked protein genes",
    "RNA pathways tested",
    "Protein pathways tested",
    "RNA FDR-significant pathways",
    "Protein FDR-significant pathways",
    "RNA nominal pathways",
    "Protein nominal pathways",
    "Pathways tested in both omics",
    "Nominal pathways in both omics",
    "RNA-only nominal pathways",
    "Protein-only nominal pathways",
    "Concordant shared nominal pathways",
    "Opposite-direction shared nominal pathways"
  ),
  Value = c(
    length(rna_rank),
    length(protein_rank),
    nrow(rna_gsea),
    nrow(protein_gsea),
    sum(rna_gsea$p.adjust < fdr_threshold),
    sum(protein_gsea$p.adjust < fdr_threshold),
    sum(rna_gsea$pvalue < nominal_threshold),
    sum(protein_gsea$pvalue < nominal_threshold),
    nrow(integrated),
    sum(integrated$EvidencePattern == "Both nominal P"),
    sum(integrated$EvidencePattern == "RNA nominal P only"),
    sum(integrated$EvidencePattern == "Protein nominal P only"),
    sum(
      integrated$EvidencePattern == "Both nominal P" &
        integrated$DirectionPattern != "Opposite"
    ),
    sum(
      integrated$EvidencePattern == "Both nominal P" &
        integrated$DirectionPattern == "Opposite"
    )
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
  "Pathway analysis complete. Shared nominal-P pathways: ",
  sum(integrated$EvidencePattern == "Both nominal P"),
  "; concordant: ",
  sum(
    integrated$EvidencePattern == "Both nominal P" &
      integrated$DirectionPattern != "Opposite"
  ),
  "; opposite: ",
  sum(
    integrated$EvidencePattern == "Both nominal P" &
      integrated$DirectionPattern == "Opposite"
  ),
  "."
)
message("Results written to: ", normalizePath(output_dir))
