#!/usr/bin/env Rscript

# Over-representation analysis of frontal-cortex DEGs and exploratory DEPs.
# GO annotations are provided by org.Mm.eg.db; KEGG uses the local mouse cache.

required_packages <- c("AnnotationDbi", "org.Mm.eg.db", "GO.db",
                       "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
})

args <- commandArgs(trailingOnly = TRUE)
rna_all_file <- if (length(args) >= 1) args[[1]] else
  "results/01_transcriptomics/tables/deseq2_all_genes.csv"
rna_deg_file <- if (length(args) >= 2) args[[2]] else
  "results/01_transcriptomics/tables/deseq2_effect_size_filtered_genes.csv"
protein_all_file <- if (length(args) >= 3) args[[3]] else
  "results/02_proteomics/tables/limma_all_proteins.csv"
protein_dep_file <- if (length(args) >= 4) args[[4]] else
  "results/02_proteomics/tables/limma_exploratory_nominal_p_proteins.csv"
reference_dir <- if (length(args) >= 5) args[[5]] else "data/reference"
output_dir <- if (length(args) >= 6) args[[6]] else
  "results/07_differential_feature_enrichment"

min_set_size <- 10L
max_set_size <- 500L
alpha <- 0.05
max_plot_terms <- 14L
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
obsolete_protein_plots <- file.path(
  figures_dir,
  as.vector(outer(c("02_dep_enrichment_dotplot", "02_dep_enrichment_circle"),
                  c(".png", ".pdf"), paste0))
)
unlink(obsolete_protein_plots[file.exists(obsolete_protein_plots)])

read_results <- function(path, required) {
  if (!file.exists(path)) stop("Input not found: ", path, call. = FALSE)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"))
  if (!all(required %in% names(x))) {
    stop("Missing columns in ", path, ": ",
         paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  x
}

message("Reading differential-expression results")
rna_all <- read_results(rna_all_file, c("GeneID", "GeneName", "pvalue"))
rna_deg <- read_results(rna_deg_file,
                        c("GeneID", "GeneName", "log2FoldChange", "padj"))
protein_all <- read_results(protein_all_file,
                            c("Accession", "GeneID", "GeneName",
                              "log2FoldChange", "pvalue"))
protein_dep <- read_results(protein_dep_file,
                            c("Accession", "GeneID", "GeneName",
                              "log2FoldChange", "pvalue"))

if (any(rna_deg$padj >= alpha | abs(rna_deg$log2FoldChange) < 0.5,
        na.rm = TRUE)) {
  stop("RNA candidate table does not match FDR < 0.05 and |log2FC| >= 0.5.")
}
if (any(protein_dep$pvalue >= alpha |
        abs(protein_dep$log2FoldChange) < 0.5, na.rm = TRUE)) {
  stop("Protein candidate table does not match P < 0.05 and |log2FC| >= 0.5.")
}

# Exclude genes without a usable differential-test P value from each universe.
rna_tested <- rna_all[is.finite(rna_all$pvalue), , drop = FALSE]
protein_tested <- protein_all[is.finite(protein_all$pvalue), , drop = FALSE]
rna_ensembl <- sub("\\..*$", "", as.character(rna_tested$GeneID))
rna_entrez <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = unique(rna_ensembl), column = "ENTREZID",
  keytype = "ENSEMBL", multiVals = "first"
))
rna_tested$EntrezID <- unname(rna_entrez[rna_ensembl])
protein_tested$EntrezID <- as.character(protein_tested$GeneID)
rna_candidates <- rna_tested[
  rna_tested$GeneID %in% rna_deg$GeneID, , drop = FALSE
]
protein_candidates <- protein_tested[
  protein_tested$Accession %in% protein_dep$Accession, , drop = FALSE
]
if (!nrow(protein_candidates)) {
  protein_candidates <- protein_tested[
    protein_tested$GeneID %in% protein_dep$GeneID, , drop = FALSE
  ]
}

valid_ids <- function(x) unique(as.character(x[!is.na(x) & x != ""]))
rna_universe <- valid_ids(rna_tested$EntrezID)
protein_universe <- valid_ids(protein_tested$EntrezID)
rna_gene_ids <- intersect(valid_ids(rna_candidates$EntrezID), rna_universe)
protein_gene_ids <- intersect(valid_ids(protein_candidates$EntrezID),
                              protein_universe)
if (!length(rna_gene_ids) || !length(protein_gene_ids)) {
  stop("No candidate genes mapped to Entrez IDs.", call. = FALSE)
}

write.csv(
  data.frame(Omics = c(rep("Transcriptomics", nrow(rna_tested)),
                       rep("Proteomics", nrow(protein_tested))),
             FeatureID = c(rna_tested$GeneID, protein_tested$Accession),
             GeneName = c(rna_tested$GeneName, protein_tested$GeneName),
             EntrezID = c(rna_tested$EntrezID, protein_tested$EntrezID),
             Candidate = c(rna_tested$GeneID %in% rna_deg$GeneID,
                           protein_tested$Accession %in% protein_dep$Accession)),
  file.path(tables_dir, "feature_to_entrez_mapping.csv"), row.names = FALSE
)

message("Building GO BP/MF/CC and KEGG annotation sets")
all_ids <- union(rna_universe, protein_universe)
go <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db, keys = all_ids, columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "ENTREZID"
))
go <- unique(go[!is.na(go$GOALL) &
                  go$ONTOLOGYALL %in% c("BP", "MF", "CC"),
                c("GOALL", "ENTREZID", "ONTOLOGYALL"), drop = FALSE])
names(go) <- c("TermID", "EntrezID", "Ontology")
go_names <- suppressMessages(AnnotationDbi::select(
  GO.db, keys = unique(go$TermID), columns = "TERM", keytype = "GOID"
))
go_names <- go_names[!duplicated(go_names$GOID), , drop = FALSE]
go$Description <- go_names$TERM[match(go$TermID, go_names$GOID)]
go$Category <- paste("GO", go$Ontology)
go <- go[, c("Category", "TermID", "Description", "EntrezID")]

kegg_gene_path <- file.path(reference_dir, "kegg_mouse_gene_pathways.tsv")
kegg_name_path <- file.path(reference_dir, "kegg_mouse_pathway_names.tsv")
if (!file.exists(kegg_gene_path) || !file.exists(kegg_name_path)) {
  stop("Local KEGG mouse reference files are required in ", reference_dir,
       ".", call. = FALSE)
}
kegg <- read.delim(kegg_gene_path, check.names = FALSE)
kegg_names <- read.delim(kegg_name_path, check.names = FALSE)
if (!all(c("Term", "Gene") %in% names(kegg)) ||
    !all(c("Term", "Name") %in% names(kegg_names))) {
  stop("KEGG reference files lack expected columns.", call. = FALSE)
}
kegg <- unique(data.frame(
  Category = "KEGG", TermID = kegg$Term,
  Description = kegg_names$Name[match(kegg$Term, kegg_names$Term)],
  EntrezID = as.character(kegg$Gene), stringsAsFactors = FALSE
))
annotation <- unique(rbind(go, kegg))
annotation <- annotation[
  !is.na(annotation$Description) & !is.na(annotation$EntrezID),
  , drop = FALSE
]

ora <- function(candidate_ids, universe_ids, omics_label) {
  background <- intersect(universe_ids, annotation$EntrezID)
  candidates <- intersect(candidate_ids, background)
  by_term <- split(annotation, paste(annotation$Category, annotation$TermID,
                                    sep = "::"))
  rows <- lapply(by_term, function(term) {
    members <- intersect(unique(term$EntrezID), background)
    m <- length(members)
    if (m < min_set_size || m > max_set_size) return(NULL)
    hits <- intersect(members, candidates)
    k <- length(hits)
    p <- phyper(k - 1, m, length(background) - m, length(candidates),
                lower.tail = FALSE)
    data.frame(
      Omics = omics_label, Category = term$Category[1],
      TermID = term$TermID[1], Description = term$Description[1],
      Overlap = k, SetSize = m, CandidateSize = length(candidates),
      BackgroundSize = length(background),
      GeneRatio = paste0(k, "/", length(candidates)),
      BgRatio = paste0(m, "/", length(background)),
      FoldEnrichment = if (k) (k / length(candidates)) /
        (m / length(background)) else 0,
      Pvalue = p, OverlapEntrez = paste(hits, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) stop("No GO/KEGG gene sets passed size filtering.")
  result <- do.call(rbind, rows)
  # Multiple testing is corrected independently within each annotation class.
  result$padj <- ave(result$Pvalue, result$Category,
                     FUN = function(p) p.adjust(p, method = "BH"))
  symbols <- suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unique(unlist(strsplit(result$OverlapEntrez[
      result$Overlap > 0], "/", fixed = TRUE))),
    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
  ))
  result$OverlapGenes <- vapply(result$OverlapEntrez, function(ids) {
    if (!nzchar(ids)) return("")
    keys <- strsplit(ids, "/", fixed = TRUE)[[1]]
    paste(ifelse(is.na(symbols[keys]), keys, symbols[keys]), collapse = "/")
  }, character(1))
  result <- result[order(result$Pvalue, result$Category,
                         result$Description), , drop = FALSE]
  rownames(result) <- NULL
  result
}

message("Testing DEG and exploratory DEP over-representation")
rna_ora <- ora(rna_gene_ids, rna_universe, "Transcriptomics")
protein_ora <- ora(protein_gene_ids, protein_universe, "Proteomics")
write.csv(rna_ora, file.path(tables_dir, "deg_ora_all_terms.csv"),
          row.names = FALSE)
write.csv(protein_ora, file.path(tables_dir, "dep_exploratory_ora_all_terms.csv"),
          row.names = FALSE)
write.csv(rna_ora[rna_ora$padj < alpha, , drop = FALSE],
          file.path(tables_dir, "deg_ora_fdr_significant.csv"),
          row.names = FALSE)
write.csv(protein_ora[protein_ora$padj < alpha, , drop = FALSE],
          file.path(tables_dir, "dep_exploratory_ora_fdr_significant.csv"),
          row.names = FALSE)

plot_data <- function(result) {
  x <- result[result$Overlap > 0, , drop = FALSE]
  if (!nrow(x)) return(x)
  # Show FDR-significant terms when available; otherwise use uncorrected hits.
  fdr <- x[x$padj < alpha, , drop = FALSE]
  uncorrected <- x[x$Pvalue < alpha, , drop = FALSE]
  score_column <- if (nrow(fdr)) "padj" else "Pvalue"
  x <- if (nrow(fdr)) fdr else if (nrow(uncorrected)) uncorrected else x
  x <- x[order(x[[score_column]], -x$FoldEnrichment), , drop = FALSE]
  x <- head(x, max_plot_terms)
  x$Label <- paste0("[", x$Category, "] ", x$Description)
  x$Label <- vapply(x$Label, function(z)
    paste(strwrap(z, width = 43), collapse = "\n"), character(1))
  x$Label <- factor(make.unique(x$Label), levels = rev(make.unique(x$Label)))
  x$PlotScore <- -log10(pmax(x[[score_column]], 1e-300))
  attr(x, "score_label") <- if (score_column == "padj")
    "-log10(FDR)" else "-log10(P)"
  attr(x, "selection_label") <- if (score_column == "padj")
    "enrichment FDR < 0.05" else "exploratory enrichment P < 0.05"
  x
}

theme_enrichment <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 19, face = "bold", color = "#17243A"),
      plot.subtitle = element_text(size = 11, color = "#526174"),
      axis.title = element_text(size = 12, face = "bold", color = "#263449"),
      axis.text.y = element_text(size = 10.3, color = "#344155"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(20, 30, 20, 20)
    )
}

save_figure <- function(plot, basename, width, height) {
  ggsave(file.path(figures_dir, paste0(basename, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figures_dir, paste0(basename, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

rna_plot_data <- plot_data(rna_ora)
protein_plot_data <- plot_data(protein_ora)
if (nrow(rna_plot_data)) {
  rna_plot <- ggplot(rna_plot_data, aes(x = PlotScore, y = Label)) +
    geom_segment(aes(x = 0, xend = PlotScore, yend = Label),
                 color = "#C9BDDF", linewidth = 1.0) +
    geom_point(aes(size = Overlap), shape = 21, fill = "#60418D",
               color = "white", stroke = 0.7) +
    scale_size_continuous(range = c(4, 9)) +
    labs(title = "Differentially expressed gene enrichment",
         subtitle = paste("Asphyxia vs Control | GO and KEGG |",
                          attr(rna_plot_data, "selection_label")),
         x = attr(rna_plot_data, "score_label"), y = NULL,
         size = "DEGs") +
    theme_enrichment()
  save_figure(rna_plot, "01_deg_enrichment_lollipop", 11.5, 8.1)
}
if (nrow(protein_plot_data)) {
  protein_plot_data$BarLabel <- vapply(
    protein_plot_data$Description,
    function(z) paste(strwrap(z, width = 44), collapse = "\n"),
    character(1)
  )
  protein_plot_data$BarLabel <- factor(
    make.unique(protein_plot_data$BarLabel),
    levels = rev(make.unique(protein_plot_data$BarLabel))
  )
  category_colors <- c("GO BP" = "#E47CB7", "GO MF" = "#3F80B8",
                       "GO CC" = "#9BD14D", "KEGG" = "#934DA2")
  protein_plot <- ggplot(protein_plot_data,
                         aes(x = PlotScore, y = BarLabel, fill = Category)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.4) +
    geom_text(aes(label = paste0(Overlap, " proteins")),
              hjust = -0.18, size = 3.2, color = "#344155") +
    scale_fill_manual(values = category_colors, name = "Pathway class") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Differential protein pathway enrichment",
         subtitle = paste("Asphyxia vs Control |",
                          attr(protein_plot_data, "selection_label")),
         x = attr(protein_plot_data, "score_label"), y = NULL) +
    theme_enrichment() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  save_figure(protein_plot, "02_dep_enrichment_barplot", 11.5, 8.1)
}

summary <- data.frame(
  Metric = c("DEG input features", "DEG mapped Entrez genes",
             "DEG pathway-annotated genes", "RNA tested Entrez genes",
             "RNA pathway-annotated background genes",
             "RNA enrichment terms with uncorrected P < 0.05",
             "RNA FDR enriched terms", "Exploratory DEP input features",
             "DEP mapped Entrez genes", "DEP pathway-annotated genes",
             "Protein tested Entrez genes",
             "Protein pathway-annotated background genes",
             "Protein enrichment terms with uncorrected P < 0.05",
             "Protein FDR enriched terms"),
  Value = c(nrow(rna_deg), length(rna_gene_ids),
            length(intersect(rna_gene_ids, annotation$EntrezID)),
            length(rna_universe),
            length(intersect(rna_universe, annotation$EntrezID)),
            sum(rna_ora$Pvalue < alpha), sum(rna_ora$padj < alpha),
            nrow(protein_dep), length(protein_gene_ids),
            length(intersect(protein_gene_ids, annotation$EntrezID)),
            length(protein_universe),
            length(intersect(protein_universe, annotation$EntrezID)),
            sum(protein_ora$Pvalue < alpha),
            sum(protein_ora$padj < alpha))
)
write.table(summary, file.path(output_dir, "analysis_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
parameters <- data.frame(
  Parameter = c("RNA candidate definition", "Protein candidate definition",
                "RNA background", "Protein background", "Annotation",
                "Gene-set size", "Multiple-testing correction"),
  Value = c("FDR < 0.05 and |shrunken log2FC| >= 0.5",
            "P value < 0.05 and |log2FC| >= 0.5 (exploratory)",
            "Genes with finite DESeq2 P value and Entrez mapping",
            "Proteins with finite limma P value and Entrez mapping",
            "org.Mm.eg.db GOALL; local mouse KEGG cache",
            paste(min_set_size, "to", max_set_size),
            "BH separately within GO BP, GO MF, GO CC, and KEGG")
)
write.table(parameters, file.path(output_dir, "analysis_parameters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Enrichment complete: RNA FDR terms = ", sum(rna_ora$padj < alpha),
        "; protein FDR terms = ", sum(protein_ora$padj < alpha), ".")
