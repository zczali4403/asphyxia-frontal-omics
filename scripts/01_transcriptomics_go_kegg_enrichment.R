#!/usr/bin/env Rscript

# Standalone GO and KEGG ORA for the Step 1 transcriptomic DEG set.
# Run from the repository root. This does not execute differential analysis.
required <- c("AnnotationDbi", "org.Mm.eg.db", "GO.db", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                          call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))
mouse_db <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
go_db <- getExportedValue("GO.db", "GO.db")

args <- commandArgs(trailingOnly = TRUE)
transcriptomics_dir <- if (length(args) >= 1) args[[1]] else
  "results/01_transcriptomics"
reference_dir <- if (length(args) >= 2) args[[2]] else "data/reference"
output_dir <- if (length(args) >= 3) args[[3]] else
  file.path(transcriptomics_dir, "enrichment")
alpha <- 0.05
min_size <- 10L
max_size <- 500L
go_top_n <- 8L
kegg_top_n <- 15L

read_input <- function(path, required_columns, tab = FALSE) {
  if (!file.exists(path)) stop("Input not found: ", path, call. = FALSE)
  x <- if (tab) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE) else
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  absent <- setdiff(required_columns, names(x))
  if (length(absent)) stop("Missing columns in ", path, ": ",
                           paste(absent, collapse = ", "), call. = FALSE)
  x
}
all_file <- file.path(transcriptomics_dir, "tables", "deseq2_all_genes.csv")
deg_file <- file.path(transcriptomics_dir, "tables",
                      "deseq2_effect_size_filtered_genes.csv")
all_genes <- read_input(all_file, c("GeneID", "GeneName", "pvalue"))
degs <- read_input(deg_file, c("GeneID", "log2FoldChange", "padj"))
if (!nrow(degs)) stop("No DEGs in the effect-size-filtered table.")
if (any(!is.finite(degs$padj) | !is.finite(degs$log2FoldChange) |
        degs$padj >= alpha | abs(degs$log2FoldChange) < 0.5)) {
  stop("DEG input must meet FDR < 0.05 and |shrunken log2FC| >= 0.5.")
}
tested <- all_genes[is.finite(all_genes$pvalue), , drop = FALSE]
if (any(!degs$GeneID %in% tested$GeneID)) {
  stop("DEG table contains genes absent from the tested-gene background.")
}
ensembl <- sub("\\..*$", "", as.character(tested$GeneID))
mapped <- suppressMessages(AnnotationDbi::mapIds(
  mouse_db, keys = unique(ensembl), column = "ENTREZID",
  keytype = "ENSEMBL", multiVals = "first"
))
tested$EntrezID <- unname(mapped[ensembl])
tested$DEG <- tested$GeneID %in% degs$GeneID
valid_ids <- function(x) unique(as.character(x[!is.na(x) & nzchar(x)]))
universe <- valid_ids(tested$EntrezID)
candidates <- valid_ids(tested$EntrezID[tested$DEG])
if (!length(candidates)) stop("No DEGs mapped to Entrez IDs.")

message("Loading local GO and KEGG annotations")
go <- suppressMessages(AnnotationDbi::select(
  mouse_db, keys = universe, keytype = "ENTREZID",
  columns = c("GOALL", "ONTOLOGYALL")
))
go <- unique(go[!is.na(go$GOALL) & go$ONTOLOGYALL %in% c("BP", "MF", "CC"),
                c("GOALL", "ENTREZID", "ONTOLOGYALL"), drop = FALSE])
names(go) <- c("TermID", "EntrezID", "Category")
go_names <- suppressMessages(AnnotationDbi::select(
  go_db, keys = unique(go$TermID), keytype = "GOID", columns = "TERM"
))
go_names <- go_names[!duplicated(go_names$GOID), , drop = FALSE]
go$Description <- go_names$TERM[match(go$TermID, go_names$GOID)]

kegg <- read_input(file.path(reference_dir, "kegg_mouse_gene_pathways.tsv"),
                   c("Term", "Gene"), tab = TRUE)
kegg_names <- read_input(file.path(reference_dir, "kegg_mouse_pathway_names.tsv"),
                         c("Term", "Name"), tab = TRUE)
kegg <- unique(data.frame(
  TermID = as.character(kegg$Term), EntrezID = as.character(kegg$Gene),
  Category = "KEGG",
  Description = kegg_names$Name[match(kegg$Term, kegg_names$Term)],
  stringsAsFactors = FALSE
))

# Each ontology/database uses its own annotated measured-gene background.
# Include zero-overlap eligible sets in the BH correction family.
ora <- function(annotation) {
  annotation <- annotation[!is.na(annotation$Description) &
                             !is.na(annotation$EntrezID), , drop = FALSE]
  rows <- list()
  counts <- list()
  for (category in unique(annotation$Category)) {
    ann <- annotation[annotation$Category == category, , drop = FALSE]
    background <- intersect(universe, ann$EntrezID)
    query <- intersect(candidates, background)
    counts[[category]] <- data.frame(Category = category,
      BackgroundGenes = length(background), DEGGenes = length(query))
    if (!length(query)) next
    sets <- split(ann, ann$TermID)
    category_rows <- lapply(sets, function(term) {
      members <- intersect(unique(term$EntrezID), background)
      m <- length(members)
      if (m < min_size || m > max_size) return(NULL)
      hits <- intersect(query, members)
      k <- length(hits)
      data.frame(Category = category, TermID = term$TermID[1],
        Description = term$Description[1], Count = k, SetSize = m,
        DEGSize = length(query), BackgroundSize = length(background),
        GeneRatio = k / length(query),
        GeneRatioFraction = paste0(k, "/", length(query)),
        BgRatio = paste0(m, "/", length(background)),
        FoldEnrichment = (k / length(query)) / (m / length(background)),
        Pvalue = phyper(k - 1, m, length(background) - m, length(query),
                        lower.tail = FALSE),
        OverlapEntrez = paste(hits, collapse = "/"),
        stringsAsFactors = FALSE)
    })
    category_rows <- Filter(Negate(is.null), category_rows)
    if (!length(category_rows)) next
    result <- do.call(rbind, category_rows)
    result$FDR <- p.adjust(result$Pvalue, method = "BH")
    rows[[category]] <- result
  }
  if (!length(rows)) stop("No eligible gene sets for this database.")
  result <- do.call(rbind, rows)
  hit_ids <- valid_ids(unlist(strsplit(result$OverlapEntrez, "/", fixed = TRUE)))
  symbols <- if (length(hit_ids)) suppressMessages(AnnotationDbi::mapIds(
    mouse_db, keys = hit_ids, column = "SYMBOL", keytype = "ENTREZID",
    multiVals = "first"
  )) else setNames(character(), character())
  result$OverlapGenes <- vapply(result$OverlapEntrez, function(ids) {
    if (!nzchar(ids)) return("")
    keys <- strsplit(ids, "/", fixed = TRUE)[[1]]
    names <- unname(symbols[keys])
    paste(ifelse(is.na(names), keys, names), collapse = "/")
  }, character(1))
  result <- result[order(result$FDR, result$Pvalue, result$TermID), , drop = FALSE]
  rownames(result) <- NULL
  list(results = result, counts = do.call(rbind, counts))
}
message("Testing GO BP/MF/CC and KEGG enrichment")
go_analysis <- ora(go)
kegg_analysis <- ora(kegg)
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tested, file.path(tables_dir, "gene_to_entrez_mapping.csv"),
          row.names = FALSE)
write.csv(rbind(go_analysis$counts, kegg_analysis$counts),
          file.path(tables_dir, "annotation_background_counts.csv"),
          row.names = FALSE)

theme_enrichment <- function() {
  theme_minimal(base_size = 12) + theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5,
                              color = "#17243A"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10, color = "#344155"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#EEF1F6", color = NA),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(15, 25, 15, 15))
}
save_plot <- function(p, filename, height) {
  ggsave(file.path(figures_dir, paste0(filename, ".png")), p,
         width = 12, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")), p,
         width = 12, height = height, bg = "white")
}
label_terms <- function(x, width) {
  x$PlotID <- factor(x$TermID, levels = rev(x$TermID))
  x$Label <- vapply(x$Description, function(z)
    paste(strwrap(z, width = width), collapse = "\n"), character(1))
  x$Score <- -log10(pmax(x$FDR, 1e-300))
  x
}
for (db in c("go", "kegg")) {
  result <- if (db == "go") go_analysis$results else kegg_analysis$results
  significant <- result[result$FDR < alpha & result$Count > 0, , drop = FALSE]
  write.csv(result, file.path(tables_dir, paste0(db, "_ora_all_terms.csv")),
            row.names = FALSE)
  write.csv(significant, file.path(tables_dir, paste0(db, "_ora_fdr_significant.csv")),
            row.names = FALSE)
  if (db == "go") {
    selected <- lapply(c("BP", "MF", "CC"), function(category)
      head(significant[significant$Category == category, , drop = FALSE], go_top_n))
    plot_terms <- do.call(rbind, selected)
  } else plot_terms <- head(significant, kegg_top_n)
  write.csv(plot_terms, file.path(tables_dir, paste0(db, "_plot_terms.csv")),
            row.names = FALSE)
  filename <- if (db == "go") "01_go_enrichment_dotplot" else
    "02_kegg_enrichment_barplot"
  if (!nrow(plot_terms)) {
    p <- ggplot() + annotate("text", x = 0, y = 0,
      label = "No enriched terms at FDR < 0.05", size = 5, color = "#526174") +
      labs(title = paste(toupper(db), "enrichment")) + theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 19))
    save_plot(p, filename, 4)
    next
  }
  x <- label_terms(plot_terms, 48)
  label_map <- setNames(x$Label, x$TermID)
  if (db == "go") {
    x$Category <- factor(x$Category, levels = c("BP", "MF", "CC"))
    p <- ggplot(x, aes(GeneRatio, PlotID)) +
      geom_point(aes(size = Count, color = Score), alpha = 0.95) +
      scale_color_gradient(low = "#C7B5E0", high = "#62368D",
                            name = expression(-log[10](FDR))) +
      scale_size_continuous(range = c(3.5, 8), name = "DEGs") +
      scale_y_discrete(labels = label_map) +
      scale_x_continuous(labels = function(z) paste0(round(z * 100), "%")) +
      facet_grid(Category ~ ., scales = "free_y", space = "free_y",
        labeller = as_labeller(c(BP = "Biological Process",
          MF = "Molecular Function", CC = "Cellular Component"))) +
      labs(title = "GO enrichment of differentially expressed genes",
            x = "Gene ratio", y = NULL) + theme_enrichment()
    save_plot(p, filename, max(5, 2 + nrow(x) * 0.43))
  } else {
    p <- ggplot(x, aes(Score, PlotID, fill = Score)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = paste0(Count, " genes")), hjust = -0.15,
                 size = 3.4, color = "#344155") +
      scale_fill_gradient(low = "#A7D9C8", high = "#247A66",
                           name = expression(-log[10](FDR))) +
      scale_y_discrete(labels = label_map) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title = "KEGG enrichment of differentially expressed genes",
            x = expression(-log[10](FDR)), y = NULL) + theme_enrichment()
    save_plot(p, filename, max(4.5, 2 + nrow(x) * 0.48))
  }
}
summary <- data.frame(
  Metric = c("DEG input features", "Mapped tested Entrez genes", "Mapped DEG Entrez genes",
    "GO BP FDR enriched terms", "GO MF FDR enriched terms", "GO CC FDR enriched terms",
    "KEGG FDR enriched terms"),
  Value = c(nrow(degs), length(universe), length(candidates),
    vapply(c("BP", "MF", "CC"), function(category)
      sum(go_analysis$results$Category == category & go_analysis$results$FDR < alpha),
      integer(1)), sum(kegg_analysis$results$FDR < alpha)))
write.table(summary, file.path(output_dir, "analysis_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
parameters <- data.frame(
  Parameter = c("DEG input", "DEG definition", "Background", "ID mapping", "GO annotation",
    "KEGG annotation", "Test", "Gene-set size", "Correction", "Plot selection"),
  Value = c(deg_file, "FDR < 0.05 and |shrunken log2FC| >= 0.5",
    "Finite DESeq2 P value; measured Entrez genes annotated within each ontology/database",
    "Ensembl to Entrez, first mapping, unique Entrez IDs",
    "org.Mm.eg.db GOALL; BP, MF, CC", reference_dir,
    "One-sided hypergeometric over-representation test", "10-500 measured genes",
    "BH within BP, MF, CC and KEGG; includes eligible zero-overlap terms",
    "FDR < 0.05 only; top 8 per GO ontology and top 15 KEGG pathways"))
write.table(parameters, file.path(output_dir, "analysis_parameters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("GO/KEGG enrichment complete. Outputs: ", output_dir)
