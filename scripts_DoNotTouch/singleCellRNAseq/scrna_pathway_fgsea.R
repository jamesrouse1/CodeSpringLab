#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: scrna_pathway_fgsea.R <out_dir> <params.tsv>")
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1]])) y else x
out_dir <- normalizePath(args[[1]], mustWork = TRUE)
params <- utils::read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
values <- setNames(as.character(params$value), as.character(params$key))
get <- function(key, default = "") {
  if (!key %in% names(values)) return(default)
  value <- trimws(values[[key]] %||% default)
  if (!nzchar(value)) default else value
}
if (!requireNamespace("fgsea", quietly = TRUE)) stop("Ranked pathway analysis requires the fgsea R package in the CodeSpringLab R runtime.")
tables <- file.path(out_dir, "tables")
figures <- file.path(out_dir, "figures")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
library_name <- get("pathway_library", "MSigDB_Hallmark_2020")
library_slug <- gsub("[^A-Za-z0-9_.-]+", "_", library_name)
library_slug <- gsub("^_+|_+$", "", library_slug)
if (!nzchar(library_slug)) library_slug <- "pathways"
gmt <- path.expand(get("pathway_gmt_file", ""))
library_url <- ""
if (!nzchar(gmt)) {
  cache_dir <- file.path(out_dir, "pathway_libraries")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gmt <- file.path(cache_dir, paste0(library_slug, ".gmt"))
  library_url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", utils::URLencode(library_name, reserved = TRUE))
  if (!file.exists(gmt) || file.info(gmt)$size <= 0) {
    temporary <- tempfile(pattern = paste0(library_slug, "_"), tmpdir = cache_dir, fileext = ".download")
    error <- tryCatch({
      utils::download.file(library_url, temporary, mode = "wb", quiet = FALSE)
      ""
    }, error = function(e) conditionMessage(e))
    if (nzchar(error) || !file.exists(temporary) || file.info(temporary)$size <= 0) stop("Could not download pathway database '", library_name, "' from Enrichr. ", error)
    downloaded <- fgsea::gmtPathways(temporary)
    if (!length(downloaded)) stop("The downloaded pathway database contained no usable gene sets: ", library_name)
    if (!file.rename(temporary, gmt)) stop("Could not save the downloaded pathway database cache: ", gmt)
  }
}
if (!file.exists(gmt) || file.info(gmt)$size <= 0) stop("The selected pathway database is unavailable: ", library_name)
de_path <- file.path(tables, "pseudobulk_differential_expression.tsv")
if (!file.exists(de_path)) stop("Best-practice pathway analysis requires the pseudobulk DESeq2 result. Complete pseudobulk differential expression first.")
de <- utils::read.delim(de_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene", "stat") %in% names(de))) stop("The pseudobulk DE result does not contain gene and Wald-statistic columns.")
de <- de[is.finite(de$stat) & nzchar(as.character(de$gene)), , drop = FALSE]
de <- de[order(abs(de$stat), decreasing = TRUE), , drop = FALSE]
de <- de[!duplicated(de$gene), , drop = FALSE]
ranks <- setNames(as.numeric(de$stat), as.character(de$gene))
species <- tolower(get("pathway_species", "human"))
mapping_note <- "Input ranks already use human gene symbols."
mapped_genes <- length(ranks)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else ""
ortholog_path <- trimws(get("pathway_ortholog_file", ""))
if (!nzchar(ortholog_path)) ortholog_path <- normalizePath(file.path(dirname(script_path), "..", "reference", "mouse_human_orthologs_MGI.tsv"), mustWork = FALSE)
read_orthologs <- function(path) {
  if (!file.exists(ortholog_path)) stop("Mouse pathway analysis requires the bundled MGI mouse-human ortholog table: ", ortholog_path)
  reader <- if (grepl("\\.csv$", path, ignore.case = TRUE)) utils::read.csv else utils::read.delim
  orthologs <- reader(path, check.names = FALSE, stringsAsFactors = FALSE)
  normalized <- setNames(names(orthologs), gsub("[^a-z0-9]+", "", tolower(names(orthologs))))
  pick <- function(candidates) {
    hits <- unname(normalized[gsub("[^a-z0-9]+", "", tolower(candidates))])
    hits <- hits[!is.na(hits)]
    if (length(hits)) hits[[1]] else ""
  }
  mouse_col <- pick(c("mouse_gene_symbol", "mouse_symbol", "mgi_symbol", "marker_symbol", "external_gene_name", "mouse_gene_name"))
  human_col <- pick(c("human_gene_symbol", "human_symbol", "hgnc_symbol", "human_gene_name", "hsapiens_homolog_associated_gene_name"))
  if (!nzchar(mouse_col) || !nzchar(human_col)) stop("Ortholog table needs recognizable mouse and human gene-symbol columns.")
  names(orthologs)[names(orthologs) == mouse_col] <- "mouse_gene_symbol"
  names(orthologs)[names(orthologs) == human_col] <- "human_gene_symbol"
  orthologs
}
orthologs <- NULL
if (species %in% c("auto", "mouse")) orthologs <- read_orthologs(ortholog_path)
if (identical(species, "auto")) {
  mouse_overlap <- sum(names(ranks) %in% orthologs$mouse_gene_symbol)
  human_overlap <- sum(names(ranks) %in% orthologs$human_gene_symbol)
  if (mouse_overlap == human_overlap) stop("Could not infer mouse versus human differential-expression symbols. Choose the species explicitly in CodeSpringApp.")
  species <- if (mouse_overlap > human_overlap) "mouse" else "human"
}
if (identical(species, "mouse")) {
  orthologs <- orthologs[nzchar(orthologs$mouse_gene_symbol) & nzchar(orthologs$human_gene_symbol), , drop = FALSE]
  required_ortholog <- c("mouse_gene_symbol", "human_gene_symbol")
  mouse_counts <- table(orthologs$mouse_gene_symbol)
  orthologs <- orthologs[mouse_counts[orthologs$mouse_gene_symbol] == 1L, , drop = FALSE]
  rank_table <- data.frame(mouse_gene_symbol = names(ranks), stat = as.numeric(ranks), stringsAsFactors = FALSE)
  mapped <- merge(rank_table, orthologs[, required_ortholog], by = "mouse_gene_symbol", all = FALSE, sort = FALSE)
  mapped <- mapped[is.finite(mapped$stat) & nzchar(mapped$human_gene_symbol), , drop = FALSE]
  mapped <- mapped[order(abs(mapped$stat), decreasing = TRUE), , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$human_gene_symbol), , drop = FALSE]
  if (NROW(mapped) < 100L) stop("Too few mouse genes mapped to human ortholog symbols for reliable pathway analysis (", NROW(mapped), ").")
  utils::write.table(mapped, file.path(tables, paste0("pathway_mouse_to_human_mapping__", library_slug, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  ranks <- setNames(mapped$stat, mapped$human_gene_symbol)
  mapped_genes <- NROW(mapped)
  mapping_note <- paste0("Mouse ranks mapped to human ortholog symbols with bundled MGI table; ambiguous mouse-to-many mappings excluded; duplicate human targets retained by largest absolute Wald statistic. Mapped genes: ", mapped_genes, ".")
}
ranks <- sort(ranks, decreasing = TRUE)
pathways <- fgsea::gmtPathways(gmt)
if (!length(pathways)) stop("The selected pathway database contained no usable gene sets: ", library_name)
result <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15L, maxSize = 500L, eps = 0)
result <- as.data.frame(result)
if ("leadingEdge" %in% names(result)) result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
result <- result[order(result$padj, -abs(result$NES), na.last = TRUE), , drop = FALSE]
specific_result <- file.path(tables, paste0("pathway_fgsea_ranked__", library_slug, ".tsv"))
utils::write.table(result, specific_result, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(result, file.path(tables, "pathway_fgsea_ranked.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
provenance <- data.frame(
  pathway_library = library_name,
  source_url = library_url,
  cached_gmt = normalizePath(gmt, winslash = "/", mustWork = TRUE),
  pathway_count = length(pathways),
  species = species,
  ranked_genes = length(ranks),
  mapping_note = mapping_note,
  generated = as.character(Sys.time()),
  stringsAsFactors = FALSE
)
utils::write.table(provenance, file.path(tables, paste0("pathway_source__", library_slug, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE) && NROW(result)) {
  shown <- head(result[is.finite(result$padj) & is.finite(result$NES), , drop = FALSE], 20L)
  if (NROW(shown)) {
    shown$pathway <- factor(shown$pathway, levels = rev(shown$pathway))
    jpplot_colors <- c("#90C3DD", "#C2E4EF", "#ECF7E1", "#FEF4AF", "#FDD484", "#FBA25B", "#F0653F", "#D42D26", "#A50026")
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = NES, y = pathway, color = -log10(pmax(padj, 1e-300)))) + ggplot2::geom_point(size = 3) + ggplot2::scale_color_gradientn(colors = jpplot_colors, name = "−log10 FDR") + ggplot2::labs(x = "Normalized enrichment score", y = NULL, title = paste("Ranked pseudobulk pathway analysis —", library_name)) + ggplot2::theme_classic(base_size = 11) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(face = "bold"))
    ggplot2::ggsave(file.path(figures, paste0("pathway_fgsea_top20__", library_slug, ".png")), plot, width = 9, height = 6, dpi = 180)
    ggplot2::ggsave(file.path(figures, "pathway_fgsea_top20.png"), plot, width = 9, height = 6, dpi = 180)
  }
}
