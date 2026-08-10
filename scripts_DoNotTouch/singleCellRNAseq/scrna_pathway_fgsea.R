#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: scrna_pathway_fgsea.R <out_dir> <params.tsv>")
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1]])) y else x
out_dir <- normalizePath(args[[1]], mustWork = TRUE)
params <- utils::read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
values <- setNames(as.character(params$value), as.character(params$key))
get <- function(key, default = "") trimws(values[[key]] %||% default)
gmt <- path.expand(get("pathway_gmt_file"))
if (!nzchar(gmt) || !file.exists(gmt)) stop("Choose a readable GMT gene-set file for pathway analysis.")
if (!requireNamespace("fgsea", quietly = TRUE)) stop("Ranked pathway analysis requires the fgsea R package in the CodeSpringLab R runtime.")
tables <- file.path(out_dir, "tables")
de_path <- file.path(tables, "pseudobulk_differential_expression.tsv")
if (!file.exists(de_path)) stop("Best-practice pathway analysis requires the pseudobulk DESeq2 result. Complete pseudobulk differential expression first.")
de <- utils::read.delim(de_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene", "stat") %in% names(de))) stop("The pseudobulk DE result does not contain gene and Wald-statistic columns.")
de <- de[is.finite(de$stat) & nzchar(as.character(de$gene)), , drop = FALSE]
de <- de[order(abs(de$stat), decreasing = TRUE), , drop = FALSE]
de <- de[!duplicated(de$gene), , drop = FALSE]
ranks <- setNames(as.numeric(de$stat), as.character(de$gene)); ranks <- sort(ranks, decreasing = TRUE)
pathways <- fgsea::gmtPathways(gmt)
result <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15L, maxSize = 500L, eps = 0)
result <- as.data.frame(result)
if ("leadingEdge" %in% names(result)) result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
result <- result[order(result$padj, -abs(result$NES), na.last = TRUE), , drop = FALSE]
utils::write.table(result, file.path(tables, "pathway_fgsea_ranked.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE) && NROW(result)) {
  shown <- head(result[is.finite(result$padj) & is.finite(result$NES), , drop = FALSE], 20L)
  if (NROW(shown)) {
    shown$pathway <- factor(shown$pathway, levels = rev(shown$pathway))
    jpplot_colors <- c("#90C3DD", "#C2E4EF", "#ECF7E1", "#FEF4AF", "#FDD484", "#FBA25B", "#F0653F", "#D42D26", "#A50026")
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = NES, y = pathway, color = -log10(pmax(padj, 1e-300)))) + ggplot2::geom_point(size = 3) + ggplot2::scale_color_gradientn(colors = jpplot_colors, name = "−log10 FDR") + ggplot2::labs(x = "Normalized enrichment score", y = NULL, title = "Ranked pseudobulk pathway analysis") + ggplot2::theme_classic(base_size = 11) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(face = "bold"))
    ggplot2::ggsave(file.path(out_dir, "figures", "pathway_fgsea_top20.png"), plot, width = 9, height = 6, dpi = 180)
  }
}
