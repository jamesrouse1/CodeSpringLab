#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: scrna_pseudobulk_deseq2.R <out_dir> <params.tsv>")
out_dir <- normalizePath(args[[1]], mustWork = TRUE)
params <- utils::read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
values <- setNames(as.character(params$value), as.character(params$key))
get <- function(key, default = "") trimws(values[[key]] %||% default)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1]])) y else x
method <- tolower(get("de_method", "both"))
if (method %in% c("cell", "cell_level")) quit(save = "no", status = 0L)
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Pseudobulk differential expression requires the DESeq2 R package in the CodeSpringLab R runtime.")
save_pseudobulk_de_plots <- function(result, normalized, design, comparison, reference, file_slug, output_dir) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    effect <- suppressWarnings(as.numeric(result$log2FoldChange)); pvalue <- suppressWarnings(as.numeric(result$padj))
    volcano_data <- data.frame(effect = effect, significance = -log10(pmax(pvalue, 1e-300)), stringsAsFactors = FALSE)
    volcano_data$status <- ifelse(is.finite(pvalue) & pvalue < 0.05 & is.finite(effect) & effect >= 0.5, paste0("Higher in ", comparison), ifelse(is.finite(pvalue) & pvalue < 0.05 & is.finite(effect) & effect <= -0.5, paste0("Higher in ", reference), "Not significant"))
    volcano_data <- volcano_data[is.finite(volcano_data$effect) & is.finite(volcano_data$significance), , drop = FALSE]
    colors <- setNames(c("#B8C2CC", "#C43C39", "#2B6CB0"), c("Not significant", paste0("Higher in ", comparison), paste0("Higher in ", reference)))
    plot <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$effect, y = .data$significance, color = .data$status)) + ggplot2::geom_point(size = 1.15, alpha = 0.65) + ggplot2::scale_color_manual(values = colors, breaks = names(colors), name = NULL) + ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#6B7280", linewidth = 0.35) + ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6B7280", linewidth = 0.35) + ggplot2::labs(title = paste0(comparison, " vs ", reference), subtitle = "Pseudobulk DESeq2 differential expression", x = "Log2 fold change", y = "-log10 adjusted P value") + ggplot2::theme_classic(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")
    ggplot2::ggsave(file.path(output_dir, paste0("volcano__pseudobulk_DESeq2__", file_slug, ".png")), plot, width = 8, height = 6, dpi = 180)
  }
  genes <- head(as.character(result$gene[order(result$padj, -abs(result$log2FoldChange), na.last = TRUE)]), 30L)
  genes <- intersect(genes[nzchar(genes)], rownames(normalized))
  if (length(genes) && requireNamespace("pheatmap", quietly = TRUE)) {
    matrix <- as.matrix(normalized[genes, as.character(design$sample_id), drop = FALSE]); matrix <- log2(matrix + 1)
    matrix <- t(scale(t(matrix))); matrix[!is.finite(matrix)] <- 0; matrix <- pmax(-2.5, pmin(2.5, matrix))
    annotation <- data.frame(Group = factor(design$group, levels = c(reference, comparison)), row.names = as.character(design$sample_id))
    path <- file.path(output_dir, paste0("heatmap__pseudobulk_DESeq2__", file_slug, ".png"))
    grDevices::png(path, width = 1800, height = 1350, res = 180, type = "cairo")
    pheatmap::pheatmap(matrix, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = annotation, color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#C51B29"))(101), border_color = NA, main = paste0("Top differential genes: ", comparison, " vs ", reference))
    grDevices::dev.off()
  }
  invisible(TRUE)
}
tables <- file.path(out_dir, "tables")
manifest_path <- file.path(tables, "differential_jobs.tsv")
if (!file.exists(manifest_path)) stop("Differential job manifest is missing; rerun the differential-expression stage.")
manifest <- utils::read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("comparison", "reference", "annotation_value", "output_dir", "file_slug", "status", "count_file", "design_file", "cell_result_file")
if (!all(required %in% names(manifest))) stop("Differential job manifest has an invalid schema.")
covariates <- Filter(nzchar, trimws(strsplit(get("de_covariates", ""), ",", fixed = TRUE)[[1]]))
manifest$pseudobulk_status <- ifelse(manifest$status == "prepared", "pending", "not prepared")
manifest$pseudobulk_note <- manifest$note
manifest$pseudobulk_result_file <- ""; manifest$normalized_counts_file <- ""
successes <- 0L
for (i in seq_len(NROW(manifest))) {
  if (!identical(manifest$status[[i]], "prepared") || !nzchar(manifest$count_file[[i]]) || !nzchar(manifest$design_file[[i]])) next
  error <- tryCatch({
    count_table <- utils::read.delim(manifest$count_file[[i]], check.names = FALSE, stringsAsFactors = FALSE)
    design <- utils::read.delim(manifest$design_file[[i]], check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("sample_id", "group", "cells") %in% names(design)) || !"gene" %in% names(count_table)) stop("Pseudobulk input tables have an invalid schema.")
    reference <- manifest$reference[[i]]; comparison <- manifest$comparison[[i]]
    design <- design[design$group %in% c(reference, comparison), , drop = FALSE]
    replicates <- table(design$group)
    if (!all(c(reference, comparison) %in% names(replicates)) || any(replicates[c(reference, comparison)] < 2L)) stop("At least two independent samples are required in each group for this population.")
    count_columns <- intersect(as.character(design$sample_id), names(count_table))
    if (length(count_columns) != NROW(design)) stop("Pseudobulk count columns do not match the design.")
    counts <- as.matrix(count_table[, count_columns, drop = FALSE]); storage.mode(counts) <- "numeric"
    rownames(counts) <- make.unique(as.character(count_table$gene)); counts <- round(counts)
    keep <- rowSums(counts >= 10) >= 2L & rowSums(counts) >= 20L
    if (sum(keep) < 10L) stop("Too few expressed genes remain after low-count filtering.")
    counts <- counts[keep, , drop = FALSE]; design <- design[match(colnames(counts), design$sample_id), , drop = FALSE]
    design$group <- stats::relevel(factor(design$group), ref = reference); rownames(design) <- design$sample_id
    missing_covariates <- setdiff(covariates, names(design))
    if (length(missing_covariates)) stop("Missing covariate(s): ", paste(missing_covariates, collapse = ", "))
    for (covariate in covariates) design[[covariate]] <- factor(design[[covariate]])
    design_formula <- stats::reformulate(c(covariates, "group"))
    model_matrix <- stats::model.matrix(design_formula, design)
    if (qr(model_matrix)$rank < ncol(model_matrix)) stop("Model is not full rank; remove a covariate confounded with group.")
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = design, design = design_formula)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    result <- as.data.frame(DESeq2::results(dds, contrast = c("group", comparison, reference), independentFiltering = TRUE))
    result <- data.frame(gene = rownames(result), comparison = paste0(comparison, "_vs_", reference), population = manifest$annotation_value[[i]], result, analysis_level = "sample-level pseudobulk DESeq2", check.names = FALSE)
    result <- result[order(result$padj, -abs(result$log2FoldChange), na.last = TRUE), , drop = FALSE]
    result_file <- file.path(manifest$output_dir[[i]], paste0("pseudobulk_DESeq2__", manifest$file_slug[[i]], ".tsv"))
    normalized_file <- file.path(manifest$output_dir[[i]], paste0("pseudobulk_normalized_counts__", manifest$file_slug[[i]], ".tsv"))
    utils::write.table(result, result_file, sep = "\t", row.names = FALSE, quote = FALSE)
    normalized <- as.data.frame(DESeq2::counts(dds, normalized = TRUE), check.names = FALSE)
    utils::write.table(data.frame(gene = rownames(normalized), normalized, check.names = FALSE), normalized_file, sep = "\t", row.names = FALSE, quote = FALSE)
    plot_error <- tryCatch({ save_pseudobulk_de_plots(result, normalized, design, comparison, reference, manifest$file_slug[[i]], manifest$output_dir[[i]]); "" }, error = function(e) conditionMessage(e))
    writeLines(c("Primary inference: sample-level pseudobulk DESeq2", paste0("Comparison: ", comparison, " vs ", reference), paste0("Population: ", manifest$annotation_value[[i]]), paste0("Replicates: ", paste(names(replicates), as.integer(replicates), sep = "=", collapse = ", ")), paste0("Model: ", paste(deparse(design_formula), collapse = "")), "Cell-level Wilcoxon results compare cells directly; cells from the same biological sample are not independent replicates."), file.path(manifest$output_dir[[i]], paste0("methods__", manifest$file_slug[[i]], ".txt")))
    manifest$pseudobulk_status[[i]] <- "complete"; manifest$pseudobulk_note[[i]] <- if (nzchar(plot_error)) paste0("DE completed; figure generation skipped: ", plot_error) else ""
    manifest$pseudobulk_result_file[[i]] <- result_file; manifest$normalized_counts_file[[i]] <- normalized_file
    if (tolower(manifest$annotation_value[[i]]) %in% c("all", "all cells")) {
      utils::write.table(result, file.path(tables, "pseudobulk_differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      utils::write.table(data.frame(gene = rownames(normalized), normalized, check.names = FALSE), file.path(tables, "pseudobulk_normalized_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
    successes <- successes + 1L
    ""
  }, error = function(e) conditionMessage(e))
  if (nzchar(error)) { manifest$pseudobulk_status[[i]] <- "skipped"; manifest$pseudobulk_note[[i]] <- error }
}
utils::write.table(manifest, file.path(tables, "differential_results_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!successes && identical(method, "pseudobulk")) stop("No population met the pseudobulk requirements. Review differential_results_manifest.tsv.")
if (!successes) warning("No population met pseudobulk requirements; cell-level Wilcoxon outputs may still be available.")
