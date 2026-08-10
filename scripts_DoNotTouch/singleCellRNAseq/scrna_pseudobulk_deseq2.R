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
    writeLines(c("Primary inference: sample-level pseudobulk DESeq2", paste0("Comparison: ", comparison, " vs ", reference), paste0("Population: ", manifest$annotation_value[[i]]), paste0("Replicates: ", paste(names(replicates), as.integer(replicates), sep = "=", collapse = ", ")), paste0("Model: ", paste(deparse(design_formula), collapse = "")), "Cell-level Wilcoxon results are exploratory because cells from the same biological sample are not independent replicates."), file.path(manifest$output_dir[[i]], paste0("methods__", manifest$file_slug[[i]], ".txt")))
    manifest$pseudobulk_status[[i]] <- "complete"; manifest$pseudobulk_note[[i]] <- ""
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
if (!successes) warning("No population met pseudobulk requirements; exploratory cell-level outputs may still be available.")
