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
count_path <- file.path(tables, "pseudobulk_counts.tsv")
design_path <- file.path(tables, "pseudobulk_design.tsv")
if (!file.exists(count_path) || !file.exists(design_path)) stop("Pseudobulk count or design table is missing; rerun the differential-expression stage.")
count_table <- utils::read.delim(count_path, check.names = FALSE, stringsAsFactors = FALSE)
design <- utils::read.delim(design_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample_id", "group", "cells") %in% names(design)) || !"gene" %in% names(count_table)) stop("Pseudobulk input tables have an invalid schema.")
reference <- get("de_reference"); comparison <- get("de_comparison")
design <- design[design$group %in% c(reference, comparison), , drop = FALSE]
replicates <- table(design$group)
if (!all(c(reference, comparison) %in% names(replicates)) || any(replicates[c(reference, comparison)] < 2L)) {
  stop("Pseudobulk DE requires at least two independent sample-level replicates in each comparison group. Cell-level output, if requested, remains exploratory.")
}
count_columns <- intersect(as.character(design$sample_id), names(count_table))
if (length(count_columns) != NROW(design)) stop("Pseudobulk count columns do not match the sample-level design.")
counts <- as.matrix(count_table[, count_columns, drop = FALSE]); storage.mode(counts) <- "numeric"
rownames(counts) <- make.unique(as.character(count_table$gene))
counts <- round(counts)
keep <- rowSums(counts >= 10) >= 2L & rowSums(counts) >= 20L
if (sum(keep) < 10L) stop("Too few expressed genes remain for pseudobulk DE after independent low-count filtering.")
counts <- counts[keep, , drop = FALSE]
design <- design[match(colnames(counts), design$sample_id), , drop = FALSE]
design$group <- stats::relevel(factor(design$group), ref = reference)
rownames(design) <- design$sample_id
covariates <- Filter(nzchar, trimws(strsplit(get("de_covariates", ""), ",", fixed = TRUE)[[1]]))
missing_covariates <- setdiff(covariates, names(design))
if (length(missing_covariates)) stop("Pseudobulk design is missing covariate(s): ", paste(missing_covariates, collapse = ", "))
for (covariate in covariates) design[[covariate]] <- factor(design[[covariate]])
design_formula <- stats::reformulate(c(covariates, "group"))
model_matrix <- stats::model.matrix(design_formula, design)
if (qr(model_matrix)$rank < ncol(model_matrix)) stop("The pseudobulk model is not full rank. Remove a covariate that is confounded with the comparison group.")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = design, design = design_formula)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
result <- as.data.frame(DESeq2::results(dds, contrast = c("group", comparison, reference), independentFiltering = TRUE))
result <- data.frame(gene = rownames(result), comparison = paste0(comparison, "_vs_", reference), result, analysis_level = "sample-level pseudobulk DESeq2", check.names = FALSE)
result <- result[order(result$padj, -abs(result$log2FoldChange), na.last = TRUE), , drop = FALSE]
utils::write.table(result, file.path(tables, "pseudobulk_differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
normalized <- as.data.frame(DESeq2::counts(dds, normalized = TRUE), check.names = FALSE)
utils::write.table(data.frame(gene = rownames(normalized), normalized, check.names = FALSE), file.path(tables, "pseudobulk_normalized_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c(
  "Primary inference: sample-level pseudobulk DESeq2",
  paste0("Comparison: ", comparison, " vs ", reference),
  paste0("Replicates: ", paste(names(replicates), as.integer(replicates), sep = "=", collapse = ", ")),
  paste0("Model: ", paste(deparse(design_formula), collapse = "")),
  "Cell-level Wilcoxon results are exploratory because cells from the same biological sample are not independent replicates."
), file.path(tables, "differential_expression_methods.txt"))
