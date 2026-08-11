#!/usr/bin/env Rscript

if (!requireNamespace("Matrix", quietly = TRUE) || !requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("scDblFinder", quietly = TRUE)) {
  message("Capture-aware scDblFinder smoke test skipped: runtime unavailable.")
  quit(save = "no", status = 0L)
}

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
pipeline <- file.path(repo, "scripts_DoNotTouch", "singleCellRNAseq", "scrna_pipeline_seurat.R")
work <- tempfile("codespring-capture-doublets-")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

set.seed(84)
genes <- c(paste0("Gene", seq_len(296)), "MT-Co1", "MT-Co2", "Rpl1", "Rps1")
write_10x <- function(sample_id, shift) {
  path <- file.path(work, sample_id)
  dir.create(path)
  counts <- matrix(rpois(length(genes) * 150L, lambda = 1.2), nrow = length(genes), dimnames = list(genes, paste0("BC", seq_len(150L), "-1")))
  counts[seq_len(30L), ] <- counts[seq_len(30L), ] + shift
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), file.path(path, "matrix.mtx"))
  utils::write.table(data.frame(id = paste0("ENSG", seq_along(genes)), gene = genes, type = "Gene Expression"), file.path(path, "features.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  utils::write.table(colnames(counts), file.path(path, "barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  path
}

manifest <- data.frame(
  sample_id = c("donor_A", "donor_B"), capture_id = c("capture_1", "capture_1"),
  input_path = c(write_10x("donor_A", 0), write_10x("donor_B", 1)), condition = c("control", "treated")
)
manifest_path <- file.path(work, "samples.tsv")
utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
params <- data.frame(
  key = c("normalization", "integration", "min_features", "min_counts", "max_features", "max_percent_mt", "min_cells_per_gene", "doublet_method", "doublet_rate", "remove_doublets", "seed"),
  value = c("lognormalize", "none", "10", "0", "0", "100", "1", "scdblfinder", "0", "false", "1234")
)
params_path <- file.path(work, "params.tsv")
utils::write.table(params, params_path, sep = "\t", row.names = FALSE, quote = FALSE)
out <- file.path(work, "output")

inspect_status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out, params_path, "inspect")))
qc_status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out, params_path, "qc")))
stopifnot(inspect_status == 0L, qc_status == 0L)
summary <- utils::read.delim(file.path(out, "tables", "doublet_summary_by_capture.tsv"), check.names = FALSE)
calls <- utils::read.delim(file.path(out, "tables", "doublet_calls.tsv"), check.names = FALSE)
stopifnot(
  NROW(summary) == 1L, summary$capture_id[[1]] == "capture_1", summary$rate_source[[1]] == "automatic",
  identical(sort(unique(calls$sample_id)), c("donor_A", "donor_B")), unique(calls$capture_id) == "capture_1",
  NROW(calls) == 300L
)
cat("SCRNA_CAPTURE_DOUBLETS_OK\n")
