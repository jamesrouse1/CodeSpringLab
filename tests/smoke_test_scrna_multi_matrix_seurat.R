#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
pipeline <- file.path(repo, "scripts_DoNotTouch", "singleCellRNAseq", "scrna_pipeline_seurat.R")
stopifnot(file.exists(pipeline), requireNamespace("Seurat", quietly = TRUE), requireNamespace("harmony", quietly = TRUE))

work <- tempfile("codespring-scrna-multimatrix-")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

set.seed(42)
genes <- c(paste0("Gene", seq_len(296)), "MT-Co1", "MT-Co2", "Rpl1", "Rps1")
write_10x <- function(sample_id, batch_shift) {
  path <- file.path(work, sample_id)
  dir.create(path)
  n <- 120L
  counts <- matrix(rpois(length(genes) * n, lambda = 0.8), nrow = length(genes), dimnames = list(genes, paste0("BC", seq_len(n), "-1")))
  counts[1:20, 1:60] <- counts[1:20, 1:60] + rpois(20 * 60, 2)
  counts[21:40, 61:120] <- counts[21:40, 61:120] + rpois(20 * 60, 2)
  counts[41:70, ] <- counts[41:70, ] + batch_shift
  Matrix::writeMM(Matrix(counts, sparse = TRUE), file.path(path, "matrix.mtx"))
  utils::write.table(data.frame(id = paste0("ENSG", seq_along(genes)), gene = genes, type = "Gene Expression"), file.path(path, "features.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  utils::write.table(colnames(counts), file.path(path, "barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  path
}

paths <- c(write_10x("sample_A", 0), write_10x("sample_B", 0), write_10x("sample_C", 1))
manifest <- data.frame(sample_id = c("sample_A", "sample_B", "sample_C"), capture_id = c("capture_1", "capture_1", "capture_2"), input_path = paths, condition = c("control", "treated", "treated"), technical_batch = c("run_1", "run_1", "run_2"))
manifest_path <- file.path(work, "samples.tsv")
utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)

params <- data.frame(
  key = c("normalization", "integration", "batch_column", "cluster_resolution", "min_features", "min_counts", "max_features", "max_percent_mt", "n_pcs", "min_cells_per_gene", "doublet_method", "doublet_rate", "remove_doublets", "seed", "harmony_theta", "harmony_lambda", "harmony_max_iter"),
  value = c("lognormalize", "harmony", "technical_batch", "0.4", "10", "0", "0", "100", "10", "1", "none", "0", "false", "1234", "2", "1", "20")
)
params_path <- file.path(work, "params.tsv")
utils::write.table(params, params_path, sep = "\t", row.names = FALSE, quote = FALSE)
out <- file.path(work, "output")

status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out, params_path, "all")))
stopifnot(status == 0L)
expected <- c(
  "_COMPLETE", "_STAGE_PREPROCESS_COMPLETE", "_STAGE_CLUSTER_COMPLETE",
  "figures/02_preintegration_umap_sample.png", "figures/02_preintegration_umap_batch.png",
  "tables/preintegration_umap_coordinates.tsv", "tables/umap_coordinates.tsv", "tables/qc_recommended_thresholds.tsv",
  "objects/processed_seurat.rds"
)
stopifnot(all(file.exists(file.path(out, expected))))
qc_recommended <- utils::read.delim(file.path(out, "tables", "qc_recommended_thresholds.tsv"), check.names = FALSE)
stopifnot(
  NROW(qc_recommended) == 4L,
  identical(qc_recommended$sample_id[[1]], "Recommended global"),
  all(c("min_features", "min_counts", "max_features", "max_percent_mt", "recommendation_basis") %in% names(qc_recommended)),
  all(qc_recommended$max_percent_mt >= 5 & qc_recommended$max_percent_mt <= 50)
)
pre <- utils::read.delim(file.path(out, "tables", "preintegration_umap_coordinates.tsv"), check.names = FALSE)
final <- utils::read.delim(file.path(out, "tables", "umap_coordinates.tsv"), check.names = FALSE)
stopifnot(NROW(pre) == 360L, NROW(final) == 360L, all(c("sample_id", "condition", "technical_batch") %in% names(pre)))
doublets <- utils::read.delim(file.path(out, "tables", "doublet_summary_by_capture.tsv"), check.names = FALSE)
stopifnot(NROW(doublets) == 2L, identical(sort(doublets$capture_id), c("capture_1", "capture_2")), is.na(doublets$expected_doublet_rate[doublets$capture_id == "capture_1"]))
obj <- readRDS(file.path(out, "objects", "processed_seurat.rds"))
stopifnot("harmony" %in% names(obj@reductions), "umap" %in% names(obj@reductions), NROW(obj) == length(genes))

# Anchor integration must group the two inputs from run_1 into one integration
# unit instead of incorrectly correcting sample_A and sample_B apart.
params$value[params$key == "integration"] <- "rpca"
utils::write.table(params, params_path, sep = "\t", row.names = FALSE, quote = FALSE)
out_rpca <- file.path(work, "output_rpca")
status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out_rpca, params_path, "all")))
stopifnot(status == 0L, file.exists(file.path(out_rpca, "_COMPLETE")))
pre_state <- readRDS(file.path(out_rpca, "checkpoints", "03_preprocessed_seurat.rds"))
stopifnot(length(pre_state$objects) == 2L, identical(sort(names(pre_state$objects)), c("run_1", "run_2")))
rpca_obj <- readRDS(file.path(out_rpca, "objects", "processed_seurat.rds"))
stopifnot(all(c("counts", "data") %in% SeuratObject::Layers(rpca_obj[["RNA"]])))

# Exercise the real capture-aware automatic-rate path when scDblFinder is
# available (the official cluster Seurat module includes it).
if (requireNamespace("scDblFinder", quietly = TRUE)) {
  params$value[params$key == "doublet_method"] <- "scdblfinder"
  params$value[params$key == "doublet_rate"] <- "0"
  params$value[params$key == "remove_doublets"] <- "false"
  utils::write.table(params, params_path, sep = "\t", row.names = FALSE, quote = FALSE)
  out_doublet <- file.path(work, "output_doublet")
  inspect_status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out_doublet, params_path, "inspect")))
  qc_status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(pipeline, manifest_path, out_doublet, params_path, "qc")))
  stopifnot(inspect_status == 0L, qc_status == 0L)
  capture_summary <- utils::read.delim(file.path(out_doublet, "tables", "doublet_summary_by_capture.tsv"), check.names = FALSE)
  capture_calls <- utils::read.delim(file.path(out_doublet, "tables", "doublet_calls.tsv"), check.names = FALSE)
  stopifnot(NROW(capture_summary) == 2L, all(capture_summary$rate_source == "automatic"), all(c("sample_id", "capture_id") %in% names(capture_calls)))
}
cat("MULTI_MATRIX_SEURAT_HARMONY_OK\n")
