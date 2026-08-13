#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
pipeline <- file.path(repo, "scripts_DoNotTouch", "singleCellRNAseq", "scrna_pipeline_seurat.R")
stopifnot(file.exists(pipeline))

# Load only the two helpers under test without executing the command-line
# pipeline. This protects compatibility with published Seurat v3/v4 reference
# objects whose RNA assay remains a legacy Assay after UpdateSeuratObject().
expressions <- as.list(parse(file = pipeline))
helper_environment <- new.env(parent = globalenv())
for (name in c("assay_layers_safe", "ensure_log_normalized_rna")) {
  definition <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) && identical(as.character(expression[[2]]), name)
  }, expressions)
  stopifnot(length(definition) == 1L)
  eval(definition[[1]], envir = helper_environment)
}

set.seed(42)
counts <- matrix(
  stats::rpois(120L * 50L, lambda = 2), nrow = 120L,
  dimnames = list(paste0("Gene", seq_len(120L)), paste0("Cell", seq_len(50L)))
)
reference <- Seurat::CreateSeuratObject(counts = counts)
reference[["RNA"]] <- SeuratObject::CreateAssayObject(counts = counts)
stopifnot(inherits(reference[["RNA"]], "Assay"), !inherits(reference[["RNA"]], "Assay5"))

reference <- helper_environment$ensure_log_normalized_rna(reference, "Legacy reference")
stopifnot(
  inherits(reference[["RNA"]], "Assay"),
  "data" %in% SeuratObject::Layers(reference[["RNA"]]),
  nrow(Seurat::GetAssayData(reference, assay = "RNA", layer = "data")) == nrow(counts)
)

cat("SCRNA_LEGACY_ASSAY_REFERENCE_OK\n")
