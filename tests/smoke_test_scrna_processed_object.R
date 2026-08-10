#!/usr/bin/env Rscript

required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  message("Processed-object continuation smoke test skipped: Seurat runtime unavailable.")
  quit(save = "no", status = 0L)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "tests"), ".."), mustWork = TRUE)
runner <- file.path(repo, "scripts_DoNotTouch", "singleCellRNAseq", "scrna_pipeline_seurat.R")
work <- tempfile("codespring_processed_object_")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

suppressPackageStartupMessages(library(Seurat))
set.seed(17)
counts <- matrix(
  stats::rpois(240L * 80L, lambda = 2), nrow = 240L,
  dimnames = list(c(paste0("MT-", seq_len(8L)), paste0("Gene", 9:240)), paste0("cell", seq_len(80L)))
)
object <- CreateSeuratObject(counts)
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, nfeatures = 120L, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, npcs = 8L, verbose = FALSE)
coordinates <- matrix(stats::rnorm(ncol(object) * 2L), ncol = 2L, dimnames = list(colnames(object), c("UMAP_1", "UMAP_2")))
object[["umap"]] <- CreateDimReducObject(embeddings = coordinates, key = "UMAP_", assay = "RNA")
object$seurat_clusters <- factor(rep(c("0", "1"), each = ncol(object) / 2L))
Seurat::Idents(object) <- object$seurat_clusters
object$cell_type <- ifelse(object$seurat_clusters == "0", "Type A", "Type B")

source_object <- file.path(work, "processed.rds")
saveRDS(object, source_object)
source_md5 <- unname(tools::md5sum(source_object))
manifest <- file.path(work, "samples.tsv")
utils::write.table(data.frame(sample_id = "fixture", input_path = source_object), manifest, sep = "\t", row.names = FALSE, quote = FALSE)
params <- file.path(work, "params.tsv")
utils::write.table(
  data.frame(key = c("normalization", "integration", "n_pcs", "min_features", "min_cells_per_gene", "doublet_method", "annotation_name"), value = c("lognormalize", "none", "8", "0", "1", "none", "cell_type")),
  params, sep = "\t", row.names = FALSE, quote = FALSE
)
output <- file.path(work, "output")
run <- function(stage) {
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(runner, manifest, output, params, stage)), stdout = TRUE, stderr = TRUE)
  code <- attr(status, "status") %||% 0L
  if (code != 0L) stop("Processed-object ", stage, " failed:\n", paste(status, collapse = "\n"))
}
run("inspect")
stopifnot(file.exists(file.path(output, "objects", "processed_seurat.rds")))
run("annotate")
stopifnot(
  identical(unname(tools::md5sum(source_object)), source_md5),
  file.exists(file.path(output, "_STAGE_ANNOTATE_COMPLETE")),
  file.exists(file.path(output, "tables", "cell_metadata.tsv")),
  file.exists(file.path(output, "figures", "05_umap_cell_type.png"))
)
message("Processed Seurat object continuation smoke test passed.")
