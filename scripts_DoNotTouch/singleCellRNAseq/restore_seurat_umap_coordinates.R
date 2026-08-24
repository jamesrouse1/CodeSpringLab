#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: restore_seurat_umap_coordinates.R <scrna_output_dir>")

out_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

candidates <- c(
  file.path(out_dir, "objects", "processed_seurat.rds"),
  file.path(out_dir, "checkpoints", "04_clustered_seurat.rds")
)
candidates <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
if (!length(candidates)) stop("No processed or clustered Seurat object was found in ", out_dir, ".")

loaded <- readRDS(candidates[[1]])
obj <- if (inherits(loaded, "Seurat")) loaded else if (is.list(loaded) && inherits(loaded$object, "Seurat")) loaded$object else NULL
if (is.null(obj)) stop("The saved file does not contain a recognizable Seurat object: ", candidates[[1]])

reduction_names <- names(obj@reductions)
final_reduction <- if ("umap" %in% reduction_names) "umap" else {
  hits <- grep("umap", reduction_names, value = TRUE, ignore.case = TRUE)
  hits <- setdiff(hits, grep("unintegrated|pre", hits, value = TRUE, ignore.case = TRUE))
  if (length(hits)) hits[[1]] else ""
}
if (!nzchar(final_reduction)) stop("The saved Seurat object has no final UMAP reduction. Available reductions: ", paste(reduction_names, collapse = ", "))

write_embedding <- function(reduction, destination) {
  coordinates <- as.data.frame(Seurat::Embeddings(obj, reduction = reduction), check.names = FALSE)
  if (NCOL(coordinates) < 2L) stop("UMAP reduction '", reduction, "' has fewer than two coordinates.")
  names(coordinates)[1:2] <- c("UMAP_1", "UMAP_2")
  metadata <- obj@meta.data[rownames(coordinates), , drop = FALSE]
  if ("cell" %in% names(metadata)) names(metadata)[names(metadata) == "cell"] <- "input_cell"
  output <- data.frame(cell = rownames(coordinates), coordinates[, c("UMAP_1", "UMAP_2"), drop = FALSE], metadata, check.names = FALSE)
  temporary <- paste0(destination, ".tmp.", Sys.getpid())
  utils::write.table(output, temporary, sep = "\t", row.names = FALSE, quote = FALSE)
  if (!file.rename(temporary, destination)) stop("Could not publish restored UMAP coordinates: ", destination)
  invisible(destination)
}

write_embedding(final_reduction, file.path(tables_dir, "umap_coordinates.tsv"))
pre_reduction <- grep("umap.*unintegrated|unintegrated.*umap|pre.*umap", reduction_names, value = TRUE, ignore.case = TRUE)
if (length(pre_reduction)) write_embedding(pre_reduction[[1]], file.path(tables_dir, "preintegration_umap_coordinates.tsv"))

writeLines(
  c(
    paste0("source_object\t", normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)),
    paste0("final_reduction\t", final_reduction),
    paste0("cells\t", ncol(obj))
  ),
  file.path(tables_dir, "umap_coordinate_restore_summary.tsv")
)
writeLines("complete", file.path(out_dir, "_STAGE_RESTORE_EMBEDDING_COMPLETE"))
message("Restored interactive UMAP coordinates for ", format(ncol(obj), big.mark = ","), " cells from reduction '", final_reduction, "'.")
