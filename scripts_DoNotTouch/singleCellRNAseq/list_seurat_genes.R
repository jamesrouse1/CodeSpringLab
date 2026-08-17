#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: list_seurat_genes.R <processed.rds> <output.tsv>")
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required to read the processed object.")

object_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- path.expand(args[[2]])
obj <- readRDS(object_path)
if (is.list(obj) && inherits(obj$object, "Seurat")) obj <- obj$object
if (!inherits(obj, "Seurat")) stop("Processed RDS is not a Seurat object.")

# RNA retains the complete normalized feature space in CodeSpringLab outputs.
# SCT is used only as a fallback for compatible external/legacy objects.
assays <- unique(c("RNA", "SCT", names(obj@assays)))
genes <- character(0)
for (assay in assays) {
  if (!assay %in% names(obj@assays)) next
  values <- rownames(obj[[assay]])
  values <- unique(trimws(as.character(values)))
  values <- values[nzchar(values)]
  if (length(values)) {
    genes <- values
    break
  }
}
if (!length(genes)) stop("No gene names were found in the processed Seurat object.")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(data.frame(gene = sort(genes), stringsAsFactors = FALSE), output_path,
  sep = "\t", row.names = FALSE, quote = FALSE)
