#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: extract_seurat_marker_expression.R <processed.rds> <gene> <output.tsv.gz>")
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required to extract normalized marker expression.")
object_path <- normalizePath(args[[1]], mustWork = TRUE)
gene <- args[[2]]
output_path <- path.expand(args[[3]])
obj <- readRDS(object_path)
if (!inherits(obj, "Seurat")) stop("Processed RDS is not a Seurat object.")

# FetchData(layer = "data") is not portable across Seurat v4/v5 and can fail
# for an otherwise valid legacy Assay. Read the normalized expression matrix
# directly, trying the current v5 layer interface and the v4 slot interface.
# Prefer SCT when it has normalized data, then reliably fall back to RNA.
assays <- unique(c(if ("SCT" %in% names(obj@assays)) "SCT", "RNA", names(obj@assays)))
read_normalized_matrix <- function(object, assay_name) {
  attempts <- list(
    function() Seurat::GetAssayData(object, assay = assay_name, layer = "data"),
    function() Seurat::GetAssayData(object, assay = assay_name, slot = "data"),
    function() SeuratObject::LayerData(object[[assay_name]], layer = "data")
  )
  for (attempt in attempts) {
    value <- tryCatch(attempt(), error = function(e) NULL)
    if (!is.null(value) && nrow(value) > 0L && ncol(value) > 0L) return(value)
  }
  NULL
}

expression <- NULL
cells <- character(0)
used_assay <- ""
for (assay in assays) {
  if (!assay %in% names(obj@assays) || !gene %in% rownames(obj[[assay]])) next
  matrix <- read_normalized_matrix(obj, assay)
  if (is.null(matrix) || !gene %in% rownames(matrix)) next
  expression <- as.numeric(matrix[gene, , drop = TRUE])
  cells <- colnames(matrix)
  used_assay <- assay
  break
}
if (is.null(expression) || !length(cells)) {
  stop("Normalized data could not be read for gene ", gene, " from the available Seurat assays.")
}
if (length(expression) != length(cells)) stop("Normalized expression length does not match the cell barcodes.")
result <- data.frame(cell = cells, expression = expression, stringsAsFactors = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_path, open = "wt")
utils::write.table(result, connection, sep = "\t", row.names = FALSE, quote = FALSE)
close(connection)
