#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: extract_seurat_marker_expression.R <processed.rds> <gene> <output.tsv.gz>")
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required to extract normalized marker expression.")
object_path <- normalizePath(args[[1]], mustWork = TRUE)
gene <- args[[2]]
output_path <- path.expand(args[[3]])
obj <- readRDS(object_path)
if (!inherits(obj, "Seurat")) stop("Processed RDS is not a Seurat object.")
assay <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
Seurat::DefaultAssay(obj) <- assay
if (!gene %in% rownames(obj[[assay]])) stop("Gene is absent from the processed Seurat object: ", gene)
fetched <- Seurat::FetchData(obj, vars = gene, layer = "data", clean = FALSE)
result <- data.frame(cell = rownames(fetched), expression = as.numeric(fetched[, 1]), stringsAsFactors = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_path, open = "wt")
utils::write.table(result, connection, sep = "\t", row.names = FALSE, quote = FALSE)
close(connection)
