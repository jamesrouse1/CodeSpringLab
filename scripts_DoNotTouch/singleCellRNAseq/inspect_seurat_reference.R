args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: inspect_seurat_reference.R reference.rda|rds output.tsv")

reference_path <- args[[1]]
output_path <- args[[2]]
if (!file.exists(reference_path)) stop("Seurat reference file was not found: ", reference_path)

extension <- tolower(tools::file_ext(reference_path))
if (identical(extension, "rds")) {
  reference <- readRDS(reference_path)
  object_name <- basename(reference_path)
} else if (identical(extension, "rda")) {
  environment <- new.env(parent = emptyenv())
  loaded_names <- load(reference_path, envir = environment)
  seurat_names <- loaded_names[vapply(loaded_names, function(name) inherits(environment[[name]], "Seurat"), logical(1))]
  if (!length(seurat_names)) stop("The .rda file does not contain a Seurat object.")
  object_name <- seurat_names[[1]]
  reference <- environment[[object_name]]
} else {
  stop("Seurat reference inspection accepts only .rda or .rds files.")
}

if (!inherits(reference, "Seurat")) stop("The selected file is not a Seurat object.")
reference <- SeuratObject::UpdateSeuratObject(reference)

summarize_labels <- function(values) {
  values <- trimws(as.character(values))
  valid <- !is.na(values) & nzchar(values)
  labels <- sort(unique(values[valid]))
  list(labels = labels, non_missing = sum(valid), usable = sum(valid) >= 20L && length(labels) >= 2L && length(labels) <= 500L)
}

active <- summarize_labels(SeuratObject::Idents(reference))
rows <- list(data.frame(
  value = "",
  source = "Active identities",
  label_count = length(active$labels),
  non_missing_cells = active$non_missing,
  preview = paste(utils::head(active$labels, 8L), collapse = " | "),
  stringsAsFactors = FALSE
))

for (column in colnames(reference@meta.data)) {
  summary <- summarize_labels(reference@meta.data[[column]])
  if (!summary$usable) next
  rows[[length(rows) + 1L]] <- data.frame(
    value = column,
    source = column,
    label_count = length(summary$labels),
    non_missing_cells = summary$non_missing,
    preview = paste(utils::head(summary$labels, 8L), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, rows)
result$reference_object <- object_name
result$reference_cells <- ncol(reference)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(output_path, ".tmp.", Sys.getpid())
utils::write.table(result, temporary, sep = "\t", row.names = FALSE, quote = FALSE)
if (!file.rename(temporary, output_path)) stop("Could not finalize reference metadata table: ", output_path)
message("Reference label choices written: ", output_path)
