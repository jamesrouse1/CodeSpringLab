args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: inspect_seurat_reference.R reference.rda|rds output.tsv")

reference_path <- args[[1]]
output_path <- args[[2]]
if (!file.exists(reference_path)) stop("Seurat reference file was not found: ", reference_path)
message("Loading reference (", format(file.info(reference_path)$size / 1024^3, digits = 3), " GiB): ", reference_path)

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
message("Reference loaded; reading identities and metadata without rewriting the object")

# Inspection only needs cell labels. Updating a legacy Seurat object here can
# traverse and rewrite every assay/reduction, which is both slow and unrelated
# to listing labels. The annotation job performs any compatibility work it
# actually needs later.
metadata <- tryCatch(reference@meta.data, error = function(e) NULL)
if (is.null(metadata) || !is.data.frame(metadata)) stop("The Seurat reference has no readable meta.data table.")
active_ident <- tryCatch(reference@active.ident, error = function(e) NULL)
if (is.null(active_ident) || length(active_ident) != nrow(metadata)) {
  active_ident <- tryCatch(SeuratObject::Idents(reference), error = function(e) rep(NA_character_, nrow(metadata)))
}

summarize_labels <- function(values) {
  values <- trimws(as.character(values))
  valid <- !is.na(values) & nzchar(values)
  labels <- sort(unique(values[valid]))
  list(labels = labels, non_missing = sum(valid), usable = sum(valid) >= 20L && length(labels) >= 2L && length(labels) <= 500L)
}

active <- summarize_labels(active_ident)
rows <- list(data.frame(
  value = "",
  source = "Active identities",
  label_count = length(active$labels),
  non_missing_cells = active$non_missing,
  preview = paste(utils::head(active$labels, 8L), collapse = " | "),
  stringsAsFactors = FALSE
))

for (column in colnames(metadata)) {
  summary <- summarize_labels(metadata[[column]])
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
result$reference_cells <- nrow(metadata)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(output_path, ".tmp.", Sys.getpid())
utils::write.table(result, temporary, sep = "\t", row.names = FALSE, quote = FALSE)
if (!file.rename(temporary, output_path)) stop("Could not finalize reference metadata table: ", output_path)
message("Reference label choices written: ", output_path)
