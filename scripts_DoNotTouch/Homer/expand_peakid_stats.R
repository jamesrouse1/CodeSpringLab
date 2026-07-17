args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: expand_peakid_stats.R INPUT [OUTPUT]")

input <- args[[1]]
output <- if (length(args) >= 2L) args[[2]] else input
if (!file.exists(input) || file.info(input)$size <= 0) stop("Annotation table is missing or empty: ", input)

table <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
peak_columns <- names(table)[grepl("^PeakID", names(table), ignore.case = TRUE)]
if (!length(peak_columns) || !nrow(table)) {
  if (!identical(normalizePath(input, mustWork = TRUE), normalizePath(output, mustWork = FALSE))) {
    if (!file.copy(input, output, overwrite = TRUE)) stop("Could not copy annotation table to: ", output)
  }
  quit(save = "no", status = 0L)
}

peak_column <- peak_columns[[1]]
parsed <- lapply(strsplit(as.character(table[[peak_column]]), "|", fixed = TRUE), function(parts) {
  fields <- parts[grepl("=", parts, fixed = TRUE)]
  if (!length(fields)) return(character(0))
  positions <- regexpr("=", fields, fixed = TRUE)
  keys <- trimws(substr(fields, 1L, positions - 1L))
  values <- substr(fields, positions + 1L, nchar(fields))
  keep <- nzchar(keys) & !duplicated(keys)
  stats::setNames(values[keep], keys[keep])
})

stat_names <- unique(unlist(lapply(parsed, names), use.names = FALSE))
stat_names <- stat_names[nzchar(stat_names)]
added_columns <- FALSE
for (name in stat_names) {
  if (name %in% names(table)) next
  values <- vapply(parsed, function(row) if (name %in% names(row)) row[[name]] else NA_character_, character(1))
  table[[name]] <- type.convert(values, as.is = TRUE, na.strings = c("NA", "NaN", ""))
  added_columns <- TRUE
}
if (!added_columns) quit(save = "no", status = 0L)

destination_dir <- dirname(output)
dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
temporary <- tempfile("expanded_peak_stats_", tmpdir = destination_dir)
on.exit(unlink(temporary), add = TRUE)
write.table(table, temporary, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
if (!file.rename(temporary, output)) {
  if (!file.copy(temporary, output, overwrite = TRUE)) stop("Could not write expanded annotation table: ", output)
  unlink(temporary)
}
