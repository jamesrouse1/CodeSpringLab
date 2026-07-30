#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: cutrun_project_summary.R <data_dir> <out_dir> <peak_calling|differential>")

data_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
kind <- match.arg(tolower(args[[3]]), c("peak_calling", "differential"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(data.frame())
  tryCatch(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = ""), error = function(e) data.frame())
}

read_kv <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- tryCatch(read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = ""), error = function(e) data.frame())
  if (!NROW(x) || NCOL(x) < 2L) return(character(0))
  values <- as.character(x[[2]])
  names(values) <- as.character(x[[1]])
  values
}

line_count <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  value <- tryCatch(system2("wc", c("-l", "--", path), stdout = TRUE, stderr = FALSE), error = function(e) character(0))
  if (!length(value)) return(NA_integer_)
  number <- suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", value[[1]])))
  if (is.na(number)) NA_integer_ else number
}

tabular_rows <- function(path) {
  n <- line_count(path)
  if (is.na(n)) NA_integer_ else max(n - 1L, 0L)
}

first_value <- function(x, name, default = "") {
  if (!name %in% names(x) || !NROW(x)) return(default)
  value <- as.character(x[[name]][[1]])
  if (is.na(value)) default else value
}

as_number <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(value), value, NA_real_)
}

bind_rows_fill <- function(rows) {
  rows <- Filter(NROW, rows)
  if (!length(rows)) return(data.frame())
  fields <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(fields, names(x))
    for (field in missing) x[[field]] <- NA
    x[fields]
  })
  do.call(rbind, rows)
}

write_workbook <- function(path, sheets) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    if (!requireNamespace("writexl", quietly = TRUE)) return(FALSE)
    writexl::write_xlsx(sheets, path)
    return(file.exists(path))
  }
  wb <- openxlsx::createWorkbook()
  title <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold", fontSize = 14)
  header <- openxlsx::createStyle(fontColour = "#1F1F1F", fgFill = "#D9EAF7", textDecoration = "bold", wrapText = TRUE, valign = "center")
  for (sheet_name in names(sheets)) {
    value <- sheets[[sheet_name]]
    openxlsx::addWorksheet(wb, sheet_name, gridLines = FALSE)
    openxlsx::writeData(wb, sheet_name, sheet_name, startRow = 1, startCol = 1, colNames = FALSE)
    openxlsx::addStyle(wb, sheet_name, title, rows = 1, cols = 1)
    openxlsx::writeData(wb, sheet_name, value, startRow = 3, startCol = 1, headerStyle = header)
    openxlsx::freezePane(wb, sheet_name, firstActiveRow = 4)
    if (NCOL(value)) openxlsx::setColWidths(wb, sheet_name, cols = seq_len(NCOL(value)), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  TRUE
}

metadata_from_diffbind_sheets <- function() {
  root <- file.path(data_dir, "cutrun_diffbind")
  sheets <- if (dir.exists(root)) list.files(root, pattern = "^diffbind_sample_sheet\\.tsv$", recursive = TRUE, full.names = TRUE) else character(0)
  rows <- lapply(sheets, function(path) {
    x <- read_tsv(path)
    keep <- intersect(c("SampleID", "Condition", "Replicate", "Factor", "Tissue"), names(x))
    if (!length(keep)) return(data.frame())
    out <- x[keep]
    names(out) <- c("Sample", "Condition", "Replicate", "Mark", "Cell type")[match(names(out), c("SampleID", "Condition", "Replicate", "Factor", "Tissue"))]
    out
  })
  meta <- bind_rows_fill(rows)
  if (!NROW(meta) || !"Sample" %in% names(meta)) return(data.frame())
  meta[!duplicated(meta$Sample), , drop = FALSE]
}

make_peak_calling_summary <- function() {
  meta <- metadata_from_diffbind_sheets()
  seacr_root <- file.path(data_dir, "seacr")
  seacr_files <- if (dir.exists(seacr_root)) list.files(seacr_root, pattern = "_seacr_summary\\.txt$", recursive = TRUE, full.names = TRUE) else character(0)
  seacr_rows <- lapply(seacr_files, function(path) {
    values <- read_kv(path)
    if (!length(values)) return(data.frame())
    sample <- basename(dirname(path)); method <- basename(dirname(dirname(path)))
    if (!grepl("^((spikein|cpm)_non|raw_norm|norm|non)_(stringent|relaxed)$", method)) method <- "legacy"
    data.frame(Sample = sample, `Peak caller / setting` = paste("SEACR", method), `Peak count` = values[["peak_count"]] %||% NA_character_,
      FRiP = values[["frip"]] %||% NA_character_, Normalization = values[["normalization"]] %||% "", Stringency = values[["stringency"]] %||% "",
      `Peak file` = values[["peak_bed"]] %||% "", stringsAsFactors = FALSE, check.names = FALSE)
  })
  macs_root <- file.path(data_dir, "macs2")
  macs_files <- if (dir.exists(macs_root)) list.files(macs_root, pattern = "_macs2_summary\\.txt$", recursive = TRUE, full.names = TRUE) else character(0)
  macs_rows <- lapply(macs_files, function(path) {
    values <- read_kv(path); sample <- basename(dirname(path))
    data.frame(Sample = sample, `Peak caller / setting` = paste("MACS2", values[["peak_type"]] %||% ""), `Peak count` = values[["peak_count"]] %||% NA_character_,
      FRiP = NA_character_, Normalization = "BAMPE internal depth normalization", Stringency = values[["qvalue"]] %||% "",
      `Peak file` = values[["peak_file"]] %||% "", stringsAsFactors = FALSE, check.names = FALSE)
  })
  overlap_root <- file.path(data_dir, "peak_overlap")
  overlap_files <- if (dir.exists(overlap_root)) list.files(overlap_root, pattern = "_summary\\.txt$", recursive = TRUE, full.names = TRUE) else character(0)
  overlap_rows <- lapply(overlap_files, function(path) {
    values <- read_kv(path)
    sample <- values[["sample"]] %||% basename(dirname(path))
    setting <- values[["overlap_name"]] %||% basename(dirname(dirname(path)))
    data.frame(Sample = sample, `Peak caller / setting` = paste("Shared overlap", setting), `Peak count` = values[["overlap_peaks"]] %||% NA_character_,
      FRiP = NA_character_, Normalization = "", Stringency = values[["minimum_reciprocal_overlap"]] %||% "",
      `Peak file` = values[["overlap_bed"]] %||% sub("_summary\\.txt$", ".bed", path), stringsAsFactors = FALSE, check.names = FALSE)
  })
  long <- bind_rows_fill(c(seacr_rows, macs_rows, overlap_rows))
  if (NROW(meta) && NROW(long)) long <- merge(meta, long, by = "Sample", all.y = TRUE, sort = FALSE)
  if (!NROW(long)) long <- data.frame(Message = "No completed CUT&RUN peak-calling summaries were found.", stringsAsFactors = FALSE)
  wide <- if (NROW(long) && all(c("Sample", "Peak caller / setting", "Peak count") %in% names(long))) {
    keys <- unique(long$Sample)
    base <- if (NROW(meta)) meta[match(keys, meta$Sample), , drop = FALSE] else data.frame(Sample = keys, stringsAsFactors = FALSE)
    for (setting in unique(long[["Peak caller / setting"]])) {
      column <- paste0("Peaks: ", setting)
      base[[column]] <- vapply(keys, function(sample) {
        hit <- long[long$Sample == sample & long[["Peak caller / setting"]] == setting, "Peak count"]
        if (!length(hit)) NA_character_ else as.character(tail(hit, 1))
      }, character(1))
    }
    base
  } else long
  list(wide = wide, long = long)
}

make_differential_summary <- function() {
  root <- file.path(data_dir, "cutrun_diffbind")
  summary_files <- if (dir.exists(root)) list.files(root, pattern = "^cutrun_diffbind_summary\\.tsv$", recursive = TRUE, full.names = TRUE) else character(0)
  rows <- lapply(summary_files, function(summary_file) {
    summary <- read_tsv(summary_file)
    if (!NROW(summary)) return(data.frame())
    summary <- summary[1, , drop = FALSE]
    out_dir <- first_value(summary, "directory", dirname(summary_file))
    if (!dir.exists(out_dir)) out_dir <- dirname(summary_file)
    sheet <- read_tsv(file.path(out_dir, "diffbind_sample_sheet.tsv"))
    norm <- read_tsv(file.path(out_dir, "normalization_factors.tsv"))
    qc <- read_tsv(file.path(out_dir, "sample_count_qc.tsv"))
    all_table <- first_value(summary, "annotated_differential_peaks_file", file.path(out_dir, "all_differential_peaks.tsv"))
    if (!file.exists(all_table)) all_table <- file.path(out_dir, "all_differential_peaks.tsv")
    sig_table <- first_value(summary, "significant_annotated_peaks_file", file.path(out_dir, "significant_differential_peaks.tsv"))
    if (!file.exists(sig_table)) sig_table <- file.path(out_dir, "significant_differential_peaks.tsv")
    row <- data.frame(
      `Peak source / method` = first_value(summary, "peak_source"),
      `Cell type` = first_value(summary, "cell_type"), Mark = first_value(summary, "mark"),
      Comparison = paste(first_value(summary, "comparison"), "vs", first_value(summary, "reference")),
      Status = first_value(summary, "status"), Normalization = first_value(summary, "normalization"),
      `Minimum peaks/sample` = first_value(summary, "minimum_peaks_per_sample"),
      `Consensus peaks` = first_value(summary, "consensus_peaks"), `Tested peaks` = tabular_rows(all_table),
      `Significant peaks (FDR <= 0.05)` = if (nzchar(first_value(summary, "significant_peaks"))) first_value(summary, "significant_peaks") else tabular_rows(sig_table),
      `Higher in comparison` = first_value(summary, "increased_in_comparison"), `Higher in reference` = first_value(summary, "increased_in_reference"),
      `Annotated all-peaks file` = all_table, `Annotated significant-peaks file` = sig_table,
      `Annotation status` = first_value(summary, "annotation_status"), stringsAsFactors = FALSE, check.names = FALSE
    )
    if (NROW(sheet)) for (i in seq_len(NROW(sheet))) {
      sample <- as.character(sheet$SampleID[[i]])
      label <- paste(as.character(sheet$Condition[[i]]), paste0("rep", as.character(sheet$Replicate[[i]])), sample, sep = " ")
      peak_file <- as.character(sheet$Peaks[[i]])
      row[[paste0(label, " | caller peaks")]] <- line_count(peak_file)
      qc_hit <- if ("ID" %in% names(qc)) qc[as.character(qc$ID) == sample, , drop = FALSE] else data.frame()
      norm_hit <- if ("SampleID" %in% names(norm)) norm[as.character(norm$SampleID) == sample, , drop = FALSE] else data.frame()
      row[[paste0(label, " | counted reads")]] <- if (NROW(qc_hit) && "Reads" %in% names(qc_hit)) qc_hit$Reads[[1]] else NA
      row[[paste0(label, " | FRiP")]] <- if (NROW(qc_hit) && "FRiP" %in% names(qc_hit)) qc_hit$FRiP[[1]] else NA
      row[[paste0(label, " | library size")]] <- if (NROW(norm_hit) && "library_size" %in% names(norm_hit)) norm_hit$library_size[[1]] else NA
      row[[paste0(label, " | normalization factor")]] <- if (NROW(norm_hit) && "normalization_factor" %in% names(norm_hit)) norm_hit$normalization_factor[[1]] else NA
    }
    row
  })
  bind_rows_fill(rows)
}

if (identical(kind, "peak_calling")) {
  result <- make_peak_calling_summary()
  write.table(result$wide, file.path(out_dir, "peak_calling_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result$long, file.path(out_dir, "peak_calling_by_caller.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write_workbook(file.path(out_dir, "peak_calling_summary.xlsx"), list(`All Target Samples` = result$wide, `By Caller and Setting` = result$long))
} else {
  result <- make_differential_summary()
  write.table(result, file.path(out_dir, "differential_peak_comparison_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write_workbook(file.path(out_dir, "differential_peak_comparison_summary.xlsx"), list(`Differential Comparisons` = result))
}

writeLines(as.character(Sys.time()), file.path(out_dir, paste0(kind, "_summary_COMPLETE")))
