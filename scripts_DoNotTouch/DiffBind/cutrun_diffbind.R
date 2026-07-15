#!/usr/bin/env Rscript

.libPaths(c("/grid/bsr/data/data/oeldemer/4.3", .libPaths()))

suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop("Usage: cutrun_diffbind.R <resolved_samples.tsv> <out_dir> <reference_condition> <min_replicates> <genome>")
}

sample_file <- normalizePath(args[[1]], mustWork = TRUE)
out_root <- normalizePath(args[[2]], mustWork = FALSE)
reference_condition <- args[[3]]
min_replicates <- suppressWarnings(as.integer(args[[4]]))
genome <- tolower(args[[5]])
if (!is.finite(min_replicates) || min_replicates < 1L) min_replicates <- 2L
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

slug <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unspecified")
}

write_tsv <- function(x, path, row.names = FALSE) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = row.names, col.names = TRUE)
}

read_seacr <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE)
  if (ncol(x) < 3L) stop("SEACR file has fewer than three columns: ", path)
  start <- suppressWarnings(as.integer(x[[2]]))
  end <- suppressWarnings(as.integer(x[[3]]))
  keep <- !is.na(start) & !is.na(end) & end > start & nzchar(as.character(x[[1]]))
  GRanges(seqnames = as.character(x[[1]][keep]), ranges = IRanges(start = start[keep] + 1L, end = end[keep]))
}

write_bed <- function(gr, path, extra = NULL) {
  out <- data.frame(
    chrom = as.character(seqnames(gr)),
    start = start(gr) - 1L,
    end = end(gr),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra)) out <- cbind(out, extra)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

condition_consensus <- function(peaks, support_required) {
  merged <- reduce(do.call(c, unname(peaks)), ignore.strand = TRUE)
  support <- Reduce(`+`, lapply(peaks, function(gr) countOverlaps(merged, gr, ignore.strand = TRUE) > 0L))
  keep <- support >= support_required
  list(peaks = merged[keep], support = support[keep])
}

safe_png <- function(path, expr, width = 1100, height = 850) {
  png(path, width = width, height = height, res = 140)
  on.exit(dev.off(), add = TRUE)
  value <- try(force(expr), silent = TRUE)
  if (!inherits(value, "try-error") && inherits(value, "trellis")) print(value)
  invisible(value)
}

samples <- read.delim(sample_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
required <- c("SampleID", "CellType", "Mark", "Condition", "Replicate", "bamReads", "Peaks", "Spikein", "normalization_mode")
missing_columns <- setdiff(required, names(samples))
if (length(missing_columns)) stop("Resolved CUT&RUN sample sheet is missing columns: ", paste(missing_columns, collapse = ", "))

for (field in required) samples[[field]] <- trimws(as.character(samples[[field]]))
samples$CellType[!nzchar(samples$CellType)] <- "all"
samples$analysis_group <- paste(samples$CellType, samples$Mark, sep = "__")

run_rows <- list()
completed <- 0L

record_run <- function(cell_type, mark, comparison, status, samples_n, peaks_n = NA_integer_, normalization = "", message = "", directory = "") {
  run_rows[[length(run_rows) + 1L]] <<- data.frame(
    cell_type = cell_type,
    mark = mark,
    comparison = comparison,
    reference = reference_condition,
    status = status,
    samples = samples_n,
    consensus_peaks = peaks_n,
    normalization = normalization,
    directory = directory,
    message = message,
    stringsAsFactors = FALSE
  )
}

for (group_name in unique(samples$analysis_group)) {
  group <- samples[samples$analysis_group == group_name, , drop = FALSE]
  cell_type <- group$CellType[[1]]
  mark <- group$Mark[[1]]
  conditions <- unique(group$Condition)

  if (!reference_condition %in% conditions) {
    record_run(cell_type, mark, "", "skipped", nrow(group), message = paste("Reference condition", reference_condition, "is absent"))
    next
  }

  comparisons <- setdiff(conditions, reference_condition)
  if (!length(comparisons)) {
    record_run(cell_type, mark, "", "skipped", nrow(group), message = "Only the reference condition is present")
    next
  }

  for (comparison in comparisons) {
    analysis_samples <- group[group$Condition %in% c(reference_condition, comparison), , drop = FALSE]
    group_sizes <- table(analysis_samples$Condition)
    required_replicates <- max(2L, min_replicates)
    if (any(!c(reference_condition, comparison) %in% names(group_sizes)) || any(group_sizes[c(reference_condition, comparison)] < required_replicates)) {
      record_run(cell_type, mark, comparison, "skipped", nrow(analysis_samples), message = paste("At least", required_replicates, "biological replicates are required in both conditions"))
      next
    }

    run_slug <- paste(slug(cell_type), slug(mark), paste0(slug(comparison), "_vs_", slug(reference_condition)), sep = "__")
    out_dir <- file.path(out_root, run_slug)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    result <- tryCatch({
      peaksets <- setNames(lapply(analysis_samples$Peaks, read_seacr), analysis_samples$SampleID)
      consensus_by_condition <- list()
      for (condition in c(reference_condition, comparison)) {
        ids <- analysis_samples$SampleID[analysis_samples$Condition == condition]
        support_required <- min_replicates
        consensus <- condition_consensus(peaksets[ids], support_required)
        if (!length(consensus$peaks)) stop("No reproducible SEACR peaks remained for condition ", condition, " at replicate threshold ", support_required)
        consensus_by_condition[[condition]] <- consensus$peaks
        write_bed(
          consensus$peaks,
          file.path(out_dir, paste0(slug(condition), ".consensus_peaks.bed")),
          data.frame(replicate_support = consensus$support)
        )
      }

      master <- reduce(do.call(c, consensus_by_condition), ignore.strand = TRUE)
      write_bed(master, file.path(out_dir, "master_consensus_peaks.bed"))

      sheet <- data.frame(
        SampleID = analysis_samples$SampleID,
        Condition = analysis_samples$Condition,
        Replicate = analysis_samples$Replicate,
        Factor = analysis_samples$Mark,
        Tissue = analysis_samples$CellType,
        bamReads = analysis_samples$bamReads,
        Peaks = analysis_samples$Peaks,
        PeakCaller = "raw",
        ScoreCol = 4L,
        stringsAsFactors = FALSE
      )

      use_spikein <- all(tolower(analysis_samples$normalization_mode) == "spikein") && all(nzchar(analysis_samples$Spikein)) && all(file.exists(analysis_samples$Spikein))
      mixed_spikein <- any(tolower(analysis_samples$normalization_mode) == "spikein") && !use_spikein
      if (mixed_spikein) stop("Spike-in normalization is incomplete: every sample in a comparison must have a readable E. coli BAM")
      if (use_spikein) sheet$Spikein <- analysis_samples$Spikein
      write_tsv(sheet, file.path(out_dir, "diffbind_sample_sheet.tsv"))

      db <- dba(sampleSheet = sheet)
      if (genome %in% c("human", "hg38", "grch38")) {
        db <- dba.blacklist(db, blacklist = DBA_BLACKLIST_HG38, greylist = FALSE)
      }

      db <- dba.count(db, peaks = master, summits = FALSE, bSubControl = FALSE)
      info <- dba.show(db)
      write_tsv(info, file.path(out_dir, "sample_count_qc.tsv"), row.names = FALSE)

      db <- if (use_spikein) {
        dba.normalize(db, method = DBA_DESEQ2, normalize = DBA_NORM_LIB, spikein = TRUE)
      } else {
        dba.normalize(db, method = DBA_DESEQ2, normalize = DBA_NORM_LIB, background = TRUE)
      }
      normalization_label <- if (use_spikein) "E. coli spike-in" else "genomic background bins"
      norm <- dba.normalize(db, method = DBA_DESEQ2, bRetrieve = TRUE)
      norm_table <- data.frame(
        SampleID = info$ID,
        library_size = norm$lib.sizes,
        normalization_factor = norm$norm.factors,
        normalization = normalization_label,
        stringsAsFactors = FALSE
      )
      write_tsv(norm_table, file.path(out_dir, "normalization_factors.tsv"))

      db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2L, reorderMeta = list(Condition = reference_condition))
      db <- dba.analyze(db, method = DBA_DESEQ2)
      write_tsv(dba.show(db, bContrasts = TRUE), file.path(out_dir, "contrasts.tsv"))

      report <- dba.report(db, contrast = 1L, method = DBA_DESEQ2, th = 1, bCounts = TRUE)
      report_df <- as.data.frame(report)
      write_tsv(report_df, file.path(out_dir, "all_differential_peaks.tsv"))
      significant <- report_df[!is.na(report_df$FDR) & report_df$FDR <= 0.05, , drop = FALSE]
      write_tsv(significant, file.path(out_dir, "significant_differential_peaks.tsv"))
      if (nrow(significant)) {
        sig_gr <- report[!is.na(report$FDR) & report$FDR <= 0.05]
        sig_name <- paste0(seqnames(sig_gr), ":", start(sig_gr), "-", end(sig_gr), "|Fold=", signif(sig_gr$Fold, 5), "|FDR=", signif(sig_gr$FDR, 5))
        write_bed(sig_gr, file.path(out_dir, "significant_differential_peaks.bed"), data.frame(name = sig_name, Fold = sig_gr$Fold))
      } else {
        file.create(file.path(out_dir, "significant_differential_peaks.bed"))
      }

      counted <- try(dba.peakset(db, bRetrieve = TRUE), silent = TRUE)
      if (!inherits(counted, "try-error")) write_tsv(as.data.frame(counted), file.path(out_dir, "consensus_peak_counts.tsv"))

      safe_png(file.path(out_dir, "correlation_heatmap.png"), plot(db))
      safe_png(file.path(out_dir, "pca_normalized_counts.png"), dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID))
      safe_png(file.path(out_dir, "volcano_differential_peaks.png"), dba.plotVolcano(db, contrast = 1L, method = DBA_DESEQ2))
      safe_png(file.path(out_dir, "ma_differential_peaks.png"), dba.plotMA(db, contrast = 1L, method = DBA_DESEQ2))
      safe_png(file.path(out_dir, "differential_peak_heatmap.png"), dba.plotHeatmap(db, contrast = 1L, method = DBA_DESEQ2, th = 1, correlations = FALSE))

      saveRDS(db, file.path(out_dir, "diffbind_object.rds"))
      list(peaks = length(master), normalization = normalization_label)
    }, error = function(e) e)

    if (inherits(result, "error")) {
      writeLines(conditionMessage(result), file.path(out_dir, "ERROR.txt"))
      record_run(cell_type, mark, comparison, "failed", nrow(analysis_samples), message = conditionMessage(result), directory = out_dir)
    } else {
      completed <- completed + 1L
      record_run(cell_type, mark, comparison, "complete", nrow(analysis_samples), peaks_n = result$peaks, normalization = result$normalization, directory = out_dir)
    }
  }
}

summary <- if (length(run_rows)) do.call(rbind, run_rows) else data.frame(
  cell_type = character(), mark = character(), comparison = character(), reference = character(), status = character(),
  samples = integer(), consensus_peaks = integer(), normalization = character(), directory = character(), message = character()
)
write_tsv(summary, file.path(out_root, "cutrun_diffbind_summary.tsv"))

if (completed < 1L) stop("No CUT&RUN differential comparison completed. See cutrun_diffbind_summary.tsv for skipped or failed groups.")
writeLines(as.character(Sys.time()), file.path(out_root, "_COMPLETE"))
