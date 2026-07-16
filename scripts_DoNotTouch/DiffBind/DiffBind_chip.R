args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) stop("Usage: DiffBind_chip.R OUTDIR SAMPLE_SHEET REFERENCE COMPARISON GENOME BLACKLIST")

outdir <- args[[1]]
sample_sheet_path <- args[[2]]
reference <- args[[3]]
comparison <- args[[4]]
genome <- tolower(args[[5]])
blacklist <- args[[6]]

.libPaths(c("/grid/bsr/data/data/oeldemer/4.3", .libPaths()))
suppressPackageStartupMessages(library(DiffBind))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
samples <- read.delim(sample_sheet_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("SampleID", "Condition", "Replicate", "bamReads", "Peaks", "PeakCaller")
missing_columns <- setdiff(required, names(samples))
if (length(missing_columns)) stop("ChIP DiffBind sample sheet is missing: ", paste(missing_columns, collapse = ", "))
if (anyDuplicated(samples$SampleID)) stop("ChIP DiffBind SampleID values must be unique.")
if (any(!file.exists(samples$bamReads))) stop("Missing ChIP BAM files: ", paste(samples$bamReads[!file.exists(samples$bamReads)], collapse = ", "))
if (any(!file.exists(samples$Peaks))) stop("Missing ChIP peak files: ", paste(samples$Peaks[!file.exists(samples$Peaks)], collapse = ", "))
group_counts <- table(samples$Condition)
if (!all(c(reference, comparison) %in% names(group_counts)) || any(group_counts[c(reference, comparison)] < 2L)) {
  stop("ChIP DiffBind requires at least two target replicates in both selected conditions.")
}

db <- dba(sampleSheet = samples)
png(file.path(outdir, "diffbind_CorrHeatmap_byPeaks.png"), width = 1000, height = 800, res = 150)
try(plot(db), silent = TRUE)
dev.off()

if (identical(genome, "human")) {
  db <- dba.blacklist(db, blacklist = DBA_BLACKLIST_HG38, greylist = FALSE)
} else if (identical(genome, "mouse") && nzchar(blacklist) && file.exists(blacklist)) {
  db <- dba.blacklist(db, blacklist = blacklist, greylist = FALSE)
}

db <- dba.count(db)
png(file.path(outdir, "diffbind_CorrHeatmap_byCounts.png"), width = 1000, height = 800, res = 150)
try(plot(db), silent = TRUE)
dev.off()

db <- dba.normalize(db)
db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2L, reorderMeta = list(Condition = reference))
db <- dba.analyze(db)
report <- as.data.frame(dba.report(db, contrast = 1))
prefix <- paste0("DifferentialPeaks_", comparison, "_vs_", reference, "_ref")
result_path <- file.path(outdir, paste0(prefix, ".txt"))
write.table(report, result_path, sep = "\t", quote = FALSE, row.names = FALSE)

bed_path <- file.path(outdir, paste0(prefix, ".with_stats.bed"))
if (nrow(report) && all(c("seqnames", "start", "end") %in% names(report))) {
  stat_columns <- setdiff(names(report), c("seqnames", "start", "end"))
  peak_ids <- vapply(seq_len(nrow(report)), function(i) {
    location <- paste0(report$seqnames[[i]], ":", report$start[[i]], "-", report$end[[i]])
    stats <- vapply(stat_columns, function(column) paste0(column, "=", as.character(report[[column]][[i]])), character(1))
    paste(c(location, stats), collapse = "|")
  }, character(1))
  score <- if ("Fold" %in% names(report)) report$Fold else rep(0, nrow(report))
  bed <- data.frame(report$seqnames, pmax(as.integer(report$start) - 1L, 0L), report$end, peak_ids, score, check.names = FALSE)
  write.table(bed, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  file.create(bed_path)
}

plot_png <- function(filename, expression) {
  png(file.path(outdir, filename), width = 1000, height = 800, res = 150)
  try(force(expression), silent = TRUE)
  dev.off()
}
plot_png("diffbind_CorrHeatmap_byDiffPeaks.png", plot(db, contrast = 1))
plot_png("diffbind_pca_byNormCounts.png", dba.plotPCA(db, DBA_CONDITION, label = DBA_CONDITION))
plot_png("diffbind_pca_byDiffPeaks.png", dba.plotPCA(db, contrast = 1, label = DBA_CONDITION))
plot_png("diffbind_volcano_byDiffPeaks.png", dba.plotVolcano(db))

writeLines(c(
  paste("status", "complete", sep = "\t"),
  paste("reference", reference, sep = "\t"),
  paste("comparison", comparison, sep = "\t"),
  paste("sample_sheet", sample_sheet_path, sep = "\t"),
  paste("result", result_path, sep = "\t"),
  paste("significant_peaks", nrow(report), sep = "\t")
), file.path(outdir, "_COMPLETE"))
