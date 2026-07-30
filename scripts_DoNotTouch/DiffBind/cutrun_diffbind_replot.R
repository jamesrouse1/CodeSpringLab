#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: cutrun_diffbind_replot.R <diffbind_object.rds> <out_dir>")

.libPaths(c("/grid/bsr/data/data/oeldemer/4.3", .libPaths()))
suppressPackageStartupMessages(library(DiffBind))

object_path <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = TRUE)
target <- file.path(out_dir, "pca_differential_peaks.png")
temporary <- paste0(target, ".tmp.png")

db <- readRDS(object_path)
png(temporary, width = 1100, height = 850, res = 140)
result <- try(dba.plotPCA(db, contrast = 1L, label = DBA_ID), silent = TRUE)
if (!inherits(result, "try-error") && inherits(result, "trellis")) print(result)
dev.off()

if (inherits(result, "try-error")) {
  unlink(temporary, force = TRUE)
  stop("DiffBind could not create the differential-peak PCA: ", as.character(result))
}
if (!file.exists(temporary) || file.info(temporary)$size <= 0) stop("DiffBind did not create a non-empty differential-peak PCA.")
if (!file.rename(temporary, target)) stop("Could not move the completed PCA into place: ", target)
writeLines(as.character(Sys.time()), file.path(out_dir, "_PCA_DIFFERENTIAL_COMPLETE"))
