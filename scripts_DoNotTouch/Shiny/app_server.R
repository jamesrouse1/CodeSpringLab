library(shiny)
DT_AVAILABLE <- requireNamespace("DT", quietly = TRUE)

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
config_path <- Sys.getenv("RNASEQ_SHINY_CONFIG", unset = file.path(app_dir, "shiny_results_config.R"))
config_env <- new.env(parent = baseenv())

if (file.exists(config_path)) {
  sys.source(config_path, envir = config_env)
}

config_value <- function(name, default = NULL) {
  if (exists(name, envir = config_env, inherits = FALSE)) {
    get(name, envir = config_env, inherits = FALSE)
  } else {
    default
  }
}

app_dir <- normalizePath(
  path.expand(config_value("app_dir", app_dir)),
  winslash = "/",
  mustWork = FALSE
)
project_name <- config_value("project_name", Sys.getenv("RNASEQ_SHINY_PROJECT", unset = "example_dataset"))
results_root <- normalizePath(
  path.expand(config_value("results_root", Sys.getenv("RNASEQ_SHINY_RESULTS_ROOT", unset = "~/csl_results"))),
  winslash = "/",
  mustWork = FALSE
)
data_dir <- normalizePath(
  path.expand(config_value("data_dir", file.path(results_root, project_name, "data"))),
  winslash = "/",
  mustWork = FALSE
)
fastqc_dir <- file.path(data_dir, "fastqc")
fastqc_cutadapt_dir <- file.path(data_dir, "fastqc_cutadapt")
design_matrix_path <- normalizePath(
  path.expand(config_value("design_matrix_path", file.path(data_dir, "design_matrix", "design_matrix.txt"))),
  winslash = "/",
  mustWork = FALSE
)
star_summary_path <- file.path(data_dir, "star_summary", "summary_matrix.txt")
counts_dir <- file.path(data_dir, "counts")
count_matrix_path <- file.path(counts_dir, "count_matrix.txt")
featurecounts_summary_path <- file.path(counts_dir, "featurecounts_summary.txt")
kallisto_dir <- file.path(data_dir, "kallisto")
rsem_dir <- file.path(data_dir, "rsem")
gtf_dir <- file.path(data_dir, "gtf")
pick_first_existing <- function(filename, dirs) {
  for (d in dirs) {
    p <- file.path(d, filename)
    if (file.exists(p)) return(p)
  }
  NULL
}
mouse_gtf_path <- pick_first_existing("gencode.vM29.annotation.gtf", c(
  gtf_dir,
  "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode"
))
human_gtf_path <- pick_first_existing("gencode.v42.chr_patch_hapl_scaff.annotation.gtf", c(
  gtf_dir,
  "/grid/bsr/data/data/utama/genome/hg38_p13_gencode"
))
deseq2_dir <- file.path(data_dir, "deseq2")
deseq2_gene_name_dir <- file.path(data_dir, "deseq2_gene_name")
gseapy_dir <- file.path(data_dir, "gseapy")
logo_search_dirs <- unique(c(
  path.expand(unlist(config_value("logo_search_dirs", character(0)))),
  path.expand("~/CodeSpringLab/scripts_DoNotTouch"),
  dirname(app_dir),
  app_dir
))

find_first_existing <- function(filename, dirs) {
  pick_first_existing(filename, dirs)
}

logo_csl_path <- find_first_existing("Logo_CSL.png", logo_search_dirs)
logo_path <- find_first_existing("Logo.png", logo_search_dirs)

safe_read_delim <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  args <- list(...)
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE) && is.null(args$row.names)) {
      args$file <- path
      args$data.table <- FALSE
      args$showProgress <- FALSE
      if (is.null(args$check.names)) args$check.names <- FALSE
      return(do.call(data.table::fread, args))
    }
    do.call(read.delim, c(list(file = path), args))
  }, error = function(e) NULL)
}

first_or_null <- function(x) {
  if (length(x)) x[1] else NULL
}

safe_tabular <- function(df, cols) {
  if (is.null(df)) return(as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols)))
  df
}

status_box <- function(message, tone = c("info", "warning", "error")) {
  tone <- match.arg(tone)
  styles <- list(
    info = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
    warning = "margin-bottom: 12px; padding: 10px 12px; background: #fff7e6; border: 1px solid #f0c36d; border-radius: 6px;",
    error = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;"
  )
  tags$div(style = styles[[tone]], message)
}

available_qc_samples <- function(base_dir) {
  if (!dir.exists(base_dir)) return(character(0))
  files <- list.files(base_dir, pattern = "_(fastqc|screen)\\.html$", full.names = FALSE)
  ids <- sub("_(fastqc|screen)\\.html$", "", files)
  ids <- sub("([._-]R)[12]([._-]?[0-9]*)$", "", ids, ignore.case = TRUE)
  sort(unique(ids[ids != files]))
}

build_star_summary_from_outputs <- function(star_dir, out_path) {
  files <- if (dir.exists(star_dir)) list.files(star_dir, pattern = "Log\\.final\\.out$", recursive = TRUE, full.names = TRUE) else character(0)
  if (!length(files)) return(NULL)
  read_one <- function(path) {
    sample <- basename(dirname(path))
    if (!nzchar(sample) || identical(sample, ".")) sample <- sub("Log\\.final\\.out$", "", basename(path))
    x <- tryCatch(read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(x) || ncol(x) < 2 || !nrow(x)) return(NULL)
    metric <- trimws(gsub("\\|$", "", trimws(as.character(x[[1]]))))
    value <- trimws(as.character(x[[2]]))
    data.frame(metric = metric, value = value, sample = sample, stringsAsFactors = FALSE)
  }
  parts <- Filter(Negate(is.null), lapply(files, read_one))
  if (!length(parts)) return(NULL)
  out <- Reduce(function(a, b) merge(a, b, by = "metric", all = TRUE), lapply(parts, function(x) {
    y <- x[, c("metric", "value")]
    names(y)[2] <- unique(x$sample)[1]
    y
  }))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}

build_featurecounts_outputs_from_samples <- function(feature_dir, count_path, summary_path) {
  files <- if (dir.exists(feature_dir)) list.files(feature_dir, pattern = "_counts\\.txt$", recursive = TRUE, full.names = TRUE) else character(0)
  if (!length(files)) return(list(count_matrix = NULL, summary = NULL))
  read_counts <- function(path) {
    x <- tryCatch(read.table(path, sep = "\t", header = TRUE, quote = "\"", comment.char = "#", check.names = FALSE), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) return(NULL)
    gene_col <- intersect(c("Geneid", "gene_id", "gene_name", "GeneID"), names(x))[1]
    if (is.na(gene_col)) gene_col <- names(x)[1]
    count_col <- tail(names(x), 1)
    sample <- basename(dirname(path))
    if (!nzchar(sample) || identical(sample, "featurecounts")) sample <- sub("_counts\\.txt$", "", basename(path))
    data.frame(Geneid = x[[gene_col]], value = suppressWarnings(as.numeric(x[[count_col]])), sample = sample, stringsAsFactors = FALSE)
  }
  count_parts <- Filter(Negate(is.null), lapply(files, read_counts))
  count_matrix <- NULL
  if (length(count_parts)) {
    count_matrix <- Reduce(function(a, b) merge(a, b, by = "Geneid", all = TRUE), lapply(count_parts, function(x) {
      y <- x[, c("Geneid", "value")]
      names(y)[2] <- unique(x$sample)[1]
      y
    }))
    count_matrix[is.na(count_matrix)] <- 0
    dir.create(dirname(count_path), recursive = TRUE, showWarnings = FALSE)
    write.table(count_matrix, count_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  read_summary <- function(path) {
    summary_file <- paste0(path, ".summary")
    sample <- basename(dirname(path))
    if (!nzchar(sample) || identical(sample, "featurecounts")) sample <- sub("_counts\\.txt$", "", basename(path))
    if (!file.exists(summary_file)) {
      counts <- read_counts(path)
      if (is.null(counts)) return(NULL)
      return(data.frame(Status = "Assigned", value = sum(counts$value, na.rm = TRUE), sample = sample, stringsAsFactors = FALSE))
    }
    x <- tryCatch(read.table(summary_file, sep = "\t", header = TRUE, quote = "\"", comment.char = "", check.names = FALSE), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) return(NULL)
    names(x)[1] <- "Status"
    value_col <- setdiff(names(x), "Status")[1]
    data.frame(Status = x$Status, value = suppressWarnings(as.numeric(x[[value_col]])), sample = sample, stringsAsFactors = FALSE)
  }
  summary_parts <- Filter(Negate(is.null), lapply(files, read_summary))
  summary <- NULL
  if (length(summary_parts)) {
    summary <- Reduce(function(a, b) merge(a, b, by = "Status", all = TRUE), lapply(summary_parts, function(x) {
      y <- x[, c("Status", "value")]
      names(y)[2] <- unique(x$sample)[1]
      y
    }))
    summary[is.na(summary)] <- 0
    dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
    write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  list(count_matrix = count_matrix, summary = summary)
}

design_df <- safe_read_delim(design_matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
if (is.null(design_df)) {
  design_df <- data.frame(sample = character(0), filename = character(0), stringsAsFactors = FALSE)
}
if (ncol(design_df) > 0 && !("sample" %in% colnames(design_df))) {
  colnames(design_df)[1] <- "sample"
}
if (ncol(design_df) > 1 && !("filename" %in% colnames(design_df))) {
  colnames(design_df)[ncol(design_df)] <- "filename"
}
if (!("sample" %in% colnames(design_df))) {
  design_df$sample <- character(0)
}
design_df$sample <- trimws(as.character(design_df$sample))
if (!("filename" %in% colnames(design_df))) {
  design_df$filename <- design_df$sample
}
comparison_columns <- setdiff(colnames(design_df), c("sample", "filename"))
sample_sort_choices <- c("Alphabetical" = "__alpha__", setNames(comparison_columns, comparison_columns))

fastq_stem <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) x <- ""
  x <- trimws(as.character(x))
  if (grepl(",", x, fixed = TRUE)) x <- strsplit(x, ",", fixed = TRUE)[[1]][1]
  x <- basename(x)
  x <- sub("_R[12]_001\\.f(ast)?q\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("_R[12]\\.f(ast)?q\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.f(ast)?q\\.gz$", "", x, ignore.case = TRUE)
  x
}

sample_fastq_stems <- setNames(
  vapply(design_df$filename, fastq_stem, character(1)),
  design_df$sample
)

star_raw <- safe_read_delim(
  star_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE
)
if (is.null(star_raw) || !ncol(star_raw)) {
  star_raw <- build_star_summary_from_outputs(file.path(data_dir, "star"), star_summary_path)
}
if (!is.null(star_raw) && ncol(star_raw) > 0) {
  colnames(star_raw)[1] <- "metric"
}
star_metrics_keep <- c(
  "Number of input reads",
  "Uniquely mapped reads %",
  "% of reads mapped to multiple loci",
  "% of reads mapped to too many loci",
  "% of reads unmapped: too short",
  "Mismatch rate per base, %",
  "Average input read length"
)
if (!is.null(star_raw) && "metric" %in% colnames(star_raw)) {
  star_raw$metric <- trimws(gsub("\\|$", "", trimws(star_raw$metric)))
  star_summary_df <- star_raw[star_raw$metric %in% star_metrics_keep, , drop = FALSE]
} else {
  star_summary_df <- data.frame(metric = character(0), stringsAsFactors = FALSE)
}

count_matrix_df <- safe_read_delim(count_matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
featurecounts_rebuilt <- NULL
if (is.null(count_matrix_df) || !ncol(count_matrix_df)) {
  featurecounts_rebuilt <- build_featurecounts_outputs_from_samples(file.path(data_dir, "featurecounts"), count_matrix_path, featurecounts_summary_path)
  count_matrix_df <- featurecounts_rebuilt$count_matrix
}
if (is.null(count_matrix_df)) {
  count_matrix_df <- data.frame(Geneid = character(0), stringsAsFactors = FALSE)
}
if ("Geneid" %in% colnames(count_matrix_df) && ncol(count_matrix_df) > 1) {
  count_matrix_nonzero_df <- count_matrix_df[rowSums(count_matrix_df[, setdiff(colnames(count_matrix_df), "Geneid"), drop = FALSE] != 0) > 0, , drop = FALSE]
  count_matrix_nonzero_df <- count_matrix_nonzero_df[order(count_matrix_nonzero_df$Geneid), , drop = FALSE]
} else {
  count_matrix_nonzero_df <- count_matrix_df
}
featurecounts_raw <- safe_read_delim(
  featurecounts_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE
)
featurecounts_malformed <- !is.null(featurecounts_raw) &&
  all(c("sample", "total_counts") %in% tolower(colnames(featurecounts_raw)))
if (is.null(featurecounts_raw) || !ncol(featurecounts_raw) || featurecounts_malformed) {
  if (is.null(featurecounts_rebuilt)) {
    featurecounts_rebuilt <- build_featurecounts_outputs_from_samples(file.path(data_dir, "featurecounts"), count_matrix_path, featurecounts_summary_path)
  }
  featurecounts_raw <- featurecounts_rebuilt$summary
}
if (!is.null(featurecounts_raw) && ncol(featurecounts_raw) > 0) {
  colnames(featurecounts_raw)[1] <- "Status"
  featurecounts_raw$Status <- trimws(featurecounts_raw$Status)
  featurecounts_summary_df <- featurecounts_raw[
    featurecounts_raw$Status != "Status" &
      rowSums(featurecounts_raw[, setdiff(colnames(featurecounts_raw), "Status"), drop = FALSE] != 0) > 0,
    ,
    drop = FALSE
  ]
} else {
  featurecounts_summary_df <- data.frame(Status = character(0), stringsAsFactors = FALSE)
}

fastqc_raw_available <- dir.exists(fastqc_dir) && length(list.files(fastqc_dir, pattern = "_(fastqc|screen)\\.html$", recursive = FALSE)) > 0
fastqc_trim_available <- dir.exists(fastqc_cutadapt_dir) && length(list.files(fastqc_cutadapt_dir, pattern = "_(fastqc|screen)\\.html$", recursive = FALSE)) > 0
star_available <- nrow(star_summary_df) > 0
count_matrix_available <- nrow(count_matrix_df) > 0
featurecounts_available <- nrow(featurecounts_summary_df) > 0
kallisto_abundance_files <- if (dir.exists(kallisto_dir)) list.files(kallisto_dir, pattern = "abundance.tsv", recursive = TRUE, full.names = TRUE) else character(0)
kallisto_sample_names <- sort(unique(basename(dirname(kallisto_abundance_files))))
kallisto_saved_matrices <- if (dir.exists(counts_dir)) list.files(counts_dir, pattern = "^kallisto_.*_matrix.*\\.txt$", full.names = TRUE) else character(0)
rsem_saved_matrices <- if (dir.exists(counts_dir)) list.files(counts_dir, pattern = "^rsem_.*_matrix.*\\.txt$", full.names = TRUE) else character(0)
kallisto_available <- length(kallisto_sample_names) > 0 || length(kallisto_saved_matrices) > 0
rsem_available <- (dir.exists(rsem_dir) && length(list.files(rsem_dir, pattern = "\\.(genes|isoforms)\\.results$", recursive = TRUE)) > 0) || length(rsem_saved_matrices) > 0
deseq2_output_dirs <- c(deseq2_dir, deseq2_gene_name_dir)
deseq2_available <- any(vapply(
  deseq2_output_dirs,
  function(d) dir.exists(d) && length(list.files(d, pattern = "^normalized_counts_.*\\.txt$", recursive = FALSE)) > 0,
  logical(1)
))
gseapy_results_available <- function() {
  dir.exists(gseapy_dir) && length(list.files(gseapy_dir, pattern = "report\\.gseapy\\..*\\.csv$", recursive = TRUE, full.names = TRUE)) > 0
}
gseapy_available <- gseapy_results_available()

design_samples <- unique(as.character(design_df$sample))
design_samples <- design_samples[!is.na(design_samples) & nzchar(design_samples)]
fallback_sample_sources <- unique(c(
  setdiff(colnames(star_summary_df), "metric"),
  setdiff(colnames(count_matrix_df), "Geneid"),
  available_qc_samples(fastqc_dir),
  available_qc_samples(fastqc_cutadapt_dir),
  if (kallisto_available) kallisto_sample_names else character(0),
  if (rsem_available) list.files(rsem_dir, full.names = FALSE) else character(0)
))
fallback_sample_sources <- fallback_sample_sources[
  !is.na(fallback_sample_sources) &
    nzchar(fallback_sample_sources) &
    !grepl("\\.(txt|tsv|csv|html|png|pdf)$", fallback_sample_sources, ignore.case = TRUE)
]
samples <- if (length(design_samples)) design_samples else sort(unique(fallback_sample_sources))
default_show_trimmed <- if (fastqc_trim_available && !fastqc_raw_available) TRUE else if (fastqc_raw_available && !fastqc_trim_available) FALSE else fastqc_trim_available
featurecounts_sample_choices <- setdiff(colnames(featurecounts_summary_df), "Status")
if (length(samples)) {
  design_ordered_featurecounts <- samples[samples %in% featurecounts_sample_choices]
  other_featurecounts <- setdiff(featurecounts_sample_choices, design_ordered_featurecounts)
  featurecounts_sample_choices <- c(design_ordered_featurecounts, other_featurecounts)
}

if (dir.exists(fastqc_dir)) addResourcePath("fastqc_results", fastqc_dir)
if (dir.exists(fastqc_cutadapt_dir)) addResourcePath("fastqc_cutadapt_results", fastqc_cutadapt_dir)
if (dir.exists(deseq2_dir)) addResourcePath("deseq2_results", deseq2_dir)
if (dir.exists(deseq2_gene_name_dir)) addResourcePath("deseq2_gene_name_results", deseq2_gene_name_dir)
if (dir.exists(gseapy_dir)) addResourcePath("gseapy_results", gseapy_dir)
if (!is.null(logo_csl_path)) addResourcePath("logo_csl_asset", dirname(logo_csl_path))
if (!is.null(logo_path)) addResourcePath("logo_asset", dirname(logo_path))

normalize_sample_token <- function(x) {
  x <- tolower(trimws(as.character(value_or(x, ""))))
  gsub("[^a-z0-9]+", "", x)
}

find_qc_html <- function(base_dir, sample_name, read_name, suffix) {
  if (!dir.exists(base_dir)) return(file.path(base_dir, paste0(sample_name, "_", read_name, "_001_", suffix, ".html")))
  sample_key <- as.character(value_or(sample_name, ""))
  sample_stem <- unname(sample_fastq_stems[sample_key])
  if (is.null(sample_stem) || length(sample_stem) == 0 || is.na(sample_stem)) sample_stem <- sample_name
  candidates <- unique(c(sample_name, sample_stem, gsub("_", "", sample_name), gsub("_", "", sample_stem)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  expected <- file.path(base_dir, paste0(candidates, "_", read_name, "_001_", suffix, ".html"))
  hit <- expected[file.exists(expected)]
  if (length(hit)) return(hit[1])

  files <- list.files(base_dir, pattern = paste0("([._-]", read_name, "([._-]?[0-9]*)?)_", suffix, "\\.html$"), full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(first_or_null(expected))
  file_stems <- sub(paste0("_", suffix, "\\.html$"), "", basename(files), ignore.case = TRUE)
  file_stems <- sub(paste0("([._-]", read_name, "([._-]?[0-9]*)?)$"), "", file_stems, ignore.case = TRUE)
  file_keys <- normalize_sample_token(file_stems)
  cand_keys <- normalize_sample_token(candidates)
  idx <- match(cand_keys, file_keys, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx)) return(files[idx[1]])
  first_or_null(expected)
}

qc_file_map <- function(sample_name, read_name = c("R1", "R2"), trimmed = FALSE) {
  read_name <- match.arg(read_name)
  base_dir <- if (trimmed) fastqc_cutadapt_dir else fastqc_dir
  list(
    fastqc = find_qc_html(base_dir, sample_name, read_name, "fastqc"),
    screen = find_qc_html(base_dir, sample_name, read_name, "screen")
  )
}

iframe_or_message <- function(path, resource_prefix, height = "calc(100vh - 210px)", missing_message = NULL) {
  if (!file.exists(path)) {
    return(tags$p(value_or(missing_message, sprintf("Missing file: %s", basename(path)))))
  }
  rel <- basename(path)
  height_css <- if (is.numeric(height)) paste0(height, "px") else height
  tags$iframe(
    class = "qc-report-frame",
    src = file.path(resource_prefix, rel),
    style = sprintf("width: 100%%; height: %s; min-height: 680px; border: 1px solid #d7e0ea;", height_css)
  )
}

kallisto_abundance_path <- function(sample_name) {
  direct <- file.path(kallisto_dir, sample_name, "abundance.tsv")
  if (file.exists(direct)) return(direct)
  sample_key <- normalize_sample_token(sample_name)
  idx <- match(sample_key, normalize_sample_token(kallisto_sample_names), nomatch = 0)
  if (idx > 0) {
    candidate <- file.path(kallisto_dir, kallisto_sample_names[idx], "abundance.tsv")
    if (file.exists(candidate)) return(candidate)
  }
  direct
}

read_kallisto_abundance <- function(sample_name) {
  file_path <- kallisto_abundance_path(sample_name)
  if (!file.exists(file_path)) {
    return(NULL)
  }
  df <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  parts <- strsplit(df$target_id, "\\|", fixed = FALSE)
  df$transcript_id <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  df$gene_id <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  df$transcript_name <- vapply(parts, function(x) if (length(x) >= 5) x[5] else NA_character_, character(1))
  df$gene_symbol <- vapply(parts, function(x) if (length(x) >= 6) x[6] else NA_character_, character(1))
  df$biotype <- vapply(parts, function(x) if (length(x) >= 8) x[8] else NA_character_, character(1))
  df <- df[order(-df$tpm, -df$est_counts, df$gene_symbol, df$transcript_name), , drop = FALSE]
  df[, c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype", "tpm", "est_counts", "length", "eff_length")]
}

kallisto_matrix_path <- file.path(kallisto_dir, "kallisto_transcript_matrix.tsv")

read_kallisto_matrix <- function() {
  if (!file.exists(kallisto_matrix_path)) return(NULL)
  safe_read_delim(kallisto_matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_kallisto_matrix <- function(df) {
  if (is.null(df)) return(invisible(FALSE))
  write.table(df, kallisto_matrix_path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}

build_kallisto_matrix <- function(sample_names) {
  per_sample <- lapply(sample_names, function(s) {
    df <- read_kallisto_abundance(s)
    if (is.null(df)) return(NULL)
    keep <- df[, c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype", "tpm"), drop = FALSE]
    colnames(keep)[colnames(keep) == "tpm"] <- s
    keep
  })
  per_sample <- per_sample[!vapply(per_sample, is.null, logical(1))]
  if (!length(per_sample)) return(NULL)
  merged <- Reduce(
    function(x, y) merge(
      x, y,
      by = c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype"),
      all = TRUE
    ),
    per_sample
  )
  sample_cols <- setdiff(colnames(merged), c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype"))
  merged[sample_cols][is.na(merged[sample_cols])] <- 0
  merged
}

read_rsem_result <- function(sample_name, data_type = c("genes", "isoforms")) {
  data_type <- match.arg(data_type)
  suffix <- if (identical(data_type, "genes")) "genes.results" else "isoforms.results"
  file_path <- file.path(rsem_dir, sample_name, paste0(sample_name, ".", suffix))
  if (!file.exists(file_path)) {
    return(NULL)
  }
  read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
}

build_rsem_matrix <- function(sample_names, data_type = c("genes", "isoforms"), metric = "expected_count") {
  data_type <- match.arg(data_type)
  key_cols <- if (identical(data_type, "genes")) c("gene_id") else c("transcript_id", "gene_id")
  per_sample <- lapply(sample_names, function(s) {
    df <- read_rsem_result(s, data_type)
    if (is.null(df) || !(metric %in% colnames(df))) return(NULL)
    keep <- df[, c(key_cols, metric), drop = FALSE]
    colnames(keep)[colnames(keep) == metric] <- s
    keep
  })
  per_sample <- per_sample[!vapply(per_sample, is.null, logical(1))]
  if (!length(per_sample)) return(NULL)
  merged <- Reduce(
    function(x, y) merge(x, y, by = key_cols, all = TRUE),
    per_sample
  )
  sample_cols <- setdiff(colnames(merged), key_cols)
  merged[sample_cols][is.na(merged[sample_cols])] <- 0
  id_col <- if (identical(data_type, "genes")) "gene_id" else "transcript_id"
  merged[order(merged[[id_col]]), , drop = FALSE]
}

rsem_metric_token <- function(metric) {
  switch(
    as.character(metric),
    "TPM" = "tpm",
    "FPKM" = "fpkm",
    "expected_count" = "expected_count",
    tolower(as.character(metric))
  )
}

rsem_saved_matrix_path <- function(data_type = "genes", metric = "TPM", label_mode = "gene_id") {
  token <- rsem_metric_token(metric)
  stem <- if (identical(data_type, "isoforms")) {
    paste0("rsem_isoform_", token, "_matrix")
  } else {
    paste0("rsem_", token, "_matrix")
  }
  converted <- file.path(counts_dir, paste0(stem, "_gene_name.txt"))
  base <- file.path(counts_dir, paste0(stem, ".txt"))
  if (identical(label_mode, "gene_name") && file.exists(converted)) return(converted)
  if (file.exists(base)) return(base)
  if (identical(label_mode, "gene_name") && file.exists(converted)) return(converted)
  base
}

kallisto_metric_token <- function(metric) {
  switch(
    as.character(metric),
    "tpm" = "tpm",
    "est_counts" = "est_counts",
    "TPM" = "tpm",
    "Estimated counts" = "est_counts",
    tolower(as.character(metric))
  )
}

kallisto_saved_matrix_path <- function(metric = "tpm", label_mode = "target_id") {
  token <- kallisto_metric_token(metric)
  stem <- paste0("kallisto_", token, "_matrix")
  converted <- file.path(counts_dir, paste0(stem, "_gene_name.txt"))
  base <- file.path(counts_dir, paste0(stem, ".txt"))
  if (identical(label_mode, "gene_name") && file.exists(converted)) return(converted)
  if (file.exists(base)) return(base)
  legacy <- file.path(kallisto_dir, "kallisto_transcript_matrix.tsv")
  if (identical(token, "tpm") && file.exists(legacy)) return(legacy)
  base
}

read_saved_count_matrix <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  safe_read_delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

gene_name_matrix_path <- function(path) {
  sub("(\\.txt|\\.tsv)$", "_gene_name.txt", path)
}

looks_like_transcript_id <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(FALSE)
  ids <- head(ids, 500)
  mean(grepl("^(ENST|ENSMUST|ENS[A-Z]*T)[0-9]+", ids)) >= 0.5
}

parse_kallisto_target_metadata <- function(target_id) {
  target_id <- as.character(target_id)
  parts <- strsplit(target_id, "\\|", fixed = FALSE)
  data.frame(
    target_id = target_id,
    transcript_id = vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)),
    gene_id = vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)),
    transcript_name = vapply(parts, function(x) if (length(x) >= 5) x[5] else NA_character_, character(1)),
    gene_name = vapply(parts, function(x) if (length(x) >= 6) x[6] else NA_character_, character(1)),
    biotype = vapply(parts, function(x) if (length(x) >= 8) x[8] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

collect_kallisto_filter_values <- function(sample_names, column_name) {
  vals <- unlist(lapply(sample_names, function(s) {
    df <- read_kallisto_abundance(s)
    if (is.null(df) || !(column_name %in% colnames(df))) return(character(0))
    as.character(df[[column_name]])
  }), use.names = FALSE)
  vals <- sort(unique(vals))
  vals[!is.na(vals) & nzchar(vals)]
}

strip_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

detect_species_from_ids <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  x_head <- head(x, 500)
  if (sum(grepl("^ENSMUSG", x_head) | grepl("^ENSMUST", x_head)) > 0) return("mouse")
  if (sum(grepl("^ENSG", x_head) | grepl("^ENST", x_head)) > 0) return("human")
  if (sum(grepl("^Gm[0-9]+$", x_head) | grepl("^[A-Z][a-z0-9]+$", x_head) | grepl("Rik$", x_head)) > sum(grepl("^[A-Z0-9\\-]+$", x_head))) return("mouse")
  if (sum(grepl("^[A-Z0-9\\-]+$", x_head)) > 0) return("human")
  NA_character_
}

read_gtf_gene_map <- function(gtf_path) {
  if (!file.exists(gtf_path)) return(NULL)
  lines <- readLines(gtf_path, warn = FALSE)
  lines <- lines[substr(lines, 1, 1) != "#"]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1)) >= 9]
  gene_fields <- fields[vapply(fields, function(x) identical(x[3], "gene"), logical(1))]
  if (!length(gene_fields)) return(NULL)
  attrs <- vapply(gene_fields, function(x) x[9], character(1))
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attrs)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attrs)
  keep <- gene_id != attrs & gene_name != attrs & nzchar(gene_id) & nzchar(gene_name)
  map <- data.frame(
    gene_id = gene_id[keep],
    gene_id_stripped = strip_version(gene_id[keep]),
    gene_name = gene_name[keep],
    stringsAsFactors = FALSE
  )
  map <- map[!duplicated(map$gene_id), , drop = FALSE]
  map
}

looks_like_gene_id <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(FALSE)
  ids <- head(ids, 500)
  mean(grepl("^(ENSG|ENSMUSG|ENS[A-Z]*G)[0-9]+", ids)) >= 0.5
}

convert_gene_labels <- function(ids, map_df) {
  ids <- as.character(ids)
  current_species <- detect_species_from_ids(ids)
  current_type <- if (looks_like_gene_id(ids)) "gene_id" else "gene_name"
  target_type <- if (identical(current_type, "gene_id")) "gene_name" else "gene_id"
  if (is.null(map_df) || !length(ids)) {
    return(list(values = ids, from = current_type, to = target_type, species = current_species, changed = FALSE, mapped = 0))
  }
  if (identical(current_type, "gene_id")) {
    idx <- match(strip_version(ids), map_df$gene_id_stripped)
    converted <- map_df$gene_name[idx]
    mapped <- sum(!is.na(converted) & nzchar(converted))
    converted[is.na(converted) | !nzchar(converted)] <- ids[is.na(converted) | !nzchar(converted)]
    return(list(values = converted, from = "gene_id", to = "gene_name", species = current_species, changed = mapped > 0, mapped = mapped))
  }
  name_map <- map_df[!duplicated(map_df$gene_name), c("gene_name", "gene_id"), drop = FALSE]
  idx <- match(ids, name_map$gene_name)
  converted <- name_map$gene_id[idx]
  mapped <- sum(!is.na(converted) & nzchar(converted))
  converted[is.na(converted) | !nzchar(converted)] <- ids[is.na(converted) | !nzchar(converted)]
  list(values = converted, from = "gene_name", to = "gene_id", species = current_species, changed = mapped > 0, mapped = mapped)
}

value_or <- function(x, default) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) default else x
}

pvalue_columns <- function(df) {
  if (is.null(df) || !NCOL(df)) return(character(0))
  names(df)[grepl("(^p$|pvalue|p\\.value|p_value|p-val|p\\.val|padj|adj\\.p|fdr|qvalue|q\\.value|q_value|q-val)", names(df), ignore.case = TRUE)]
}

pvalue_render_js <- function(digits = 3) {
  DT::JS(sprintf(
    "function(data, type, row, meta) {
       if (type === 'display' || type === 'filter') {
         var x = parseFloat(data);
         if (!isNaN(x) && isFinite(x)) return x.toExponential(%d);
       }
       return data;
     }",
    digits
  ))
}

table_widget <- function(id) {
  tags$div(
    class = "table-scroll-shell",
    if (DT_AVAILABLE) {
      DT::DTOutput(id)
    } else {
      tableOutput(id)
    }
  )
}

simple_dt <- function(df, page_length = 50, scroll_y = "520px", dom = "lfrtip", container = NULL) {
  pvalue_cols_all <- pvalue_columns(df)
  pvalue_targets <- match(pvalue_cols_all, names(df), nomatch = 0) - 1
  pvalue_targets <- pvalue_targets[pvalue_targets >= 0]
  column_defs <- c(
    list(list(width = "118px", targets = "_all")),
    lapply(pvalue_targets, function(target) list(targets = target, render = pvalue_render_js(3)))
  )
  dt_args <- list(
    data = df,
    rownames = FALSE,
    class = "compact stripe hover cell-border",
    options = list(
      scrollX = TRUE,
      scrollY = scroll_y,
      pageLength = page_length,
      lengthMenu = list(c(25, 50, 100, -1), c("25", "50", "100", "All")),
      paging = TRUE,
      pagingType = "full_numbers",
      dom = dom,
      deferRender = FALSE,
      processing = FALSE,
      searchDelay = 350,
      autoWidth = FALSE,
      columnDefs = column_defs
    )
  )
  if (!is.null(container)) dt_args$container <- container
  dt <- do.call(DT::datatable, dt_args)
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(numeric_cols)) {
    pvalue_cols <- intersect(numeric_cols, pvalue_cols_all)
    integer_cols <- numeric_cols[vapply(df[numeric_cols], function(x) {
      finite <- x[is.finite(x) & !is.na(x)]
      length(finite) == 0 || all(abs(finite - round(finite)) < 1e-8)
    }, logical(1))]
    decimal_cols <- setdiff(numeric_cols, c(integer_cols, pvalue_cols))
    if (length(integer_cols)) dt <- DT::formatRound(dt, columns = integer_cols, digits = 0)
    if (length(decimal_cols)) dt <- DT::formatRound(dt, columns = decimal_cols, digits = 2)
  }
  dt
}

sample_annotation_col <- function(preferred = NULL) {
  preferred <- value_or(preferred, "")
  if (nzchar(preferred) && preferred %in% comparison_columns) return(preferred)
  if ("treatment" %in% comparison_columns) return("treatment")
  first_or_null(comparison_columns)
}

matrix_id_cols <- function(df) {
  sample_cols <- intersect(samples, colnames(df))
  setdiff(colnames(df), sample_cols)
}

format_numeric_commas <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  out <- df
  p_cols <- pvalue_columns(out)
  for (nm in colnames(out)) {
    if (is.numeric(out[[nm]])) {
      if (nm %in% p_cols) {
        out[[nm]] <- ifelse(is.na(out[[nm]]), NA_character_, sprintf("%.3e", out[[nm]]))
      } else {
        out[[nm]] <- ifelse(
          is.na(out[[nm]]),
          NA_character_,
          prettyNum(out[[nm]], big.mark = ",", scientific = FALSE, trim = TRUE)
        )
      }
    }
  }
  out
}

sort_df_by_col <- function(df, sort_col, sort_dir = "asc") {
  if (is.null(df) || !nrow(df) || is.null(sort_col) || !nzchar(sort_col) || !(sort_col %in% colnames(df))) {
    return(df)
  }
  vals <- df[[sort_col]]
  suppress_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(vals), fixed = TRUE)))
  ord <- if (sum(!is.na(suppress_num)) > 0) {
    if (identical(sort_dir, "desc")) order(-suppress_num, na.last = TRUE) else order(suppress_num, na.last = TRUE)
  } else {
    if (identical(sort_dir, "desc")) order(as.character(vals), decreasing = TRUE, na.last = TRUE) else order(as.character(vals), na.last = TRUE)
  }
  df[ord, , drop = FALSE]
}

order_sample_columns <- function(df, annotation_col = "__alpha__", design_df = NULL, id_cols = character(0)) {
  if (is.null(df) || !nrow(df)) return(df)
  id_cols <- intersect(id_cols, colnames(df))
  sample_cols <- setdiff(colnames(df), id_cols)
  if (!length(sample_cols)) return(df)
  if (!is.null(design_df) && !is.null(annotation_col) && nzchar(annotation_col) && !identical(annotation_col, "__alpha__") && annotation_col %in% colnames(design_df)) {
    ann_vals <- as.character(design_df[[annotation_col]][match(sample_cols, design_df$sample)])
    ann_vals[is.na(ann_vals) | !nzchar(ann_vals)] <- "~"
    ord <- order(ann_vals, sample_cols)
  } else {
    ord <- order(sample_cols)
  }
  df[, c(id_cols, sample_cols[ord]), drop = FALSE]
}

build_group_header_container <- function(df, id_cols = character(0), annotation_col = NULL, design_df = NULL) {
  if (!DT_AVAILABLE || is.null(annotation_col) || !nzchar(annotation_col) || is.null(design_df) || !(annotation_col %in% colnames(design_df))) {
    return(NULL)
  }
  cols <- colnames(df)
  sample_cols <- setdiff(cols, id_cols)
  ann_vals <- as.character(design_df[[annotation_col]][match(sample_cols, design_df$sample)])
  ann_vals[is.na(ann_vals) | !nzchar(ann_vals)] <- "Unassigned"
  run_ids <- cumsum(c(TRUE, ann_vals[-1] != ann_vals[-length(ann_vals)]))
  blocks <- split(seq_along(sample_cols), run_ids)
  palette <- c("#0f62c6", "#15936f", "#2a7fbe", "#1f9d8a", "#3c78b4", "#2b8f7d", "#4c6fb0", "#2f7f75")
  top_cells <- list()
  if (length(id_cols)) {
    top_cells[[length(top_cells) + 1]] <- tags$th(colspan = length(id_cols), "")
  }
  for (i in seq_along(blocks)) {
    idx <- blocks[[i]]
    grp <- ann_vals[idx[1]]
    bg <- palette[((i - 1) %% length(palette)) + 1]
    top_cells[[length(top_cells) + 1]] <- tags$th(
      colspan = length(idx),
      style = paste0("text-align:center;color:#fff;background:", bg, ";font-size:11px;white-space:nowrap;"),
      grp
    )
  }
  bottom_cells <- lapply(cols, function(x) tags$th(x))
  tags$table(
    class = "display compact cell-border stripe hover",
    tags$thead(
      tags$tr(top_cells),
      tags$tr(bottom_cells)
    )
  )
}

aggregate_display_matrix <- function(df, label_col, keep_cols = character(0)) {
  sample_cols <- setdiff(colnames(df), c(label_col, keep_cols))
  if (!length(sample_cols)) return(df)
  agg <- aggregate(df[, sample_cols, drop = FALSE], by = list(df[[label_col]]), FUN = sum, na.rm = TRUE)
  colnames(agg)[1] <- label_col
  if (length(keep_cols)) {
    meta <- df[, c(label_col, keep_cols), drop = FALSE]
    meta <- meta[!duplicated(meta[[label_col]]), , drop = FALSE]
    agg <- merge(meta, agg, by = label_col, all.y = TRUE, sort = FALSE)
  }
  agg[order(agg[[label_col]]), , drop = FALSE]
}

sanitize_filename <- function(x) {
  x <- trimws(value_or(x, ""))
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


download_filename <- function(prefix, ext = "csv") {
  clean_project <- sanitize_filename(project_name)
  clean_prefix <- sanitize_filename(prefix)
  if (!nzchar(clean_project)) clean_project <- "rna_seq_project"
  if (!nzchar(clean_prefix)) clean_prefix <- "table"
  paste0(clean_project, "_", clean_prefix, "_", format(Sys.Date(), "%Y%m%d"), ".", ext)
}

safe_download_df <- function(expr, empty_message = "No data available.") {
  tryCatch({
    df <- force(expr)
    if (is.null(df)) {
      data.frame(message = empty_message, stringsAsFactors = FALSE)
    } else {
      as.data.frame(df, stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    data.frame(message = empty_message, stringsAsFactors = FALSE)
  })
}

write_csv_download <- function(df, file) {
  if (is.null(df)) df <- data.frame()
  write.csv(as.data.frame(df, stringsAsFactors = FALSE), file, row.names = FALSE, na = "")
}

download_table_button <- function(id, label) {
  tags$div(
    style = "margin-top: 8px; margin-bottom: 12px;",
    downloadButton(id, label, class = "btn-sm")
  )
}

gseapy_comparison_dir <- function(treatment_value, control_value) {
  if (!dir.exists(gseapy_dir)) return(NULL)
  candidates <- unique(c(
    sprintf("%s_vs_%s", treatment_value, control_value),
    sprintf("%s_vs_%s", sanitize_filename(treatment_value), sanitize_filename(control_value)),
    sprintf("%s_vs_%s", control_value, treatment_value),
    sprintf("%s_vs_%s", sanitize_filename(control_value), sanitize_filename(treatment_value))
  ))
  paths <- file.path(gseapy_dir, candidates)
  hits <- paths[dir.exists(paths) & vapply(paths, function(p) length(list.files(p, pattern = "^report\\.gseapy\\..*\\.csv$", recursive = FALSE)) > 0, logical(1))]
  if (length(hits)) return(hits[1])
  dirs <- list.dirs(gseapy_dir, recursive = FALSE, full.names = TRUE)
  dirs_with_reports <- dirs[vapply(dirs, function(p) length(list.files(p, pattern = "^report\\.gseapy\\..*\\.csv$", recursive = FALSE)) > 0, logical(1))]
  if (length(dirs_with_reports) == 1) return(dirs_with_reports[1])
  paths[1]
}

gseapy_comparison_rel <- function(treatment_value, control_value) {
  comp_dir <- gseapy_comparison_dir(treatment_value, control_value)
  if (is.null(comp_dir)) return(sprintf("%s_vs_%s", treatment_value, control_value))
  root <- normalizePath(gseapy_dir, winslash = "/", mustWork = FALSE)
  full <- normalizePath(comp_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(full, prefix)) substring(full, nchar(prefix) + 1) else basename(full)
}

list_gseapy_collections <- function(comp_dir) {
  if (is.null(comp_dir) || !dir.exists(comp_dir)) return(character(0))
  files <- list.files(comp_dir, pattern = "^report\\.gseapy\\..*\\.csv$", full.names = FALSE)
  sub("^report\\.gseapy\\.(.*)\\.csv$", "\\1", files)
}

read_gseapy_report <- function(comp_dir, collection_name) {
  path <- file.path(comp_dir, sprintf("report.gseapy.%s.csv", collection_name))
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

format_gsea_table <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  out <- df
  p_cols <- intersect(c("NOM p-val", "FDR q-val", "FWER p-val"), colnames(out))
  for (nm in p_cols) {
    vals <- suppressWarnings(as.numeric(out[[nm]]))
    out[[nm]] <- ifelse(
      is.na(vals),
      NA_character_,
      ifelse(vals < 1e-3, formatC(vals, format = "e", digits = 3), formatC(vals, format = "f", digits = 4))
    )
  }
  out
}

pathway_pdf_relpath <- function(treatment_value, control_value, pathway_name, heatmap = FALSE) {
  stems <- unique(c(pathway_name, sanitize_filename(pathway_name)))
  suffix <- if (heatmap) ".heatmap.pdf" else ".pdf"
  for (stem in stems[nzchar(stems)]) {
    rel <- file.path(gseapy_comparison_rel(treatment_value, control_value), "gsea", paste0(stem, suffix))
    full <- file.path(gseapy_dir, rel)
    if (file.exists(full)) return(file.path("gseapy_results", rel))
  }
  NULL
}

pdf_render_relpath <- function(resource_rel) {
  if (is.null(resource_rel)) return(NULL)
  rel <- sub("^gseapy_results/", "", resource_rel)
  base_rel <- tools::file_path_sans_ext(rel)
  for (ext in c(".png", ".jpg", ".jpeg", ".svg")) {
    existing_rel <- paste0(base_rel, ext)
    if (file.exists(file.path(gseapy_dir, existing_rel))) return(file.path("gseapy_results", existing_rel))
  }

  full_pdf <- file.path(gseapy_dir, rel)
  if (!file.exists(full_pdf)) return(NULL)
  png_rel <- paste0(base_rel, ".rendered.png")
  full_png <- file.path(gseapy_dir, png_rel)
  if (!file.exists(full_png)) {
    dir.create(dirname(full_png), recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("pdftools", quietly = TRUE)) {
      tryCatch(
        pdftools::pdf_convert(full_pdf, format = "png", pages = 1, dpi = 300, filenames = full_png),
        error = function(e) NULL
      )
    }
  }
  if (!file.exists(full_png) && requireNamespace("magick", quietly = TRUE)) {
    tryCatch({
      img <- magick::image_read_pdf(full_pdf, density = 300)
      magick::image_write(img[1], path = full_png, format = "png")
    }, error = function(e) NULL)
  }
  if (!file.exists(full_png) && nzchar(Sys.which("pdftoppm"))) {
    tryCatch(
      system2("pdftoppm", c("-png", "-r", "300", "-singlefile", full_pdf, tools::file_path_sans_ext(full_png))),
      error = function(e) NULL
    )
  }
  if (!file.exists(full_png) && nzchar(Sys.which("python3"))) {
    py <- paste(
      "import sys",
      "pdf, out = sys.argv[1], sys.argv[2]",
      "def done():\n    sys.exit(0)",
      "try:\n    import fitz\n    doc = fitz.open(pdf)\n    page = doc.load_page(0)\n    pix = page.get_pixmap(matrix=fitz.Matrix(3, 3), alpha=False)\n    pix.save(out)\n    done()\nexcept Exception:\n    pass",
      "try:\n    import pypdfium2 as pdfium\n    doc = pdfium.PdfDocument(pdf)\n    page = doc[0]\n    bitmap = page.render(scale=3)\n    bitmap.to_pil().save(out)\n    done()\nexcept Exception:\n    pass",
      "try:\n    from pdf2image import convert_from_path\n    imgs = convert_from_path(pdf, dpi=300, first_page=1, last_page=1)\n    imgs[0].save(out)\n    done()\nexcept Exception:\n    pass",
      "sys.exit(1)",
      sep = "\n"
    )
    tryCatch(system2("python3", c("-c", py, full_pdf, full_png), stdout = FALSE, stderr = FALSE), error = function(e) NULL)
  }
  if (file.exists(full_png)) file.path("gseapy_results", png_rel) else NULL
}

js_string <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("'", "\\'", x, fixed = TRUE)
  paste0("'", x, "'")
}

pdf_canvas_block <- function(title, rel) {
  canvas_id <- paste0(
    "gsea_pdf_canvas_",
    gsub("[^A-Za-z0-9_]", "_", substr(paste(title, rel, sep = "_"), 1, 80)),
    "_",
    sample.int(.Machine$integer.max, 1)
  )
  status_id <- paste0(canvas_id, "_status")
  render_js <- sprintf(
    "(function() {
      var url = %s;
      var canvasId = %s;
      var statusId = %s;
      function setStatus(message) {
        var el = document.getElementById(statusId);
        if (el) el.textContent = message || '';
      }
      function renderPdf() {
        var canvas = document.getElementById(canvasId);
        if (!canvas) return;
        if (!window.pdfjsLib) {
          setTimeout(renderPdf, 200);
          return;
        }
        window.pdfjsLib.GlobalWorkerOptions.workerSrc = 'https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.worker.min.js';
        window.pdfjsLib.getDocument(url).promise.then(function(pdf) {
          return pdf.getPage(1);
        }).then(function(page) {
          var container = canvas.parentElement;
          var viewport = page.getViewport({ scale: 1 });
          var targetWidth = Math.max(700, container ? container.clientWidth : 900);
          var scale = targetWidth / viewport.width;
          var scaled = page.getViewport({ scale: scale });
          var ratio = window.devicePixelRatio || 1;
          canvas.width = Math.floor(scaled.width * ratio);
          canvas.height = Math.floor(scaled.height * ratio);
          canvas.style.width = Math.floor(scaled.width) + 'px';
          canvas.style.height = Math.floor(scaled.height) + 'px';
          var ctx = canvas.getContext('2d');
          ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
          setStatus('');
          return page.render({ canvasContext: ctx, viewport: scaled }).promise;
        }).catch(function(err) {
          setStatus('Could not render this PDF in the browser. Open the original PDF below.');
        });
      }
      renderPdf();
    })();",
    js_string(rel),
    js_string(canvas_id),
    js_string(status_id)
  )
  tags$div(
    tags$h5(title),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.min.js"),
    tags$div(
      style = "width: 100%; overflow-x: auto; border: 1px solid #ddd; background: #fff; margin-bottom: 12px;",
      tags$canvas(id = canvas_id, style = "display: block; max-width: 100%;")
    ),
    tags$p(id = status_id, style = "color: #9a3412;"),
    tags$script(HTML(render_js)),
    tags$p(tags$a(href = rel, target = "_blank", "Open original PDF"))
  )
}

deseq2_result_dir <- function(label_mode = "gene_id") {
  if (identical(value_or(label_mode, "gene_id"), "gene_name")) deseq2_gene_name_dir else deseq2_dir
}

normalized_counts_path <- function(treatment_value, control_value, label_mode = "gene_id") {
  file.path(deseq2_result_dir(label_mode), sprintf("normalized_counts_%s_vs_%s(ref).txt", treatment_value, control_value))
}

deg_table_path <- function(treatment_value, control_value, label_mode = "gene_id") {
  file.path(deseq2_result_dir(label_mode), sprintf("DEG_%s_vs_%s(ref).txt", treatment_value, control_value))
}

pca_plot_path <- function(compare_col, treatment_value, control_value) {
  file.path(deseq2_dir, sprintf("pca_%s_%s_vs_%s(ref).png", compare_col, treatment_value, control_value))
}

pca_all_plot_path <- function(compare_col) {
  file.path(deseq2_dir, sprintf("pca_all_samples_%s.png", compare_col))
}

read_normalized_counts <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  first_col <- colnames(df)[1]
  colnames(df)[1] <- "gene_label"
  if ("DESCRIPTION" %in% colnames(df)) {
    df <- df[, c("gene_label", setdiff(colnames(df), c("gene_label", "DESCRIPTION"))), drop = FALSE]
  }
  df[order(df$gene_label), , drop = FALSE]
}

first_normalized_counts_path <- function(label_mode = "gene_id") {
  search_dir <- deseq2_result_dir(label_mode)
  files <- if (dir.exists(search_dir)) {
    list.files(search_dir, pattern = "^normalized_counts_.*\\.txt$", full.names = TRUE)
  } else character(0)
  first_or_null(sort(files))
}

normalized_counts_gene_name_path <- function(treatment_value, control_value) {
  sub("(\\.txt|\\.tsv)$", "_aggregated.txt", gene_name_matrix_path(normalized_counts_path(treatment_value, control_value, "gene_id")))
}

available_deseq2_label_modes <- function(treatment_value = NULL, control_value = NULL) {
  choices <- c("Gene ID" = "gene_id", "Gene name" = "gene_name")
  keep <- vapply(unname(choices), function(mode) {
    if (!is.null(treatment_value) && !is.null(control_value) && nzchar(treatment_value) && nzchar(control_value)) {
      return(file.exists(deg_table_path(treatment_value, control_value, mode)) ||
        file.exists(normalized_counts_path(treatment_value, control_value, mode)))
    }
    search_dir <- deseq2_result_dir(mode)
    dir.exists(search_dir) && length(list.files(search_dir, pattern = "^(DEG|normalized_counts)_.*\\.txt$", full.names = TRUE)) > 0
  }, logical(1))
  choices[keep]
}

expression_matrix_from_df <- function(df, sample_ids = samples) {
  if (is.null(df) || !nrow(df)) return(NULL)
  gene_col <- intersect(c("gene_label", "Geneid", "gene_id", "target_id"), colnames(df))[1]
  if (is.na(gene_col)) gene_col <- colnames(df)[1]
  sample_cols <- intersect(sample_ids, colnames(df))
  if (length(sample_cols) < 2) return(NULL)
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(as.character(df[[gene_col]]))
  mat <- mat[rowSums(is.finite(mat), na.rm = TRUE) == ncol(mat), , drop = FALSE]
  mat <- mat[apply(mat, 1, stats::var, na.rm = TRUE) > 0, , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)
  mat
}

draw_pca_plot <- function(mat, metadata, color_col, include_values = character(0), label_samples = TRUE) {
  validate(need(!is.null(mat) && ncol(mat) >= 2 && nrow(mat) >= 2, "PCA needs at least two samples and two variable genes/transcripts."))
  metadata <- metadata[match(colnames(mat), metadata$sample), , drop = FALSE]
  keep <- !is.na(metadata$sample) & nzchar(metadata$sample)
  color_values <- rep("Samples", nrow(metadata))
  if (nzchar(color_col) && color_col %in% colnames(metadata)) {
    color_values <- as.character(metadata[[color_col]])
    color_values[is.na(color_values) | !nzchar(color_values)] <- "NA"
    selected <- value_or(include_values, character(0))
    selected <- selected[nzchar(selected)]
    if (length(selected)) keep <- keep & color_values %in% selected
  } else {
    color_col <- "Samples"
  }
  mat <- mat[, keep, drop = FALSE]
  color_values <- color_values[keep]
  metadata <- metadata[keep, , drop = FALSE]
  validate(need(ncol(mat) >= 2, "Select at least two samples for PCA."))
  log_mat <- log2(mat + 1)
  pca <- stats::prcomp(t(log_mat), center = TRUE, scale. = FALSE)
  var_pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  plot_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    color_group = factor(color_values),
    stringsAsFactors = FALSE
  )
  title <- "PCA of all selected samples"
  subtitle <- paste("Colored by", color_col)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(PC1, PC2, color = color_group)) +
      ggplot2::geom_point(size = 4.2, alpha = 0.9) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = sprintf("PC1 (%.1f%%)", var_pct[1]),
        y = sprintf("PC2 (%.1f%%)", var_pct[2]),
        color = color_col
      ) +
      ggplot2::theme_classic(base_family = "sans") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 18, color = "#17202f"),
        plot.subtitle = ggplot2::element_text(size = 12, color = "#657084"),
        axis.title = ggplot2::element_text(face = "bold", color = "#17202f"),
        axis.text = ggplot2::element_text(color = "#334155"),
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        panel.grid.major = ggplot2::element_line(color = "#e6edf5", size = 0.3),
        panel.grid.minor = ggplot2::element_blank()
      )
    if (isTRUE(label_samples)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3.2, color = "#17202f", min.segment.length = 0, box.padding = 0.35, point.padding = 0.25, seed = 8)
      } else {
        p <- p + ggplot2::geom_text(ggplot2::aes(label = sample), vjust = -0.8, size = 3)
      }
    }
    print(p)
  } else {
    plot(plot_df$PC1, plot_df$PC2, pch = 19, col = as.integer(plot_df$color_group), xlab = sprintf("PC1 (%.1f%%)", var_pct[1]), ylab = sprintf("PC2 (%.1f%%)", var_pct[2]), main = title)
    if (isTRUE(label_samples)) text(plot_df$PC1, plot_df$PC2, labels = plot_df$sample, pos = 3, cex = 0.8)
    legend("topright", legend = levels(plot_df$color_group), col = seq_along(levels(plot_df$color_group)), pch = 19, bty = "n")
  }
}

read_deg_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
  df <- cbind(gene_label = rownames(df), df, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df
}

sort_deg_table <- function(df, sort_mode = "original", p_col = "padj") {
  if (is.null(df) || !nrow(df) || identical(sort_mode, "original")) return(df)
  if (!("log2FoldChange" %in% colnames(df)) || !(p_col %in% colnames(df))) return(df)
  lfc <- suppressWarnings(as.numeric(df$log2FoldChange))
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  lfc[is.na(lfc)] <- 0
  pvals[is.na(pvals)] <- Inf
  direction_priority <- if (identical(sort_mode, "down")) ifelse(lfc < 0, 0, 1) else ifelse(lfc > 0, 0, 1)
  df[order(direction_priority, pvals, -abs(lfc)), , drop = FALSE]
}

format_deg_table <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  out <- df
  for (nm in c("pvalue", "padj")) {
    if (nm %in% colnames(out)) {
      vals <- out[[nm]]
      out[[nm]] <- ifelse(
        is.na(vals),
        NA_character_,
        ifelse(vals < 1e-3, formatC(vals, format = "e", digits = 3), formatC(vals, format = "f", digits = 4))
      )
    }
  }
  out
}

heatmap_theme_options <- function(style) {
  style <- value_or(style, "publication")
  if (identical(style, "minimal")) return(list(palette = "blue_red", border = "none"))
  if (identical(style, "compact")) return(list(palette = "gray_red", border = "none"))
  if (identical(style, "labeled")) return(list(palette = "blue_red", border = "light"))
  if (identical(style, "contrast")) return(list(palette = "green_purple", border = "light"))
  list(palette = "blue_red", border = "none")
}

heatmap_palette_colors <- function(style, theme = "publication") {
  if (identical(value_or(style, "theme"), "theme")) {
    style <- heatmap_theme_options(theme)$palette
  }
  if (identical(style, "viridis")) return(grDevices::hcl.colors(100, "viridis"))
  if (identical(style, "magma")) return(grDevices::hcl.colors(100, "magma"))
  if (identical(style, "green_purple")) return(colorRampPalette(c("#1b7837", "white", "#762a83"))(100))
  if (identical(style, "navy_orange")) return(colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100))
  if (identical(style, "gray_red")) return(colorRampPalette(c("#3b3b3b", "white", "#b2182b"))(100))
  colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
}

heatmap_color_breaks <- function(mat, colors) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(seq(-1, 1, length.out = length(colors) + 1))
  lo <- min(vals, na.rm = TRUE)
  hi <- max(vals, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || identical(lo, hi)) {
    lo <- lo - 1
    hi <- hi + 1
  }
  if (lo < 0 && hi > 0) {
    lim <- max(abs(c(lo, hi)))
    return(seq(-lim, lim, length.out = length(colors) + 1))
  }
  seq(lo, hi, length.out = length(colors) + 1)
}

heatmap_border_color <- function(style, theme = "publication") {
  if (identical(value_or(style, "theme"), "theme")) {
    style <- heatmap_theme_options(theme)$border
  }
  if (identical(style, "light")) "#d9d9d9" else NA
}

heatmap_font_family <- function(style) {
  style <- value_or(style, "sans")
  if (style %in% c("sans", "serif", "mono")) style else "sans"
}

volcano_palette_colors <- function(style) {
  style <- value_or(style, "publication")
  if (identical(style, "minimal")) {
    return(c(nonsig = "#b8bec6", up = "#d95f02", down = "#1b9e77", selected = "#111111"))
  }
  if (identical(style, "colorblind")) {
    return(c(nonsig = "#b9b9b9", up = "#d55e00", down = "#0072b2", selected = "#000000"))
  }
  if (identical(style, "dark")) {
    return(c(nonsig = "#667085", up = "#ff6b6b", down = "#4dabf7", selected = "#ffffff"))
  }
  if (identical(style, "soft")) {
    return(c(nonsig = "#c8cdd4", up = "#c44e52", down = "#4c72b0", selected = "#2f2f2f"))
  }
  c(nonsig = "#aeb6c2", up = "#b2182b", down = "#2166ac", selected = "#111111")
}

volcano_theme_options <- function(style) {
  style <- value_or(style, "publication")
  if (identical(style, "dark")) {
    return(list(bg = "#111827", fg = "#f9fafb", grid = "#374151", axis = "#e5e7eb"))
  }
  if (identical(style, "minimal")) {
    return(list(bg = "#ffffff", fg = "#202020", grid = "#eeeeee", axis = "#444444"))
  }
  list(bg = "#ffffff", fg = "#202020", grid = "#e6e9ef", axis = "#3a3a3a")
}

prepare_volcano_df <- function(df, p_col, p_cutoff, lfc_cutoff) {
  if (is.null(df) || !nrow(df) || !("log2FoldChange" %in% colnames(df)) || !(p_col %in% colnames(df))) {
    return(NULL)
  }
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  lfc <- suppressWarnings(as.numeric(df$log2FoldChange))
  genes <- as.character(df$gene_label)
  keep <- !is.na(genes) & nzchar(genes) & !is.na(pvals) & !is.na(lfc)
  out <- data.frame(
    gene_label = genes[keep],
    log2FoldChange = lfc[keep],
    pvalue_metric = pvals[keep],
    neg_log10_p = -log10(pmax(pvals[keep], .Machine$double.xmin)),
    stringsAsFactors = FALSE
  )
  out$status <- "Not significant"
  out$status[out$pvalue_metric <= p_cutoff & out$log2FoldChange >= lfc_cutoff] <- "Up"
  out$status[out$pvalue_metric <= p_cutoff & out$log2FoldChange <= -lfc_cutoff] <- "Down"
  out
}

volcano_label_genes <- function(plot_df, mode, manual_genes, max_labels) {
  if (is.null(plot_df) || !nrow(plot_df)) return(character(0))
  mode <- value_or(mode, "top")
  max_labels <- max(0, as.integer(value_or(max_labels, 20)))
  if (identical(mode, "none") || max_labels == 0) return(character(0))
  if (identical(mode, "manual")) {
    return(intersect(value_or(manual_genes, character(0)), plot_df$gene_label))
  }
  sig <- plot_df[plot_df$status %in% c("Up", "Down"), , drop = FALSE]
  if (!nrow(sig)) return(character(0))
  sig <- sig[order(-sig$neg_log10_p, -abs(sig$log2FoldChange)), , drop = FALSE]
  head(sig$gene_label, max_labels)
}

draw_volcano_plot <- function(plot_df, labels, title, subtitle, p_cutoff, lfc_cutoff, style, palette_style,
                              point_size, point_alpha, label_size, show_thresholds, show_grid, font_family) {
  validate(need(!is.null(plot_df) && nrow(plot_df) > 0, "No differential expression rows are available for this volcano plot."))
  pal <- volcano_palette_colors(palette_style)
  theme <- volcano_theme_options(style)
  label_hits <- plot_df[plot_df$gene_label %in% labels, , drop = FALSE]

  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggrepel", quietly = TRUE)) {
    plot_df$status <- factor(plot_df$status, levels = c("Down", "Not significant", "Up"))
    label_hits$status <- factor(label_hits$status, levels = c("Down", "Not significant", "Up"))
    threshold_y <- -log10(pmax(p_cutoff, .Machine$double.xmin))
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = log2FoldChange, y = neg_log10_p, color = status)
    )
    if (isTRUE(show_grid)) {
      p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = theme$grid, size = 0.25),
                              panel.grid.minor = ggplot2::element_blank())
    } else {
      p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    }
    if (isTRUE(show_thresholds)) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), color = "#7f8792", linetype = "dashed", size = 0.45) +
        ggplot2::geom_hline(yintercept = threshold_y, color = "#7f8792", linetype = "dashed", size = 0.45)
    }
    p <- p +
      ggplot2::geom_point(alpha = as.numeric(value_or(point_alpha, 0.72)), size = as.numeric(value_or(point_size, 2.5))) +
      ggplot2::scale_color_manual(
        values = c("Down" = pal[["down"]], "Not significant" = pal[["nonsig"]], "Up" = pal[["up"]]),
        breaks = c("Up", "Down", "Not significant")
      ) +
      ggplot2::labs(
        title = title,
        x = "log2 fold change",
        y = paste0("-log10(", subtitle, ")"),
        color = NULL
      ) +
      ggplot2::theme_classic(base_family = font_family) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = theme$bg, color = NA),
        panel.background = ggplot2::element_rect(fill = theme$bg, color = NA),
        axis.text = ggplot2::element_text(color = theme$axis),
        axis.title = ggplot2::element_text(color = theme$fg, face = "bold"),
        plot.title = ggplot2::element_text(color = theme$fg, face = "bold", hjust = 0, size = 14),
        legend.position = "top",
        legend.justification = "left",
        legend.text = ggplot2::element_text(color = theme$fg),
        legend.background = ggplot2::element_rect(fill = theme$bg, color = NA)
      )
    if (nrow(label_hits)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_hits,
        ggplot2::aes(label = gene_label),
        color = pal[["selected"]],
        size = max(2.5, as.numeric(value_or(label_size, 2.5))),
        box.padding = 0.55,
        point.padding = 0.35,
        min.segment.length = 0,
        segment.color = "#7f8792",
        segment.alpha = 0.75,
        segment.size = 0.3,
        max.overlaps = Inf,
        seed = 8,
        force = 2.5,
        force_pull = 0.35
      )
    }
    print(p)
    return(invisible(p))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(bg = theme$bg, fg = theme$fg, col.axis = theme$axis, col.lab = theme$fg, col.main = theme$fg, family = font_family)
  par(mar = c(5, 5, 4.5, 2) + 0.1)

  x <- plot_df$log2FoldChange
  y <- plot_df$neg_log10_p
  xlim <- range(c(x, -lfc_cutoff, lfc_cutoff), finite = TRUE)
  ylim <- range(c(0, y), finite = TRUE)
  xpad <- max(0.5, diff(xlim) * 0.08)
  ypad <- max(0.5, diff(ylim) * 0.08)
  xlim <- xlim + c(-xpad, xpad)
  ylim <- c(0, ylim[2] + ypad)

  plot(x, y, type = "n", xlim = xlim, ylim = ylim, xlab = "log2 fold change",
       ylab = paste0("-log10(", subtitle, ")"), main = title, cex.main = 1.15, font.main = 2, bty = "l")
  if (isTRUE(show_grid)) abline(h = pretty(ylim), v = pretty(xlim), col = theme$grid, lwd = 0.8)
  if (isTRUE(show_thresholds)) {
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = "#7f8792", lty = 2, lwd = 1.1)
    abline(h = -log10(pmax(p_cutoff, .Machine$double.xmin)), col = "#7f8792", lty = 2, lwd = 1.1)
  }
  cols <- rep(pal[["nonsig"]], nrow(plot_df))
  cols[plot_df$status == "Up"] <- pal[["up"]]
  cols[plot_df$status == "Down"] <- pal[["down"]]
  points(x, y, pch = 16, col = adjustcolor(cols, alpha.f = as.numeric(value_or(point_alpha, 0.72))), cex = as.numeric(value_or(point_size, 2.5)))
  if (nrow(label_hits)) {
    text(label_hits$log2FoldChange, label_hits$neg_log10_p, labels = label_hits$gene_label,
         pos = 3, cex = as.numeric(value_or(label_size, 2.5)), col = pal[["selected"]], xpd = NA)
  }
  legend("topright", legend = c("Up", "Down", "Not significant"),
         col = c(pal[["up"]], pal[["down"]], pal[["nonsig"]]), pch = 16, bty = "n", text.col = theme$fg, cex = 0.85)
}

make_enrichr_gene_lists <- function(df, p_col, p_cutoff, lfc_cutoff) {
  if (is.null(df) || !nrow(df) || !(p_col %in% colnames(df)) || !("log2FoldChange" %in% colnames(df))) {
    return(list(up = character(0), down = character(0)))
  }
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  lfc <- suppressWarnings(as.numeric(df$log2FoldChange))
  genes <- as.character(df$gene_label)
  keep <- !is.na(genes) & nzchar(genes) & !is.na(pvals) & !is.na(lfc)
  up <- genes[keep & pvals <= p_cutoff & lfc >= lfc_cutoff]
  down <- genes[keep & pvals <= p_cutoff & lfc <= (-lfc_cutoff)]
  list(up = unique(up), down = unique(down))
}

qc_subtabs <- list(
  tabPanel(
    "Initial QC",
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectInput("sample", "Sample", choices = samples, selected = first_or_null(samples), selectize = FALSE),
        checkboxInput("show_trimmed", "Show cutadapt-trimmed QC", value = default_show_trimmed),
        tags$hr(),
        helpText("This tab renders FastQC and FastQ Screen HTML reports for both reads.")
      ),
      mainPanel(
        width = 10,
        uiOutput("qc_status_ui"),
        tabsetPanel(
          tabPanel("R1 FastQC", uiOutput("r1_fastqc_ui")),
          tabPanel("R1 Screen", uiOutput("r1_screen_ui")),
          tabPanel("R2 FastQC", uiOutput("r2_fastqc_ui")),
          tabPanel("R2 Screen", uiOutput("r2_screen_ui"))
        )
      )
    )
  ),
  tabPanel(
    "STAR",
    sidebarLayout(
      sidebarPanel(
        selectInput("star_sample", "STAR sample", choices = samples, selected = first_or_null(samples), selectize = FALSE),
        selectInput("star_sample_sort_col", "Sort sample columns by", choices = comparison_columns, selected = if ("treatment" %in% comparison_columns) "treatment" else first_or_null(comparison_columns), selectize = FALSE),
        tags$hr(),
        helpText("Compact alignment metrics pulled from STAR summary output.")
      ),
      mainPanel(
        uiOutput("star_status_ui"),
        h4("Alignment Summary Across Samples"),
        table_widget("star_summary_table"),
        download_table_button("download_star_summary", "Download STAR summary"),
        tags$hr(style = "margin: 10px 0 8px 0;"),
        h4(style = "margin-top: 0;", "Selected Sample"),
        tableOutput("star_sample_table"),
        download_table_button("download_star_sample", "Download selected STAR sample")
      )
    )
  ),
  tabPanel(
    "FeatureCounts QC",
    sidebarLayout(
      sidebarPanel(
        selectInput("featurecounts_qc_sample", "FeatureCounts sample", choices = featurecounts_sample_choices, selected = first_or_null(featurecounts_sample_choices), selectize = FALSE),
        selectInput("featurecounts_qc_sample_sort_col", "Sort sample columns by", choices = comparison_columns, selected = if ("treatment" %in% comparison_columns) "treatment" else first_or_null(comparison_columns), selectize = FALSE)
      ),
      mainPanel(
        uiOutput("featurecounts_qc_status_ui"),
        h4("FeatureCounts Summary Across Samples"),
        table_widget("featurecounts_summary_table"),
        download_table_button("download_featurecounts_summary", "Download FeatureCounts summary"),
        tags$hr(style = "margin: 10px 0 8px 0;"),
        h4(style = "margin-top: 0;", "Selected Sample"),
        tableOutput("featurecounts_sample_table"),
        download_table_button("download_featurecounts_sample", "Download selected FeatureCounts sample")
      )
    )
  )
)

counts_subtabs <- list(
  tabPanel(
    "Raw Counts",
    sidebarLayout(
      sidebarPanel(
        selectizeInput("gene_query", "Search gene", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body")),
        selectInput("featurecounts_label_mode", "Gene label", choices = c("Gene ID" = "gene_id", "Gene name" = "gene_name"), selected = "gene_id", selectize = FALSE),
        uiOutput("featurecounts_convert_ui"),
        tags$hr(),
        helpText("Search the raw featureCounts matrix. Gene-name mode saves a converted table in the counts folder and combines duplicate gene names by summing counts.")
      ),
      mainPanel(
        uiOutput("featurecounts_status_ui"),
        table_widget("gene_search_table"),
        download_table_button("download_raw_counts_table", "Download raw counts table")
      )
    )
  )
)

kallisto_filter_columns <- c("gene_symbol", "transcript_name", "transcript_id", "gene_id")
counts_subtabs <- c(
  counts_subtabs,
  list(
    tabPanel(
      "RSEM",
      sidebarLayout(
        sidebarPanel(
          selectInput("rsem_type", "Show", choices = c("Genes" = "genes", "Isoforms" = "isoforms"), selected = "genes", selectize = FALSE),
          selectInput("rsem_metric", "Metric", choices = c("TPM", "expected_count", "FPKM"), selected = "TPM", selectize = FALSE),
          selectInput("rsem_label_mode", "Gene label", choices = c("Gene ID" = "gene_id", "Gene name" = "gene_name"), selected = "gene_id", selectize = FALSE),
          uiOutput("rsem_convert_ui"),
          selectizeInput("rsem_query", "Select gene/transcript of interest", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body")),
          actionButton("rsem_load_matrix", "Load RSEM matrix", class = "btn-primary"),
          tags$hr(),
          helpText("Loads saved RSEM matrices from the counts folder. Gene-name conversion writes a new matrix file beside the original.")
        ),
        mainPanel(
          h4("RSEM Matrix"),
          uiOutput("rsem_status_ui"),
          table_widget("rsem_table"),
          download_table_button("download_rsem_table", "Download RSEM table")
        )
      )
    ),
    tabPanel(
      "DESeq2 Normalized Counts",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "deseq_compare_col",
            "Comparison column",
            choices = c("NA" = "NA", comparison_columns),
            selected = if ("treatment" %in% comparison_columns) "treatment" else "NA", selectize = FALSE),
          uiOutput("deseq_treatment_ui"),
          uiOutput("deseq_control_ui"),
          selectInput("deseq_label_mode", "Gene label", choices = c("Gene ID" = "gene_id", "Gene name" = "gene_name"), selected = "gene_id", selectize = FALSE),
          selectizeInput("deseq_gene_query", "Select gene", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body")),
          uiOutput("deseq_convert_ui"),
          tags$hr(),
          helpText("Shows DESeq2 normalized counts. Gene-name mode saves a duplicate-combined converted table beside the original normalized-counts file.")
        ),
        mainPanel(
          uiOutput("deseq_status_ui"),
          table_widget("deseq_counts_table"),
          download_table_button("download_deseq_counts_table", "Download DESeq2 counts table")
        )
      )
    ),
    tabPanel(
      "Kallisto",
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "kallisto_view_mode",
            "View mode",
            choices = c("Single sample matrix" = "sample", "Transcript matrix view" = "matrix"),
            selected = "matrix"
          ),
          conditionalPanel(
            "input.kallisto_view_mode == 'sample'",
            selectInput("kallisto_sample", "Kallisto sample", choices = kallisto_sample_names, selected = first_or_null(kallisto_sample_names), selectize = FALSE),
            selectInput("kallisto_sample_filter_col", "Filter by", choices = c("All transcripts" = "all", kallisto_filter_columns), selected = "all", selectize = FALSE),
            selectizeInput("kallisto_sample_filter_value", "Select value", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body"))
          ),
          conditionalPanel(
            "input.kallisto_view_mode == 'matrix'",
            selectInput("kallisto_matrix_metric", "Matrix", choices = c("TPM" = "tpm", "Estimated counts" = "est_counts"), selected = "tpm", selectize = FALSE),
            selectInput("kallisto_label_mode", "Label", choices = c("Transcript ID" = "target_id", "Gene name" = "gene_name"), selected = "target_id", selectize = FALSE),
            uiOutput("kallisto_convert_ui"),
            selectInput("kallisto_filter_col", "Filter by", choices = c("All rows" = "all", "target_id", "gene_name", "transcript_id", "gene_id", "gene_symbol"), selected = "all", selectize = FALSE),
            selectizeInput("kallisto_filter_value", "Select value", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body")),
            tags$div(class = "tiny-note", "Kallisto matrix view reads saved matrices from the counts folder.")
          ),
          actionButton("kallisto_load_table", "Load Kallisto table", class = "btn-primary"),
          tags$hr(),
          helpText("Single sample matrix shows transcript-level abundance for one sample. Transcript matrix view shows matching transcripts across all samples.")
        ),
        mainPanel(
          h4("Kallisto Transcript Abundance"),
          uiOutput("kallisto_status_ui"),
          table_widget("kallisto_table"),
          download_table_button("download_kallisto_table", "Download Kallisto table")
        )
      )
    )
  )
)

app_tabs <- list(
  tabPanel("QC", do.call(tabsetPanel, qc_subtabs)),
  tabPanel("Counts", do.call(tabsetPanel, c(list(id = "counts_subtab"), counts_subtabs))),
  tabPanel(
    "Differential Expression",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "deg_compare_col",
          "Comparison column",
          choices = c("NA" = "NA", comparison_columns),
          selected = if ("treatment" %in% comparison_columns) "treatment" else "NA", selectize = FALSE),
        uiOutput("deg_treatment_ui"),
        uiOutput("deg_control_ui"),
        selectizeInput("deg_gene_query", "Select gene", choices = NULL, selected = NULL, multiple = FALSE, options = list(dropdownParent = "body")),
        uiOutput("deg_label_mode_ui"),
        uiOutput("deg_convert_ui"),
        selectInput("deg_p_col", "P-value column", choices = c("padj", "pvalue"), selected = "padj", selectize = FALSE),
        numericInput("deg_p_cutoff", "P-value cutoff", value = 0.05, min = 0, max = 1, step = 0.001),
        numericInput("deg_lfc_cutoff", "Absolute log2FC cutoff", value = 0, min = 0, step = 0.1),
        helpText("These cutoffs only build the Enrichr up/down gene lists below; they do not filter the DEG table."),
        selectInput(
          "deg_sort_mode",
          "Sort DEGs",
          choices = c(
            "Upregulated genes at top" = "up",
            "Downregulated genes at top" = "down"
          ),
          selected = "up", selectize = FALSE)
      ),
      mainPanel(
        uiOutput("deg_status_ui"),
        tags$div(
          style = "max-height: 620px; overflow-y: auto; overflow-x: auto;",
          table_widget("deg_table")
        ),
        tags$hr(),
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center;",
          h4("Upregulated Genes For Enrichr"),
          tags$a(href = "https://maayanlab.cloud/OxEnrichr/", target = "_blank", "Open OxEnrichr")
        ),
        tags$p("Copy this list, open OxEnrichr, and paste the genes into a new submission."),
        textAreaInput("deg_up_genes", NULL, value = "", width = "100%", height = "180px"),
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center;",
          h4("Downregulated Genes For Enrichr"),
          tags$a(href = "https://maayanlab.cloud/OxEnrichr/", target = "_blank", "Open OxEnrichr")
        ),
        tags$p("Copy this list, open OxEnrichr, and paste the genes into a new submission."),
        textAreaInput("deg_down_genes", NULL, value = "", width = "100%", height = "180px")
      )
    )
  ),
  tabPanel(
    "Plots",
    sidebarLayout(
      sidebarPanel(
        conditionalPanel(
          "input.plot_subtab != 'PCA'",
          selectInput(
            "plot_compare_col",
            "Comparison column",
            choices = c("NA" = "NA", comparison_columns),
            selected = if ("treatment" %in% comparison_columns) "treatment" else "NA", selectize = FALSE),
          uiOutput("plot_treatment_ui"),
          uiOutput("plot_control_ui"),
          uiOutput("plot_deseq_label_mode_ui")
        ),
        conditionalPanel(
          "input.plot_subtab == 'PCA'",
          selectInput(
            "pca_source",
            "PCA values",
            choices = c("DESeq2 normalized counts" = "deseq", "featureCounts count matrix" = "counts"),
            selected = "deseq", selectize = FALSE),
          selectInput(
            "pca_color_col",
            "Color by",
            choices = c("None" = "__none__", setNames(comparison_columns, comparison_columns)),
            selected = if ("treatment" %in% comparison_columns) "treatment" else "__none__", selectize = FALSE),
          uiOutput("pca_include_values_ui"),
          checkboxInput("pca_label_samples", "Label samples", value = TRUE),
          numericInput("pca_plot_width", "Display width (px)", value = 800, min = 400, max = 2400, step = 50),
          numericInput("pca_plot_height", "Display height (px)", value = 600, min = 300, max = 2400, step = 50)
        ),
        conditionalPanel(
          "input.plot_subtab == 'Volcano'",
          selectInput("volcano_p_col", "P-value column", choices = c("padj", "pvalue"), selected = "padj", selectize = FALSE),
          numericInput("volcano_p_cutoff", "P-value cutoff", value = 0.05, min = 0, max = 1, step = 0.001),
          numericInput("volcano_lfc_cutoff", "Absolute log2FC cutoff", value = 1, min = 0, step = 0.1),
          selectInput(
            "volcano_style",
            "Volcano style",
            choices = c("Clean publication" = "publication", "Minimal" = "minimal", "Dark" = "dark", "Colorblind" = "colorblind", "Soft" = "soft"),
            selected = "publication", selectize = FALSE),
          selectInput(
            "volcano_palette",
            "Volcano colors",
            choices = c("Publication" = "publication", "Minimal" = "minimal", "Dark" = "dark", "Colorblind" = "colorblind", "Soft" = "soft"),
            selected = "publication", selectize = FALSE),
          selectInput(
            "volcano_label_mode",
            "Gene labels",
            choices = c("Top significant genes" = "top", "Select genes manually" = "manual", "No labels" = "none"),
            selected = "top", selectize = FALSE),
          conditionalPanel(
            "input.volcano_label_mode == 'manual'",
            selectizeInput("volcano_label_genes", "Genes to label", choices = NULL, selected = NULL, multiple = TRUE, options = list(dropdownParent = "body"))
          ),
          conditionalPanel(
            "input.volcano_label_mode == 'top'",
            numericInput("volcano_max_labels", "Maximum labels", value = 20, min = 0, max = 100, step = 1)
          ),
          numericInput("volcano_point_size", "Point size", value = 2.5, min = 0.2, max = 8, step = 0.1),
          numericInput("volcano_point_alpha", "Point opacity", value = 0.72, min = 0.1, max = 1, step = 0.05),
          numericInput("volcano_label_size", "Label size", value = 2.5, min = 0.5, max = 8, step = 0.1),
          selectInput(
            "volcano_font_family",
            "Font family",
            choices = c("Sans" = "sans", "Serif" = "serif", "Mono" = "mono"),
            selected = "sans", selectize = FALSE),
          checkboxInput("volcano_show_thresholds", "Show cutoff lines", value = TRUE),
          checkboxInput("volcano_show_grid", "Show grid", value = TRUE),
          numericInput("volcano_plot_width", "Display width (px)", value = 800, min = 400, max = 2400, step = 50),
          numericInput("volcano_plot_height", "Display height (px)", value = 600, min = 300, max = 2400, step = 50),
          selectInput(
            "volcano_save_mode",
            "Saved image size",
            choices = c("Match display" = "match", "High-res 2x" = "double", "High-res 3x" = "triple", "High-res 5x" = "quintuple", "Custom" = "custom"),
            selected = "double", selectize = FALSE),
          conditionalPanel(
            "input.volcano_save_mode == 'custom'",
            numericInput("volcano_save_width", "Custom saved width (px)", value = 1800, min = 600, max = 8000, step = 100),
            numericInput("volcano_save_height", "Custom saved height (px)", value = 1400, min = 600, max = 8000, step = 100),
            numericInput("volcano_save_res", "Custom saved DPI", value = 192, min = 72, max = 600, step = 10)
          )
        ),
        conditionalPanel(
          "input.plot_subtab == 'Heatmap'",
          radioButtons(
            "heatmap_gene_mode",
            "Heatmap genes",
            choices = c("Select genes manually" = "manual", "Top up/down DEGs" = "topdeg"),
            selected = "topdeg"
          ),
          selectInput(
            "heatmap_source",
            "Heatmap values",
            choices = c("DESeq2 normalized counts" = "deseq", "RSEM TPM" = "rsem"),
            selected = "deseq", selectize = FALSE),
          conditionalPanel(
            "input.heatmap_gene_mode == 'manual'",
            selectizeInput("heatmap_genes", "Genes to plot", choices = NULL, selected = NULL, multiple = TRUE, options = list(dropdownParent = "body"))
          ),
          conditionalPanel(
            "input.heatmap_gene_mode == 'topdeg'",
            selectInput("heatmap_rank_p_col", "Rank by p-value column", choices = c("padj", "pvalue"), selected = "padj", selectize = FALSE),
            numericInput("heatmap_total_genes", "Total genes", value = 20, min = 2, step = 2),
            checkboxInput("heatmap_equal_split", "Equal numbers of up and down", value = TRUE)
          ),
          selectizeInput(
            "heatmap_annotation_cols",
            "Column annotations",
            choices = comparison_columns,
            selected = if ("treatment" %in% comparison_columns) "treatment" else character(0),
            multiple = TRUE, options = list(dropdownParent = "body")),
          checkboxInput("heatmap_cluster_rows", "Cluster rows", value = TRUE),
          checkboxInput("heatmap_cluster_cols", "Cluster columns", value = TRUE),
          conditionalPanel(
            "!input.heatmap_cluster_cols",
            selectInput(
              "heatmap_order_col",
              "Order samples by",
              choices = comparison_columns,
              selected = if ("treatment" %in% comparison_columns) "treatment" else comparison_columns[1], selectize = FALSE)
          ),
          checkboxInput("heatmap_scale_rows", "Scale heatmap rows", value = TRUE),
          selectInput(
            "heatmap_theme",
            "Heatmap style",
            choices = c("Clean publication" = "publication", "Minimal" = "minimal", "Compact" = "compact", "Fully labeled" = "labeled", "High contrast" = "contrast"),
            selected = "publication", selectize = FALSE),
          selectInput(
            "heatmap_palette",
            "Heatmap colors",
            choices = c("Use style default" = "theme", "Blue-white-red" = "blue_red", "Viridis" = "viridis", "Magma" = "magma", "Green-white-purple" = "green_purple", "Blue-white-red vivid" = "navy_orange", "Gray-white-red" = "gray_red"),
            selected = "theme", selectize = FALSE),
          selectInput(
            "heatmap_border_style",
            "Cell borders",
            choices = c("Use style default" = "theme", "None" = "none", "Light grid" = "light"),
            selected = "theme", selectize = FALSE),
          selectInput(
            "heatmap_font_family",
            "Font family",
            choices = c("Sans" = "sans", "Serif" = "serif", "Mono" = "mono"),
            selected = "sans", selectize = FALSE),
          selectInput("heatmap_row_font_size", "Gene label size", choices = c("Small" = 7, "Medium" = 9, "Large" = 11, "XL" = 13), selected = 9, selectize = FALSE),
          selectInput("heatmap_col_font_size", "Sample label size", choices = c("Small" = 8, "Medium" = 10, "Large" = 12, "XL" = 14), selected = 10, selectize = FALSE),
          selectInput("heatmap_col_angle", "Sample label angle", choices = c("0" = "0", "45" = "45", "90" = "90", "315" = "315"), selected = "45", selectize = FALSE),
          checkboxInput("heatmap_show_row_names", "Show gene labels", value = TRUE),
          checkboxInput("heatmap_show_col_names", "Show sample labels", value = TRUE),
          checkboxInput("heatmap_show_annotation_names", "Show annotation track labels", value = FALSE),
          selectInput("heatmap_distance", "Clustering distance", choices = c("Euclidean" = "euclidean", "Correlation" = "correlation", "Manhattan" = "manhattan"), selected = "euclidean", selectize = FALSE),
          selectInput("heatmap_cluster_method", "Clustering method", choices = c("Complete" = "complete", "Average" = "average", "Single" = "single", "Ward D2" = "ward.D2"), selected = "complete", selectize = FALSE),
          numericInput("heatmap_plot_width", "Display width (px)", value = 900, min = 400, max = 2400, step = 50),
          numericInput("heatmap_plot_height", "Display height (px)", value = 700, min = 300, max = 2400, step = 50),
          selectInput(
            "heatmap_save_mode",
            "Saved image size",
            choices = c("Match display" = "match", "High-res 2x" = "double", "High-res 3x" = "triple", "High-res 5x" = "quintuple", "Custom" = "custom"),
            selected = "double", selectize = FALSE),
          conditionalPanel(
            "input.heatmap_save_mode == 'custom'",
            numericInput("heatmap_save_width", "Custom saved width (px)", value = 1800, min = 600, max = 8000, step = 100),
            numericInput("heatmap_save_height", "Custom saved height (px)", value = 1400, min = 600, max = 8000, step = 100),
            numericInput("heatmap_save_res", "Custom saved DPI", value = 96, min = 72, max = 600, step = 10)
          )
        )
      ),
      mainPanel(
        uiOutput("plot_status_ui"),
        tabsetPanel(
          id = "plot_subtab",
          tabPanel("PCA", plotOutput("pca_plot", height = "620px")),
          tabPanel("Volcano", textInput("volcano_filename", "Filename", value = ""), actionButton("save_volcano_btn", "Save volcano plot"), plotOutput("deg_volcano_plot", height = "620px")),
          tabPanel("Heatmap", textInput("heatmap_filename", "Filename", value = ""), actionButton("save_heatmap_btn", "Save heatmap"), plotOutput("deg_heatmap_plot", height = "700px"))
        )
      )
    )
  ),
  tabPanel(
    "GSEA",
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "gsea_compare_col",
              "Comparison column",
              choices = c("NA" = "NA", comparison_columns),
              selected = if ("treatment" %in% comparison_columns) "treatment" else "NA", selectize = FALSE),
            uiOutput("gsea_treatment_ui"),
            uiOutput("gsea_control_ui"),
            uiOutput("gsea_collection_ui"),
            uiOutput("gsea_pathway_ui"),
            tags$hr(),
            helpText("Select a tested comparison, then a pathway collection and pathway to view GSEA outputs.")
          ),
          mainPanel(
            uiOutput("gsea_status_ui"),
            h4("Collection Summary Plots"),
            uiOutput("gsea_summary_plots_ui"),
            tags$hr(),
            tabsetPanel(
              tabPanel(
                "All Pathways",
                tags$div(
                  style = "max-height: 520px; overflow-y: auto; overflow-x: auto;",
                  table_widget("gsea_all_pathways_table")
                )
              ),
              tabPanel(
                "Individual Pathway",
                table_widget("gsea_pathway_table"),
                tags$hr(),
                uiOutput("gsea_pathway_plots_ui")
              )
            )
          )
        )
  )
)

ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      body {
        background:
          radial-gradient(circle at top left, rgba(64, 142, 255, 0.18), transparent 34%),
          radial-gradient(circle at top right, rgba(24, 196, 140, 0.14), transparent 28%),
          linear-gradient(180deg, #f4f7fb 0%, #edf2f7 100%);
        color: #10233a;
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
      }
      .container-fluid {
        max-width: 1600px;
        padding: 24px 28px 36px 28px;
      }
      .app-shell {
        background: rgba(255, 255, 255, 0.72);
        border: 1px solid rgba(255, 255, 255, 0.9);
        box-shadow: 0 18px 45px rgba(30, 57, 90, 0.12);
        backdrop-filter: blur(14px);
        -webkit-backdrop-filter: blur(14px);
        border-radius: 28px;
        overflow: hidden;
      }
      .hero {
        position: relative;
        padding: 42px 40px 34px 40px;
        background:
          linear-gradient(135deg, rgba(8, 41, 84, 0.96) 0%, rgba(12, 69, 132, 0.92) 54%, rgba(20, 132, 107, 0.86) 100%);
        color: white;
      }
      .hero:after {
        content: '';
        position: absolute;
        inset: 0;
        background:
          radial-gradient(circle at 85% 18%, rgba(255,255,255,0.16), transparent 20%),
          radial-gradient(circle at 18% 120%, rgba(255,255,255,0.10), transparent 24%);
        pointer-events: none;
      }
      .hero-kicker {
        font-size: 14px;
        opacity: 0.88;
        margin-bottom: 10px;
      }
      .hero-title {
        font-size: 40px;
        line-height: 1.08;
        font-weight: 800;
        margin: 0 0 10px 0;
      }
      .hero-subtitle {
        max-width: 880px;
        font-size: 15px;
        line-height: 1.6;
        opacity: 0.92;
        margin: 0;
      }
      .hero-topbar {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 20px;
      }
      .hero-copy {
        min-width: 0;
      }
      .hero-logos {
        display: flex;
        align-items: center;
        gap: 28px;
        flex-shrink: 0;
      }
      .hero-logo {
        height: 168px;
        width: auto;
        object-fit: contain;
        background: transparent !important;
        border: none !important;
        padding: 0;
        border-radius: 0;
        box-shadow: none !important;
      }
      .main-tabs {
        padding: 20px 22px 26px 22px;
      }
      .main-tabs > .tabbable > .nav-tabs {
        border-bottom: none;
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        margin-bottom: 18px;
      }
      .main-tabs > .tabbable > .nav-tabs > li {
        margin-bottom: 0;
      }
      .main-tabs > .tabbable > .nav-tabs > li > a {
        border: none !important;
        border-radius: 999px !important;
        background: rgba(255,255,255,0.8);
        color: #27425f !important;
        font-weight: 700;
        padding: 11px 18px;
        box-shadow: 0 6px 16px rgba(28, 54, 88, 0.08);
        transition: all 0.18s ease;
      }
      .main-tabs > .tabbable > .nav-tabs > li.active > a,
      .main-tabs > .tabbable > .nav-tabs > li.active > a:hover,
      .main-tabs > .tabbable > .nav-tabs > li.active > a:focus {
        background: linear-gradient(135deg, #0f62c6, #19a974) !important;
        color: white !important;
        box-shadow: 0 10px 22px rgba(24, 95, 185, 0.28);
      }
      .main-tabs > .tabbable > .nav-tabs > li > a:hover {
        background: rgba(240, 247, 255, 0.98) !important;
        transform: translateY(-1px);
      }
      .tab-content .tabbable > .nav-tabs {
        border-bottom: none;
        display: inline-flex;
        flex-wrap: wrap;
        gap: 8px;
        margin: 6px 0 18px 0;
        padding: 8px;
        background: rgba(237, 244, 251, 0.95);
        border: 1px solid rgba(209, 223, 239, 0.95);
        border-radius: 18px;
      }
      .tab-content .tabbable > .nav-tabs > li {
        margin-bottom: 0;
      }
      .tab-content .tabbable > .nav-tabs > li > a {
        border: none !important;
        border-radius: 14px !important;
        background: transparent;
        color: #48627d !important;
        font-weight: 700;
        padding: 10px 14px;
        box-shadow: none;
        transition: all 0.16s ease;
      }
      .tab-content .tabbable > .nav-tabs > li.active > a,
      .tab-content .tabbable > .nav-tabs > li.active > a:hover,
      .tab-content .tabbable > .nav-tabs > li.active > a:focus {
        background: linear-gradient(135deg, #ffffff, #f4f9ff) !important;
        color: #143150 !important;
        box-shadow: 0 8px 16px rgba(35, 63, 99, 0.10);
      }
      .tab-content .tabbable > .nav-tabs > li > a:hover {
        background: rgba(255,255,255,0.65) !important;
        transform: none;
      }
      .tab-content {
        padding-top: 6px;
      }
      .well, .sidebar-panel, .main-panel, .panel, .tabbable > .tab-content {
        background: transparent;
        border: none;
        box-shadow: none;
      }
      .row > .col-sm-4 > .well,
      .row > .col-sm-3 > .well,
      .row > .col-md-4 > .well,
      .row > .col-md-3 > .well {
        background: rgba(248, 251, 255, 0.96);
        border: 1px solid rgba(197, 215, 236, 0.9);
        border-radius: 22px;
        padding: 20px 18px;
        box-shadow: inset 0 1px 0 rgba(255,255,255,0.9);
      }
      .row > .col-sm-8,
      .row > .col-sm-9,
      .row > .col-md-8,
      .row > .col-md-9 {
        background: rgba(255,255,255,0.86);
        border: 1px solid rgba(214, 225, 238, 0.95);
        border-radius: 22px;
        padding: 20px 22px;
        box-shadow: 0 10px 26px rgba(32, 56, 84, 0.08);
      }
      .form-control, .selectize-input, .selectize-control.single .selectize-input, textarea {
        border-radius: 14px !important;
        border: 1px solid #d7e0ea !important;
        box-shadow: none !important;
        min-height: 46px;
        padding-top: 10px;
        padding-bottom: 10px;
        font-size: 15px !important;
      }
      .selectize-dropdown, .selectize-input, .form-control {
        font-size: 15px;
      }
      .control-label {
        font-size: 13px;
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: #5a7088;
        margin-bottom: 8px;
      }
      .btn, .btn-default, .btn-primary {
        border: none !important;
        border-radius: 14px !important;
        background: linear-gradient(135deg, #0f62c6, #15936f) !important;
        color: white !important;
        font-weight: 700;
        padding: 10px 16px;
        box-shadow: 0 10px 20px rgba(17, 94, 177, 0.22);
        transition: transform 0.16s ease, box-shadow 0.16s ease;
      }
      .btn:hover {
        transform: translateY(-1px);
        box-shadow: 0 12px 24px rgba(17, 94, 177, 0.28);
      }
      h4, h5 {
        color: #18314e;
        font-weight: 800;
        letter-spacing: -0.02em;
        font-size: 22px;
      }
      table {
        background: white;
        border-radius: 16px;
        overflow: hidden;
      }
      table.dataTable thead th {
        background: #edf4fb !important;
        color: #304a66 !important;
      }
      table.dataTable thead th,
      table.dataTable tbody td {
        border-right: 1px solid #c7d6e8 !important;
      }
      table.dataTable thead th:first-child,
      table.dataTable tbody td:first-child {
        border-left: 1px solid #c7d6e8 !important;
      }
      table.dataTable thead tr:first-child th {
        border-bottom: 1px solid #aac2de !important;
      }
      table.dataTable tbody tr:nth-child(odd) {
        background: rgba(246, 250, 255, 0.9) !important;
      }
      table.dataTable {
        table-layout: fixed;
      }
      table.dataTable tbody td, table.dataTable thead th {
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        min-width: 118px;
        max-width: 118px;
      }
      table.dataTable tbody td:first-child, table.dataTable thead th:first-child {
        min-width: 150px;
        max-width: 190px;
      }
      .dataTables_wrapper .dataTables_scrollHead {
        position: sticky;
        top: 0;
        z-index: 10;
      }
      .table > thead > tr > th {
        background: #edf4fb;
        color: #304a66;
        border-bottom: none;
        font-size: 14px;
        text-transform: uppercase;
        letter-spacing: 0.06em;
        padding-top: 14px !important;
        padding-bottom: 14px !important;
      }
      .table > tbody > tr:nth-child(odd) {
        background: rgba(246, 250, 255, 0.9);
      }
      .table > tbody > tr > td {
        border-top: 1px solid #ecf1f6;
        font-size: 15px;
        line-height: 1.45;
        color: #22384f;
        padding-top: 12px !important;
        padding-bottom: 12px !important;
      }
      iframe, img, embed {
        border-radius: 18px;
        background: white;
        box-shadow: 0 8px 22px rgba(32, 56, 84, 0.10);
      }
      .qc-report-frame {
        display: block;
        width: 100%;
        max-width: 100%;
        min-width: 0;
      }
      .container-fluid {
        max-width: 1920px;
        padding: 10px 12px 18px 12px;
      }
      .main-tabs {
        padding: 10px 12px 14px 12px;
      }
      .main-tabs > .tabbable > .nav-tabs {
        margin-bottom: 8px;
      }
      .tab-content {
        padding-top: 0;
      }
      .row > .col-sm-8,
      .row > .col-sm-9,
      .row > .col-md-8,
      .row > .col-md-9 {
        padding: 10px 12px;
        border-radius: 14px;
      }
      .row > .col-sm-4 > .well,
      .row > .col-sm-3 > .well,
      .row > .col-md-4 > .well,
      .row > .col-md-3 > .well {
        padding: 12px 12px;
        border-radius: 14px;
      }
      iframe, img, embed {
        border-radius: 10px;
        box-shadow: 0 4px 12px rgba(32, 56, 84, 0.08);
      }
      .shiny-plot-output {
        max-width: 100%;
      }

      /* Dropdown visibility hardening */
      .selectize-control,
      .selectize-control.single,
      .selectize-control.multi {
        position: relative;
        z-index: 30;
        margin-bottom: 12px;
      }
      .form-group:has(.selectize-control.dropdown-active),
      .selectize-control.dropdown-active {
        position: relative;
        z-index: 100000 !important;
      }
      .selectize-input,
      .selectize-input input,
      select.form-control,
      .form-control option {
        color: #132033 !important;
        background: #ffffff !important;
      }
      select.form-control {
        color: #132033 !important;
        background-color: #ffffff !important;
        border: 1px solid #cfd9e6 !important;
        border-radius: 8px !important;
      }
      .selectize-dropdown {
        z-index: 100000 !important;
        background: #ffffff !important;
        color: #132033 !important;
        border: 1px solid #b9c9dc !important;
        box-shadow: 0 14px 30px rgba(20,38,64,.18) !important;
      }
      .selectize-dropdown-content {
        background: #ffffff !important;
        color: #132033 !important;
        max-height: 320px !important;
        overflow-y: auto !important;
      }
      .selectize-dropdown .option,
      .selectize-dropdown .optgroup-header,
      .selectize-dropdown [data-selectable] {
        color: #132033 !important;
        background: #ffffff !important;
        opacity: 1 !important;
        padding: 8px 12px !important;
        line-height: 1.3 !important;
      }
      .selectize-dropdown .option.active,
      .selectize-dropdown [data-selectable].active {
        color: #07111f !important;
        background: #e8f2ff !important;
      }
      .selectize-dropdown .option:hover,
      .selectize-dropdown [data-selectable]:hover {
        color: #07111f !important;
        background: #edf7f4 !important;
      }
      .tab-content,
      .tab-pane,
      .main-panel,
      .sidebar-panel,
      .well,
      .form-group {
        overflow: visible;
      }
      
      .table-scroll-shell {
        width: 100%;
        max-width: 100%;
        overflow: auto;
        -webkit-overflow-scrolling: touch;
        border-radius: 16px;
      }
      .table-scroll-shell .dataTables_wrapper,
      .dataTables_length,
      .dataTables_filter,
      .dataTables_info,
      .dataTables_paginate {
        margin: 8px 0;
        color: #48627d;
        font-size: 13px;
      }
      .dataTables_paginate .paginate_button {
        border-radius: 8px !important;
        border: 1px solid #d7e0ea !important;
        background: #ffffff !important;
        color: #27425f !important;
        margin: 0 2px !important;
      }
      .dataTables_paginate .paginate_button.current {
        background: linear-gradient(135deg, #0f62c6, #15936f) !important;
        color: white !important;
        border-color: transparent !important;
      }
      .dataTables_wrapper {
        width: 100%;
        max-width: 100%;
        overflow-x: auto;
        overflow-y: visible;
      }
      .dataTables_scroll {
        width: 100%;
        max-width: 100%;
        overflow-x: auto;
      }
      .dataTables_scrollBody {
        max-height: min(62vh, 650px) !important;
        overflow: auto !important;
      }
      .tab-pane,
      .tab-content,
      .main-panel,
      .col-sm-8,
      .col-sm-9,
      .col-md-8,
      .col-md-9 {
        min-width: 0;
        max-width: 100%;
        overflow-x: auto;
      }
      .shiny-html-output {
        max-width: 100%;
        overflow-x: auto;
      }
      .shiny-html-output > table,
      .table {
        min-width: max-content;
      }
      .shiny-notification {
        border-radius: 14px;
        box-shadow: 0 10px 30px rgba(26, 49, 78, 0.18);
      }
      .help-block, .helpText, p {
        color: #536a82;
        font-size: 15px;
        line-height: 1.6;
      }
      @media (max-width: 900px) {
        .container-fluid {
          padding: 12px;
        }
        .hero {
          padding: 28px 22px;
        }
        .hero-title {
          font-size: 32px;
        }
        .hero-topbar {
          flex-direction: column;
          align-items: flex-start;
        }
        .main-tabs {
          padding: 14px 12px 18px 12px;
        }
      }
    "))
  ),
  div(
    class = "app-shell",
    div(
      class = "hero",
      div(
        class = "hero-topbar",
        div(
          class = "hero-copy",
          h1(class = "hero-title", "RNA-Seq Results Explorer"),
          div(class = "hero-kicker", "Developed by CSHL's Bioinformatics Shared Resource")
        ),
        div(
          class = "hero-logos",
          if (!is.null(logo_csl_path)) tags$img(class = "hero-logo", src = file.path("logo_csl_asset", basename(logo_csl_path))),
          if (!is.null(logo_path)) tags$img(class = "hero-logo", src = file.path("logo_asset", basename(logo_path)))
        )
      )
    ),
    div(
      class = "main-tabs",
      do.call(tabsetPanel, c(list(id = "main_tabs"), app_tabs))
    )
  )
)

server <- function(input, output, session) {
  gtf_cache <- reactiveValues(mouse = NULL, human = NULL)

  get_gtf_map <- function(species) {
    if (identical(species, "mouse")) {
      if (is.null(gtf_cache$mouse)) gtf_cache$mouse <- read_gtf_gene_map(mouse_gtf_path)
      return(gtf_cache$mouse)
    }
    if (identical(species, "human")) {
      if (is.null(gtf_cache$human)) gtf_cache$human <- read_gtf_gene_map(human_gtf_path)
      return(gtf_cache$human)
    }
    NULL
  }

  save_normalized_counts_gene_names <- function(treatment_value, control_value) {
    base_path <- normalized_counts_path(treatment_value, control_value, "gene_id")
    out_path <- normalized_counts_gene_name_path(treatment_value, control_value)
    if (file.exists(out_path)) {
      return(list(ok = TRUE, message = sprintf("Using saved gene-name normalized counts: %s", out_path), path = out_path))
    }
    df <- read_normalized_counts(base_path)
    if (is.null(df) || !"gene_label" %in% colnames(df) || !nrow(df)) {
      return(list(ok = FALSE, message = "Normalized-counts file was not found or was empty.", path = NA_character_))
    }
    if (!looks_like_gene_id(df$gene_label)) {
      return(list(ok = TRUE, message = "DESeq2 normalized-count labels already look like gene names.", path = base_path))
    }
    species <- detect_species_from_ids(df$gene_label)
    map_df <- get_gtf_map(species)
    conv <- convert_gene_labels(df$gene_label, map_df)
    mapped <- sum(!is.na(conv$values) & nzchar(conv$values) & conv$values != df$gene_label)
    if (mapped == 0) {
      return(list(ok = FALSE, message = sprintf("No DESeq2 normalized-count gene IDs were mapped using the %s GTF.", value_or(species, "detected")), path = NA_character_))
    }
    df$gene_label <- conv$values
    df <- aggregate_display_matrix(df, "gene_label")
    colnames(df)[colnames(df) == "gene_label"] <- "gene_name"
    write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    list(ok = TRUE, message = sprintf("Saved duplicate-combined DESeq2 gene-name normalized counts (%s IDs mapped) to %s.", mapped, out_path), path = out_path)
  }

  aggregate_gene_name_df <- function(df, label_col = "gene_name", keep_cols = character(0)) {
    if (is.null(df) || !nrow(df) || !label_col %in% colnames(df)) return(df)
    df <- df[!is.na(df[[label_col]]) & nzchar(as.character(df[[label_col]])), , drop = FALSE]
    sample_cols <- setdiff(colnames(df), unique(c(label_col, keep_cols, "gene_id", "Geneid", "transcript_id", "target_id", "gene_symbol", "transcript_name", "biotype")))
    numeric_cols <- sample_cols[vapply(df[, sample_cols, drop = FALSE], function(x) any(!is.na(suppressWarnings(as.numeric(x)))), logical(1))]
    if (!length(numeric_cols)) return(df)
    for (col in numeric_cols) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    agg <- aggregate(df[, numeric_cols, drop = FALSE], by = list(df[[label_col]]), FUN = sum, na.rm = TRUE)
    colnames(agg)[1] <- label_col
    agg[order(agg[[label_col]]), , drop = FALSE]
  }

  aggregated_gene_name_matrix_path <- function(path) {
    sub("(\\.txt|\\.tsv)$", "_aggregated.txt", gene_name_matrix_path(path))
  }

  save_aggregated_gene_name_matrix <- function(path, mode = c("rsem", "kallisto", "featurecounts", "deseq")) {
    mode <- match.arg(mode)
    gene_path <- gene_name_matrix_path(path)
    if (!file.exists(gene_path)) {
      saved <- save_gene_name_matrix(path, mode)
      if (!isTRUE(saved$ok)) return(saved)
      gene_path <- saved$path
    }
    df <- read_saved_count_matrix(gene_path)
    if (is.null(df) || !nrow(df) || !"gene_name" %in% colnames(df)) {
      return(list(ok = FALSE, message = "Gene-name matrix was not found or did not contain a gene_name column.", path = NA_character_))
    }
    out <- aggregate_gene_name_df(df, "gene_name")
    out_file <- aggregated_gene_name_matrix_path(path)
    write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    list(ok = TRUE, message = sprintf("Saved duplicate-combined gene-name matrix with %s genes to %s.", nrow(out), out_file), path = out_file)
  }

  save_gene_name_matrix <- function(path, mode = c("rsem", "kallisto", "featurecounts", "deseq")) {
    mode <- match.arg(mode)
    df <- read_saved_count_matrix(path)
    if (is.null(df) || !nrow(df)) {
      return(list(ok = FALSE, message = "Matrix file was not found or was empty.", path = NA_character_))
    }
    if (grepl("_gene_name\\.txt$", basename(path))) {
      return(list(ok = TRUE, message = "This matrix is already the gene-name version.", path = path))
    }
    out_file <- gene_name_matrix_path(path)
    if (identical(mode, "kallisto")) {
      id_col <- intersect(c("target_id", "transcript_id", "gene_id"), colnames(df))[1]
      if (is.na(id_col)) id_col <- colnames(df)[1]
      meta <- parse_kallisto_target_metadata(df[[id_col]])
      if ("gene_symbol" %in% colnames(df)) meta$gene_name <- df$gene_symbol
      if ("gene_name" %in% colnames(df)) meta$gene_name <- df$gene_name
      if ("gene_id" %in% colnames(df)) meta$gene_id <- df$gene_id
      if ("transcript_id" %in% colnames(df)) meta$transcript_id <- df$transcript_id
      if ("transcript_name" %in% colnames(df)) meta$transcript_name <- df$transcript_name
      if ("biotype" %in% colnames(df)) meta$biotype <- df$biotype
      if (!any(!is.na(meta$gene_name) & nzchar(meta$gene_name)) && any(!is.na(meta$gene_id) & nzchar(meta$gene_id))) {
        species <- detect_species_from_ids(meta$gene_id)
        map_df <- get_gtf_map(species)
        conv <- convert_gene_labels(meta$gene_id, map_df)
        meta$gene_name <- conv$values
      }
      sample_cols <- setdiff(colnames(df), c(id_col, "target_id", "transcript_id", "gene_id", "gene_name", "gene_symbol", "transcript_name", "biotype"))
      keep_meta <- c("gene_name", "transcript_id", "gene_id", "transcript_name", "biotype")
      out <- cbind(meta[, keep_meta, drop = FALSE], df[, sample_cols, drop = FALSE])
      mapped <- sum(!is.na(out$gene_name) & nzchar(out$gene_name) & out$gene_name != out$gene_id, na.rm = TRUE)
      write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      return(list(ok = TRUE, message = sprintf("Saved Kallisto gene-name matrix (%s rows with gene names) to %s.", mapped, out_file), path = out_file))
    }

    gene_col <- intersect(c("gene_id", "Geneid"), colnames(df))[1]
    if (is.na(gene_col)) gene_col <- colnames(df)[1]
    ids <- as.character(df[[gene_col]])
    species <- detect_species_from_ids(ids)
    map_df <- get_gtf_map(species)
    conv <- convert_gene_labels(ids, map_df)
    mapped <- sum(!is.na(conv$values) & nzchar(conv$values) & conv$values != ids)
    if (mapped == 0) {
      return(list(ok = FALSE, message = "No Ensembl gene IDs were mapped to gene names from the local GTF.", path = NA_character_))
    }
    if ("transcript_id" %in% colnames(df)) {
      df$gene_name <- conv$values
      out <- df[, c("transcript_id", "gene_name", setdiff(colnames(df), c("transcript_id", "gene_name", gene_col))), drop = FALSE]
    } else {
      df[[gene_col]] <- conv$values
      colnames(df)[match(gene_col, colnames(df))] <- "gene_name"
      out <- df
    }
    write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    label <- switch(mode, rsem = "RSEM", kallisto = "Kallisto", featurecounts = "featureCounts", deseq = "DESeq2")
    list(ok = TRUE, message = sprintf("Saved %s gene-name matrix (%s IDs mapped) to %s.", label, mapped, out_file), path = out_file)
  }

  current_resource_prefix <- reactive({
    if (isTRUE(input$show_trimmed)) "fastqc_cutadapt_results" else "fastqc_results"
  })

  qc_mode_label <- reactive({
    if (isTRUE(input$show_trimmed)) "trimmed" else "untrimmed"
  })

  qc_mode_available <- reactive({
    if (isTRUE(input$show_trimmed)) fastqc_trim_available else fastqc_raw_available
  })

  current_files <- reactive({
    if (is.null(input$sample) || !nzchar(input$sample)) return(NULL)
    list(
      r1 = qc_file_map(input$sample, "R1", trimmed = isTRUE(input$show_trimmed)),
      r2 = qc_file_map(input$sample, "R2", trimmed = isTRUE(input$show_trimmed))
    )
  })

  output$qc_status_ui <- renderUI({
    if (!fastqc_raw_available && !fastqc_trim_available) {
      return(status_box("FastQC has not been run yet.", "warning"))
    }
    if (!qc_mode_available()) {
      return(status_box(sprintf("QC has not been run on %s reads.", qc_mode_label()), "warning"))
    }
    if (is.null(input$sample) || !nzchar(input$sample)) {
      return(status_box("No sample names are available yet for QC display.", "warning"))
    }
    NULL
  })

  output$r1_fastqc_ui <- renderUI({
    if (is.null(current_files())) return(tags$p("No sample selected."))
    if (!qc_mode_available()) return(tags$p(sprintf("QC has not been run on %s reads.", qc_mode_label())))
    iframe_or_message(current_files()$r1$fastqc, current_resource_prefix(), missing_message = sprintf("FastQC has not been generated yet for %s R1 (%s reads).", input$sample, qc_mode_label()))
  })

  output$r1_screen_ui <- renderUI({
    if (is.null(current_files())) return(tags$p("No sample selected."))
    if (!qc_mode_available()) return(tags$p(sprintf("QC has not been run on %s reads.", qc_mode_label())))
    iframe_or_message(current_files()$r1$screen, current_resource_prefix(), missing_message = sprintf("FastQ Screen has not been generated yet for %s R1 (%s reads).", input$sample, qc_mode_label()))
  })

  output$r2_fastqc_ui <- renderUI({
    if (is.null(current_files())) return(tags$p("No sample selected."))
    if (!qc_mode_available()) return(tags$p(sprintf("QC has not been run on %s reads.", qc_mode_label())))
    iframe_or_message(current_files()$r2$fastqc, current_resource_prefix(), missing_message = sprintf("FastQC has not been generated yet for %s R2 (%s reads).", input$sample, qc_mode_label()))
  })

  output$r2_screen_ui <- renderUI({
    if (is.null(current_files())) return(tags$p("No sample selected."))
    if (!qc_mode_available()) return(tags$p(sprintf("QC has not been run on %s reads.", qc_mode_label())))
    iframe_or_message(current_files()$r2$screen, current_resource_prefix(), missing_message = sprintf("FastQ Screen has not been generated yet for %s R2 (%s reads).", input$sample, qc_mode_label()))
  })

  output$star_status_ui <- renderUI({
    if (!star_available) {
      return(status_box("STAR alignment summary has not been generated yet.", "warning"))
    }
    NULL
  })

  star_summary_table_df <- reactive({
    req(star_available)
    order_sample_columns(star_summary_df, annotation_col = value_or(input$star_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("metric"))
  })

  star_sample_table_df <- reactive({
    req(star_available)
    req(input$star_sample)
    sample_col <- input$star_sample
    req(sample_col %in% colnames(star_summary_df))
    data.frame(
      Metric = star_summary_df$metric,
      Value = star_summary_df[[sample_col]],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  if (DT_AVAILABLE) {
    output$star_summary_table <- DT::renderDT({
      df <- star_summary_table_df()
      simple_dt(df, page_length = 50)
    }, server = FALSE)
  } else {
    output$star_summary_table <- renderTable({
      df <- star_summary_table_df()
      rownames(df) <- df$metric
      df$metric <- NULL
      format_numeric_commas(df)
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
  }

  output$star_sample_table <- renderTable({
    format_numeric_commas(star_sample_table_df())
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$download_star_summary <- downloadHandler(
    filename = function() download_filename("qc_star_summary"),
    content = function(file) write_csv_download(safe_download_df(star_summary_table_df()), file)
  )
  output$download_star_sample <- downloadHandler(
    filename = function() download_filename(paste0("qc_star_", value_or(input$star_sample, "sample"))),
    content = function(file) write_csv_download(safe_download_df(star_sample_table_df()), file)
  )

  output$featurecounts_qc_status_ui <- renderUI({
    if (!featurecounts_available) {
      return(status_box("FeatureCounts summary has not been generated yet.", "warning"))
    }
    NULL
  })

  featurecounts_summary_table_df <- reactive({
    req(featurecounts_available)
    order_sample_columns(featurecounts_summary_df, annotation_col = value_or(input$featurecounts_qc_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("Status"))
  })

  featurecounts_sample_table_df <- reactive({
    req(featurecounts_available)
    req(input$featurecounts_qc_sample)
    sample_col <- input$featurecounts_qc_sample
    req(sample_col %in% colnames(featurecounts_summary_df))
    data.frame(
      Status = featurecounts_summary_df$Status,
      Value = featurecounts_summary_df[[sample_col]],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  if (DT_AVAILABLE) {
    output$featurecounts_summary_table <- DT::renderDT({
      df <- featurecounts_summary_table_df()
      simple_dt(df, page_length = 50)
    }, server = FALSE)
  } else {
    output$featurecounts_summary_table <- renderTable({
      df <- featurecounts_summary_table_df()
      rownames(df) <- df$Status
      df$Status <- NULL
      format_numeric_commas(df)
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
  }

  output$featurecounts_sample_table <- renderTable({
    format_numeric_commas(featurecounts_sample_table_df())
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$download_featurecounts_summary <- downloadHandler(
    filename = function() download_filename("qc_featurecounts_summary"),
    content = function(file) write_csv_download(safe_download_df(featurecounts_summary_table_df()), file)
  )
  output$download_featurecounts_sample <- downloadHandler(
    filename = function() download_filename(paste0("qc_featurecounts_", value_or(input$featurecounts_qc_sample, "sample"))),
    content = function(file) write_csv_download(safe_download_df(featurecounts_sample_table_df()), file)
  )

  output$featurecounts_convert_ui <- renderUI({
    if (!count_matrix_available) return(NULL)
    if (identical(value_or(input$featurecounts_label_mode, "gene_id"), "gene_name")) {
      converted_path <- gene_name_matrix_path(count_matrix_path)
      aggregated_path <- aggregated_gene_name_matrix_path(count_matrix_path)
      if (file.exists(aggregated_path)) {
        tags$span(class = "tiny-note", "Using saved duplicate-combined gene-name count matrix.")
      } else if (file.exists(converted_path)) {
        tags$span(class = "tiny-note", "Using saved gene-name count matrix.")
      } else {
        tags$span(class = "tiny-note", "Gene-name count matrix will be saved automatically.")
      }
    }
  })

  featurecounts_display_df <- reactive({
    if (identical(value_or(input$featurecounts_label_mode, "gene_id"), "gene_name")) {
      aggregated_path <- aggregated_gene_name_matrix_path(count_matrix_path)
      if (!file.exists(aggregated_path)) {
        save_aggregated_gene_name_matrix(count_matrix_path, "featurecounts")
      }
      df <- read_saved_count_matrix(aggregated_path)
      if (!is.null(df) && "gene_name" %in% colnames(df)) {
        colnames(df)[colnames(df) == "gene_name"] <- "Geneid"
        return(df)
      }
      converted_path <- gene_name_matrix_path(count_matrix_path)
      df <- read_saved_count_matrix(converted_path)
      if (!is.null(df) && "gene_name" %in% colnames(df)) {
        colnames(df)[colnames(df) == "gene_name"] <- "Geneid"
        return(df)
      }
    }
    count_matrix_nonzero_df
  })

  output$featurecounts_status_ui <- renderUI({
    if (!count_matrix_available) {
      return(status_box("Raw count matrix has not been generated yet.", "warning"))
    }
    if (identical(value_or(input$featurecounts_label_mode, "gene_id"), "gene_name")) {
      aggregated_path <- aggregated_gene_name_matrix_path(count_matrix_path)
      converted_path <- gene_name_matrix_path(count_matrix_path)
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        if (file.exists(aggregated_path)) {
          sprintf("Showing saved duplicate-combined gene-name counts: %s", aggregated_path)
        } else if (file.exists(converted_path)) {
          sprintf("Showing saved gene-name counts: %s", converted_path)
        } else {
          "Gene-name counts could not be saved from the current raw count matrix."
        }
      )
    } else {
      NULL
    }
  })

  observe({
    req(identical(input$main_tabs, "Counts"))
    req(identical(input$counts_subtab, "Raw Counts"))
    if (!count_matrix_available) {
      updateSelectizeInput(session, "gene_query", choices = character(0), selected = character(0), server = TRUE)
      return()
    }
    df <- featurecounts_display_df()
    updateSelectizeInput(
      session,
      "gene_query",
      choices = df$Geneid,
      selected = character(0),
      server = TRUE
    )
  })

  raw_counts_table_df <- reactive({
    req(count_matrix_available)
    display_df <- featurecounts_display_df()
    q <- input$gene_query
    if (is.null(q) || !nzchar(trimws(q))) {
      show_df <- display_df
    } else {
      idx <- tolower(display_df$Geneid) == tolower(trimws(q))
      hits <- display_df[idx, , drop = FALSE]
      show_df <- if (nrow(hits) == 0) display_df else hits
    }
    sort_df_by_col(show_df, "Geneid", "asc")
  })

  if (DT_AVAILABLE) {
    output$gene_search_table <- DT::renderDT({
      df <- raw_counts_table_df()
      simple_dt(df, page_length = 50)
    }, server = FALSE)
  } else {
    output$gene_search_table <- renderTable({
      df <- raw_counts_table_df()
      rownames(df) <- make.unique(as.character(df$Geneid))
      df$Geneid <- NULL
      format_numeric_commas(df)
    }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
  }

  output$download_raw_counts_table <- downloadHandler(
    filename = function() download_filename("counts_raw_featurecounts"),
    content = function(file) write_csv_download(safe_download_df(raw_counts_table_df()), file)
  )

  if (rsem_available) {
    rsem_matrix_cache <- reactiveValues()
    rsem_display_cache <- reactiveValues()
    rsem_label_status <- reactiveVal(NULL)
    rsem_conversion_status <- reactiveVal(NULL)
    rsem_loaded <- reactiveVal(FALSE)

    observeEvent(input$rsem_load_matrix, {
      rsem_loaded(TRUE)
      rsem_label_status(NULL)
    }, ignoreInit = TRUE)

    observeEvent(list(input$rsem_type, input$rsem_metric, input$rsem_label_mode), {
      rsem_loaded(FALSE)
      rsem_label_status(NULL)
      rsem_conversion_status(NULL)
      updateSelectizeInput(session, "rsem_query", choices = character(0), selected = character(0), server = TRUE)
    }, ignoreInit = TRUE)

    output$rsem_convert_ui <- renderUI({
      req(input$rsem_type, input$rsem_metric)
      base_path <- rsem_saved_matrix_path(input$rsem_type, input$rsem_metric, "gene_id")
      converted_path <- gene_name_matrix_path(base_path)
      if (!file.exists(base_path)) {
        return(tags$span(class = "tiny-note", "Saved RSEM matrix not found yet. Run the RSEM step first."))
      }
      if (identical(value_or(input$rsem_label_mode, "gene_id"), "gene_name")) {
        if (!file.exists(converted_path)) save_gene_name_matrix(base_path, "rsem")
        return(tags$span(class = "tiny-note", if (file.exists(converted_path)) "Using saved RSEM gene-name matrix." else "RSEM gene-name matrix could not be saved."))
      }
      NULL
    })

    get_rsem_matrix <- reactive({
      req(input$rsem_type, input$rsem_metric)
      label_mode <- value_or(input$rsem_label_mode, "gene_id")
      if (identical(label_mode, "gene_name")) {
        base_path <- rsem_saved_matrix_path(input$rsem_type, input$rsem_metric, "gene_id")
        converted_path <- gene_name_matrix_path(base_path)
        if (!file.exists(converted_path)) {
          res <- save_gene_name_matrix(base_path, "rsem")
          rsem_conversion_status(res$message)
        }
      }
      matrix_path <- rsem_saved_matrix_path(input$rsem_type, input$rsem_metric, label_mode)
      cache_key <- matrix_path
      cached <- rsem_matrix_cache[[cache_key]]
      if (is.null(cached)) {
        cached <- read_saved_count_matrix(matrix_path)
        rsem_matrix_cache[[cache_key]] <- cached
      }
      cached
    })

    rsem_display_df <- reactive({
      req(isTRUE(rsem_loaded()))
      raw_df <- get_rsem_matrix()
      req(!is.null(raw_df))
      label_mode <- value_or(input$rsem_label_mode, "gene_id")
      matrix_path <- rsem_saved_matrix_path(input$rsem_type, input$rsem_metric, label_mode)
      cache_key <- paste(matrix_path, label_mode, sep = "::")
      cached <- rsem_display_cache[[cache_key]]
      if (!is.null(cached)) {
        rsem_label_status(attr(cached, "label_status", exact = TRUE))
        return(cached)
      }

      df <- raw_df
      status_message <- rsem_conversion_status()
      if (is.null(status_message) || !nzchar(status_message)) {
        status_message <- sprintf("Using saved RSEM matrix: %s", matrix_path)
      }
      if (!grepl("_gene_name\\.txt$", basename(matrix_path)) && "gene_id" %in% colnames(df)) {
        original_gene_id <- as.character(df$gene_id)
        if (identical(label_mode, "gene_name")) {
          if (looks_like_gene_id(original_gene_id)) {
            species <- detect_species_from_ids(original_gene_id)
            map_df <- get_gtf_map(species)
            conv <- convert_gene_labels(original_gene_id, map_df)
            df$gene_name <- conv$values
            mapped_all <- sum(!is.na(df$gene_name) & nzchar(df$gene_name) & df$gene_name != original_gene_id)
            status_message <- if (mapped_all > 0) {
              sprintf("RSEM display is showing gene names using %s GTF (%s IDs mapped).", value_or(species, "detected"), mapped_all)
            } else {
              sprintf("RSEM gene-name display requested, but no IDs were mapped using the %s GTF.", value_or(species, "detected"))
            }
          } else {
            df$gene_name <- original_gene_id
            status_message <- "RSEM gene labels already look like gene names."
          }
        }
      }

      if (identical(label_mode, "gene_name") && "gene_name" %in% colnames(df)) {
        if (identical(input$rsem_type, "genes")) {
          df <- df[, c("gene_name", setdiff(colnames(df), c("gene_name", "gene_id"))), drop = FALSE]
        } else {
          df <- df[, c("transcript_id", "gene_name", setdiff(colnames(df), c("transcript_id", "gene_name", "gene_id"))), drop = FALSE]
        }
      } else if (identical(input$rsem_type, "genes") && "gene_id" %in% colnames(df)) {
        df <- df[, c("gene_id", setdiff(colnames(df), c("gene_id", "gene_name"))), drop = FALSE]
      } else {
        df <- df[, c(intersect(c("transcript_id", "gene_id"), colnames(df)), setdiff(colnames(df), c("transcript_id", "gene_id", "gene_name"))), drop = FALSE]
      }

      attr(df, "label_status") <- status_message
      rsem_display_cache[[cache_key]] <- df
      rsem_label_status(status_message)
      df
    })

    output$rsem_status_ui <- renderUI({
      if (!isTRUE(rsem_loaded())) {
        return(status_box("Choose RSEM settings, then click 'Load RSEM matrix' to render the table.", "info"))
      }
      df <- rsem_display_df()
      if (is.null(df)) {
        return(status_box("Saved RSEM matrix has not been generated yet. Run the RSEM step to build it in the counts folder.", "warning"))
      }
      msg <- rsem_label_status()
      if (is.null(msg) || !nzchar(msg)) return(NULL)
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        msg
      )
    })

    observe({
      req(identical(input$main_tabs, "Counts"))
      req(identical(input$counts_subtab, "RSEM"))
      req(input$rsem_type, input$rsem_metric)
      req(isTRUE(rsem_loaded()))
      df <- rsem_display_df()
      req(!is.null(df))
      id_col <- if (identical(input$rsem_type, "genes")) colnames(df)[1] else if (identical(value_or(input$rsem_label_mode, "gene_id"), "gene_name") && "gene_name" %in% colnames(df)) "gene_name" else "transcript_id"
      vals <- sort(unique(df[[id_col]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      updateSelectizeInput(
        session,
        "rsem_query",
        choices = vals,
        selected = character(0),
        server = TRUE
      )
      updateSelectInput(
        session,
        "rsem_sort_col",
        choices = colnames(df),
        selected = if (id_col %in% colnames(df)) id_col else colnames(df)[1]
      )
    })

    rsem_table_df <- reactive({
      df <- rsem_display_df()
      req(!is.null(df))
      id_col <- if (identical(input$rsem_type, "genes")) colnames(df)[1] else if (identical(value_or(input$rsem_label_mode, "gene_id"), "gene_name") && "gene_name" %in% colnames(df)) "gene_name" else "transcript_id"
      q <- input$rsem_query
      if (is.null(q) || !nzchar(trimws(q))) {
        show_df <- df
      } else {
        hits <- df[tolower(df[[id_col]]) == tolower(trimws(q)), , drop = FALSE]
        show_df <- if (nrow(hits) == 0) df else hits
      }
      sort_df_by_col(show_df, id_col, "asc")
    })

    if (DT_AVAILABLE) {
      output$rsem_table <- DT::renderDT({
        if (!isTRUE(rsem_loaded())) {
          return(simple_dt(data.frame(Message = "Click 'Load RSEM matrix' to render this table.", stringsAsFactors = FALSE), page_length = 50, scroll_y = "140px"))
        }
        df <- rsem_table_df()
        simple_dt(df, page_length = 50)
      }, server = FALSE)
    } else {
      output$rsem_table <- renderTable({
        if (!isTRUE(rsem_loaded())) return(data.frame(Message = "Click 'Load RSEM matrix' to render this table.", stringsAsFactors = FALSE))
        df <- rsem_table_df()
        id_col <- if (identical(input$rsem_type, "genes")) colnames(df)[1] else if (identical(value_or(input$rsem_label_mode, "gene_id"), "gene_name") && "gene_name" %in% colnames(df)) "gene_name" else "transcript_id"
        rownames(df) <- make.unique(as.character(df[[id_col]]))
        df[[id_col]] <- NULL
        format_numeric_commas(df)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
    }
    output$download_rsem_table <- downloadHandler(
      filename = function() download_filename(paste0("counts_rsem_", value_or(input$rsem_type, "matrix"), "_", value_or(input$rsem_metric, "metric"), "_", value_or(input$rsem_label_mode, "label"))),
      content = function(file) write_csv_download(safe_download_df(if (isTRUE(rsem_loaded())) rsem_table_df() else NULL), file)
    )
  } else {
    output$rsem_convert_ui <- renderUI(NULL)
    output$rsem_status_ui <- renderUI({
      status_box("RSEM results have not been generated yet.", "warning")
    })
    output$rsem_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable(NULL)
    output$download_rsem_table <- downloadHandler(
      filename = function() download_filename("counts_rsem_unavailable"),
      content = function(file) write_csv_download(safe_download_df(NULL), file)
    )
  }

  if (deseq2_available) {
    get_deseq_values <- reactive({
      if (is.null(input$deseq_compare_col) || identical(input$deseq_compare_col, "NA") || !nzchar(input$deseq_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$deseq_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$deseq_treatment_ui <- renderUI({
      vals <- get_deseq_values()
      selectInput("deseq_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals, selectize = FALSE)
    })

    output$deseq_control_ui <- renderUI({
      vals <- get_deseq_values()
      selectInput("deseq_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals, selectize = FALSE)
    })

    deseq_counts_raw_df <- reactive({
      req(input$deseq_compare_col)
      if (identical(input$deseq_compare_col, "NA")) return(NULL)
      req(input$deseq_treatment, input$deseq_control)
      path <- normalized_counts_path(input$deseq_treatment, input$deseq_control, "gene_id")
      if (identical(value_or(input$deseq_label_mode, "gene_id"), "gene_name")) {
        gene_name_run_path <- normalized_counts_path(input$deseq_treatment, input$deseq_control, "gene_name")
        if (file.exists(gene_name_run_path)) {
          path <- gene_name_run_path
        } else {
          saved <- save_normalized_counts_gene_names(input$deseq_treatment, input$deseq_control)
          if (isTRUE(saved$ok) && file.exists(saved$path)) path <- saved$path
        }
      }
      df <- read_normalized_counts(path)
      if (is.null(df)) return(NULL)
      if ("gene_name" %in% colnames(df) && !"gene_label" %in% colnames(df)) colnames(df)[colnames(df) == "gene_name"] <- "gene_label"
      keep_samples <- design_df$sample[design_df[[input$deseq_compare_col]] %in% c(input$deseq_treatment, input$deseq_control)]
      keep_cols <- c("gene_label", intersect(keep_samples, colnames(df)))
      df[, keep_cols, drop = FALSE]
    })

    output$deseq_convert_ui <- renderUI({
      df <- deseq_counts_raw_df()
      if (is.null(df)) return(NULL)
      if (identical(value_or(input$deseq_label_mode, "gene_id"), "gene_name")) {
        path <- normalized_counts_gene_name_path(input$deseq_treatment, input$deseq_control)
        if (file.exists(path)) {
          tags$span(class = "tiny-note", "Using saved duplicate-combined gene-name normalized counts.")
        } else {
          tags$span(class = "tiny-note", "Using gene-name DESeq2 results when available.")
        }
      }
    })

    deseq_counts_df <- reactive({
      df <- deseq_counts_raw_df()
      if (is.null(df)) return(NULL)
      df
    })

    output$deseq_status_ui <- renderUI({
      if (!deseq2_available) {
        return(status_box("DESeq2 normalized counts have not been generated yet.", "warning"))
      }
      if (!nrow(design_df) || !length(comparison_columns)) {
        return(status_box("The design matrix is not available yet, so comparisons cannot be selected.", "warning"))
      }
      req(input$deseq_compare_col)
      if (identical(input$deseq_compare_col, "NA")) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fff7e6; border: 1px solid #f0c36d; border-radius: 6px;",
          "Select a comparison column to view normalized counts."
        ))
      }
      req(input$deseq_treatment, input$deseq_control)
      path <- normalized_counts_path(input$deseq_treatment, input$deseq_control, "gene_id")
      if (identical(value_or(input$deseq_label_mode, "gene_id"), "gene_name")) {
        gene_name_run_path <- normalized_counts_path(input$deseq_treatment, input$deseq_control, "gene_name")
        saved_path <- normalized_counts_gene_name_path(input$deseq_treatment, input$deseq_control)
        if (file.exists(gene_name_run_path)) path <- gene_name_run_path else if (file.exists(saved_path)) path <- saved_path
      }
      if (!file.exists(path)) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;",
          sprintf("This comparison has not been tested: %s vs %s.", input$deseq_treatment, input$deseq_control)
        ))
      }
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf(
          "Showing normalized counts for %s vs %s using %s (%s).",
          input$deseq_treatment,
          input$deseq_control,
          input$deseq_compare_col,
          if (identical(value_or(input$deseq_label_mode, "gene_id"), "gene_name")) "gene names" else "gene IDs"
        )
      )
    })

    observe({
      req(identical(input$main_tabs, "Counts"))
      req(identical(input$counts_subtab, "DESeq2 Normalized Counts"))
      df <- deseq_counts_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "deseq_gene_query", choices = character(0), selected = character(0), server = TRUE)
        updateSelectInput(session, "deseq_sort_col", choices = c("gene_label"), selected = "gene_label")
      } else {
        updateSelectizeInput(session, "deseq_gene_query", choices = df$gene_label, selected = character(0), server = TRUE)
        updateSelectInput(
          session,
          "deseq_sort_col",
          choices = colnames(df),
          selected = if ("gene_label" %in% colnames(df)) "gene_label" else colnames(df)[1]
        )
      }
    })

    deseq_counts_table_df <- reactive({
      df <- deseq_counts_df()
      req(!is.null(df))
      q <- input$deseq_gene_query
      if (is.null(q) || !nzchar(trimws(q))) {
        show_df <- df
      } else {
        hits <- df[tolower(df$gene_label) == tolower(trimws(q)), , drop = FALSE]
        show_df <- if (nrow(hits) == 0) df else hits
      }
      sort_df_by_col(show_df, "gene_label", "asc")
    })

    if (DT_AVAILABLE) {
      output$deseq_counts_table <- DT::renderDT({
        df <- deseq_counts_table_df()
        simple_dt(df, page_length = 50)
      }, server = FALSE)
    } else {
      output$deseq_counts_table <- renderTable({
        show_df <- deseq_counts_table_df()
        rownames(show_df) <- make.unique(as.character(show_df$gene_label))
        show_df$gene_label <- NULL
        format_numeric_commas(show_df)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
    }
    output$download_deseq_counts_table <- downloadHandler(
      filename = function() download_filename(paste0("counts_deseq2_normalized_", value_or(input$deseq_treatment, "treatment"), "_vs_", value_or(input$deseq_control, "control"))),
      content = function(file) write_csv_download(safe_download_df(deseq_counts_table_df()), file)
    )

    get_deg_values <- reactive({
      if (is.null(input$deg_compare_col) || identical(input$deg_compare_col, "NA") || !nzchar(input$deg_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$deg_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$deg_treatment_ui <- renderUI({
      vals <- get_deg_values()
      selectInput("deg_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals, selectize = FALSE)
    })

    output$deg_control_ui <- renderUI({
      vals <- get_deg_values()
      selectInput("deg_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals, selectize = FALSE)
    })

    output$deg_label_mode_ui <- renderUI({
      choices <- available_deseq2_label_modes(value_or(input$deg_treatment, ""), value_or(input$deg_control, ""))
      if (!length(choices)) choices <- c("Gene ID" = "gene_id")
      selected <- value_or(input$deg_label_mode, unname(choices)[1])
      if (!selected %in% unname(choices)) selected <- unname(choices)[1]
      selectInput("deg_label_mode", "DESeq2 result labels", choices = choices, selected = selected, selectize = FALSE)
    })

    get_plot_values <- reactive({
      if (is.null(input$plot_compare_col) || identical(input$plot_compare_col, "NA") || !nzchar(input$plot_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$plot_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$plot_treatment_ui <- renderUI({
      vals <- get_plot_values()
      selectInput("plot_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals, selectize = FALSE)
    })

    output$plot_control_ui <- renderUI({
      vals <- get_plot_values()
      selectInput("plot_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals, selectize = FALSE)
    })

    output$plot_deseq_label_mode_ui <- renderUI({
      treatment_value <- value_or(input$plot_treatment, "")
      control_value <- value_or(input$plot_control, "")
      choices <- available_deseq2_label_modes(treatment_value, control_value)
      if (!length(choices)) choices <- available_deseq2_label_modes()
      if (!length(choices)) choices <- c("Gene ID" = "gene_id")
      selected <- value_or(input$plot_deseq_label_mode, unname(choices)[1])
      if (!selected %in% unname(choices)) selected <- unname(choices)[1]
      selectInput(
        "plot_deseq_label_mode",
        "DESeq2 result labels",
        choices = choices,
        selected = selected,
        selectize = FALSE
      )
    })

    deg_raw_df <- reactive({
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      req(input$deg_treatment, input$deg_control)
      read_deg_table(deg_table_path(input$deg_treatment, input$deg_control, value_or(input$deg_label_mode, "gene_id")))
    })

    output$deg_convert_ui <- renderUI({
      NULL
    })

    deg_df <- reactive({
      df <- deg_raw_df()
      if (is.null(df)) return(NULL)
      df
    })

    plot_raw_df <- reactive({
      req(input$plot_compare_col)
      if (identical(input$plot_compare_col, "NA")) return(NULL)
      req(input$plot_treatment, input$plot_control)
      read_deg_table(deg_table_path(input$plot_treatment, input$plot_control, value_or(input$plot_deseq_label_mode, "gene_id")))
    })

    plot_deg_df <- reactive({
      df <- plot_raw_df()
      if (is.null(df)) return(NULL)
      df
    })

    output$deg_status_ui <- renderUI({
      if (!deseq2_available) {
        return(status_box("Differential expression results have not been generated yet.", "warning"))
      }
      if (!nrow(design_df) || !length(comparison_columns)) {
        return(status_box("The design matrix is not available yet, so comparisons cannot be selected.", "warning"))
      }
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fff7e6; border: 1px solid #f0c36d; border-radius: 6px;",
          "Select a comparison column to view differential expression results."
        ))
      }
      req(input$deg_treatment, input$deg_control)
      path <- deg_table_path(input$deg_treatment, input$deg_control, value_or(input$deg_label_mode, "gene_id"))
      if (!file.exists(path)) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;",
          sprintf("This comparison has not been tested: %s vs %s.", input$deg_treatment, input$deg_control)
        ))
      }
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf(
          "Showing DEG results for %s vs %s using %s (%s).",
          input$deg_treatment,
          input$deg_control,
          input$deg_compare_col,
          if (identical(value_or(input$deg_label_mode, "gene_id"), "gene_name")) "gene names" else "gene IDs"
        )
      )
    })

    output$plot_status_ui <- renderUI({
      if (!deseq2_available) {
        return(status_box("Differential expression outputs have not been generated yet.", "warning"))
      }
      if (identical(value_or(input$plot_subtab, "PCA"), "PCA")) {
        return(status_box("PCA uses all selected samples by default. Use the metadata controls to color and optionally filter samples.", "info"))
      }
      req(input$plot_compare_col)
      if (identical(input$plot_compare_col, "NA")) {
        return(status_box("Select a comparison column to view plots for a tested comparison.", "warning"))
      }
      req(input$plot_treatment, input$plot_control)
      path <- deg_table_path(input$plot_treatment, input$plot_control, value_or(input$plot_deseq_label_mode, "gene_id"))
      if (!file.exists(path)) {
        return(status_box(sprintf("This comparison has not been tested: %s vs %s.", input$plot_treatment, input$plot_control), "error"))
      }
      label_text <- if (identical(value_or(input$plot_deseq_label_mode, "gene_id"), "gene_name")) "gene names" else "gene IDs"
      status_box(sprintf("Plotting %s vs %s using %s (%s).", input$plot_treatment, input$plot_control, input$plot_compare_col, label_text), "info")
    })

    deg_norm_counts_df <- reactive({
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      req(input$deg_treatment, input$deg_control)
      path <- normalized_counts_path(input$deg_treatment, input$deg_control)
      df <- read_normalized_counts(path)
      if (is.null(df)) return(NULL)
      keep_samples <- design_df$sample[design_df[[input$deg_compare_col]] %in% c(input$deg_treatment, input$deg_control)]
      keep_cols <- c("gene_label", intersect(keep_samples, colnames(df)))
      df[, keep_cols, drop = FALSE]
    })

    plot_norm_counts_df <- reactive({
      req(input$plot_compare_col)
      if (identical(input$plot_compare_col, "NA")) return(NULL)
      req(input$plot_treatment, input$plot_control)
      path <- normalized_counts_path(input$plot_treatment, input$plot_control, value_or(input$plot_deseq_label_mode, "gene_id"))
      df <- read_normalized_counts(path)
      if (is.null(df)) return(NULL)
      keep_samples <- design_df$sample[design_df[[input$plot_compare_col]] %in% c(input$plot_treatment, input$plot_control)]
      keep_cols <- c("gene_label", intersect(keep_samples, colnames(df)))
      df[, keep_cols, drop = FALSE]
    })

    heatmap_rsem_df <- reactive({
      req(rsem_available)
      req(input$plot_compare_col)
      if (identical(input$plot_compare_col, "NA")) return(NULL)
      req(input$plot_treatment, input$plot_control)
      df <- read_saved_count_matrix(rsem_saved_matrix_path("genes", "TPM", "gene_id"))
      if (is.null(df)) return(NULL)
      keep_samples <- design_df$sample[design_df[[input$plot_compare_col]] %in% c(input$plot_treatment, input$plot_control)]
      keep_cols <- c("gene_id", intersect(keep_samples, colnames(df)))
      df <- df[, keep_cols, drop = FALSE]
      colnames(df)[colnames(df) == "gene_id"] <- "gene_label"
      df
    })

    heatmap_source_df <- reactive({
      if (identical(value_or(input$heatmap_source, "deseq"), "rsem")) {
        heatmap_rsem_df()
      } else {
        plot_norm_counts_df()
      }
    })

    output$pca_include_values_ui <- renderUI({
      col <- value_or(input$pca_color_col, "__none__")
      if (identical(col, "__none__") || !col %in% colnames(design_df)) return(NULL)
      vals <- sort(unique(as.character(design_df[[col]])))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      selectizeInput("pca_include_values", "Include values", choices = vals, selected = vals, multiple = TRUE, options = list(dropdownParent = "body"))
    })

    pca_source_df <- reactive({
      if (identical(value_or(input$pca_source, "deseq"), "counts")) {
        df <- count_matrix_nonzero_df
        if (is.null(df) || !nrow(df)) return(NULL)
        colnames(df)[1] <- "gene_label"
        return(df)
      }
      selected_path <- NULL
      if (!is.null(input$plot_treatment) && !is.null(input$plot_control)) {
        candidate <- normalized_counts_path(input$plot_treatment, input$plot_control, value_or(input$plot_deseq_label_mode, "gene_id"))
        if (file.exists(candidate)) selected_path <- candidate
      }
      if (is.null(selected_path)) selected_path <- first_normalized_counts_path(value_or(input$plot_deseq_label_mode, "gene_id"))
      if (is.null(selected_path)) selected_path <- first_normalized_counts_path()
      read_normalized_counts(selected_path)
    })

    output$pca_plot <- renderPlot({
      df <- pca_source_df()
      validate(need(!is.null(df), "No count matrix or DESeq2 normalized-counts file is available for PCA."))
      mat <- expression_matrix_from_df(df, samples)
      color_col <- value_or(input$pca_color_col, "__none__")
      if (identical(color_col, "__none__")) color_col <- ""
      draw_pca_plot(
        mat,
        metadata = design_df,
        color_col = color_col,
        include_values = value_or(input$pca_include_values, character(0)),
        label_samples = isTRUE(input$pca_label_samples)
      )
    }, width = function() {
      as.numeric(value_or(input$pca_plot_width, 900))
    }, height = function() {
      as.numeric(value_or(input$pca_plot_height, 700))
    })

    observe({
      req(identical(input$main_tabs, "Differential Expression"))
      df <- deg_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "deg_gene_query", choices = character(0), selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session, "deg_gene_query", choices = df$gene_label, selected = character(0), server = TRUE)
      }
    })

    observe({
      req(identical(input$main_tabs, "Differential Expression"))
      df <- deg_df()
      if (is.null(df)) {
        updateTextAreaInput(session, "deg_up_genes", value = "")
        updateTextAreaInput(session, "deg_down_genes", value = "")
      } else {
        gene_lists <- make_enrichr_gene_lists(
          df,
          p_col = input$deg_p_col,
          p_cutoff = value_or(input$deg_p_cutoff, 0.05),
          lfc_cutoff = value_or(input$deg_lfc_cutoff, 0)
        )
        updateTextAreaInput(session, "deg_up_genes", value = paste(gene_lists$up, collapse = "\n"))
        updateTextAreaInput(session, "deg_down_genes", value = paste(gene_lists$down, collapse = "\n"))
      }
    })

    observe({
      req(identical(input$main_tabs, "Plots"))
      req(identical(input$plot_subtab, "Heatmap"))
      df <- heatmap_source_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "heatmap_genes", choices = character(0), selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session, "heatmap_genes", choices = df$gene_label, selected = isolate(input$heatmap_genes), server = TRUE)
      }
    })

    observe({
      req(identical(input$main_tabs, "Plots"))
      req(identical(input$plot_subtab, "Volcano"))
      df <- plot_deg_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "volcano_label_genes", choices = character(0), selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session, "volcano_label_genes", choices = df$gene_label, selected = isolate(input$volcano_label_genes), server = TRUE)
      }
    })

    volcano_plot_state <- reactive({
      df <- plot_deg_df()
      req(!is.null(df))
      p_col <- value_or(input$volcano_p_col, "padj")
      validate(need(p_col %in% colnames(df), "Selected p-value column is not present in the DEG table."))
      plot_df <- prepare_volcano_df(
        df,
        p_col = p_col,
        p_cutoff = as.numeric(value_or(input$volcano_p_cutoff, 0.05)),
        lfc_cutoff = as.numeric(value_or(input$volcano_lfc_cutoff, 1))
      )
      validate(need(!is.null(plot_df) && nrow(plot_df) > 0, "No usable log2FoldChange and p-value rows were found."))
      labels <- volcano_label_genes(
        plot_df,
        mode = value_or(input$volcano_label_mode, "top"),
        manual_genes = value_or(input$volcano_label_genes, character(0)),
        max_labels = value_or(input$volcano_max_labels, 20)
      )
      list(plot_df = plot_df, labels = labels, p_col = p_col)
    })

    output$deg_volcano_plot <- renderPlot({
      vp <- volcano_plot_state()
      draw_volcano_plot(
        vp$plot_df,
        labels = vp$labels,
        title = sprintf("%s vs %s", value_or(input$plot_treatment, "treatment"), value_or(input$plot_control, "control")),
        subtitle = sprintf("%s, %s", vp$p_col, if (identical(value_or(input$plot_deseq_label_mode, "gene_id"), "gene_name")) "gene names" else "gene IDs"),
        p_cutoff = as.numeric(value_or(input$volcano_p_cutoff, 0.05)),
        lfc_cutoff = as.numeric(value_or(input$volcano_lfc_cutoff, 1)),
        style = value_or(input$volcano_style, "publication"),
        palette_style = value_or(input$volcano_palette, "publication"),
        point_size = value_or(input$volcano_point_size, 2.5),
        point_alpha = value_or(input$volcano_point_alpha, 0.72),
        label_size = value_or(input$volcano_label_size, 2.5),
        show_thresholds = isTRUE(input$volcano_show_thresholds),
        show_grid = isTRUE(input$volcano_show_grid),
        font_family = heatmap_font_family(input$volcano_font_family)
      )
    }, width = function() {
      as.numeric(value_or(input$volcano_plot_width, 900))
    }, height = function() {
      as.numeric(value_or(input$volcano_plot_height, 700))
    })

    observeEvent(input$save_volcano_btn, {
      vp <- volcano_plot_state()
      label_suffix <- if (identical(value_or(input$plot_deseq_label_mode, "gene_id"), "gene_name")) "gene_name" else "gene_id"
      default_name <- sprintf(
        "volcano_%s_%s_vs_%s_%s.png",
        value_or(input$plot_compare_col, "comparison"),
        value_or(input$plot_treatment, "treatment"),
        value_or(input$plot_control, "control"),
        label_suffix
      )
      filename <- sanitize_filename(input$volcano_filename)
      if (!nzchar(filename)) filename <- default_name
      if (!grepl("\\.png$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".png")
      out_path <- file.path(deseq2_result_dir(value_or(input$plot_deseq_label_mode, "gene_id")), filename)
      display_width <- as.numeric(value_or(input$volcano_plot_width, 900))
      display_height <- as.numeric(value_or(input$volcano_plot_height, 700))
      save_mode <- value_or(input$volcano_save_mode, "double")
      save_scale <- if (identical(save_mode, "quintuple")) 5 else if (identical(save_mode, "triple")) 3 else if (identical(save_mode, "double")) 2 else 1
      save_width <- display_width * save_scale
      save_height <- display_height * save_scale
      save_res <- 96 * save_scale
      if (identical(save_mode, "custom")) {
        save_width <- as.numeric(value_or(input$volcano_save_width, 1800))
        save_height <- as.numeric(value_or(input$volcano_save_height, 1400))
        save_res <- as.numeric(value_or(input$volcano_save_res, 192))
      }
      png(out_path, width = save_width, height = save_height, units = "px", res = save_res)
      draw_volcano_plot(
        vp$plot_df,
        labels = vp$labels,
        title = sprintf("%s vs %s", value_or(input$plot_treatment, "treatment"), value_or(input$plot_control, "control")),
        subtitle = sprintf("%s, %s", vp$p_col, if (identical(value_or(input$plot_deseq_label_mode, "gene_id"), "gene_name")) "gene names" else "gene IDs"),
        p_cutoff = as.numeric(value_or(input$volcano_p_cutoff, 0.05)),
        lfc_cutoff = as.numeric(value_or(input$volcano_lfc_cutoff, 1)),
        style = value_or(input$volcano_style, "publication"),
        palette_style = value_or(input$volcano_palette, "publication"),
        point_size = value_or(input$volcano_point_size, 2.5),
        point_alpha = value_or(input$volcano_point_alpha, 0.72),
        label_size = value_or(input$volcano_label_size, 2.5),
        show_thresholds = isTRUE(input$volcano_show_thresholds),
        show_grid = isTRUE(input$volcano_show_grid),
        font_family = heatmap_font_family(input$volcano_font_family)
      )
      dev.off()
      showNotification(sprintf("Saved volcano plot to %s", out_path), type = "message", duration = 6)
    })

    heatmap_selected_genes <- reactive({
      if (identical(input$heatmap_gene_mode, "manual")) {
        return(value_or(input$heatmap_genes, character(0)))
      }
      df <- plot_deg_df()
      req(!is.null(df))
      rank_col <- value_or(input$heatmap_rank_p_col, "padj")
      req(rank_col %in% colnames(df), "log2FoldChange" %in% colnames(df))
      pvals <- suppressWarnings(as.numeric(df[[rank_col]]))
      lfc <- suppressWarnings(as.numeric(df$log2FoldChange))
      genes <- as.character(df$gene_label)
      keep <- !is.na(genes) & nzchar(genes) & !is.na(pvals) & !is.na(lfc)
      work <- data.frame(gene = genes[keep], pval = pvals[keep], lfc = lfc[keep], stringsAsFactors = FALSE)
      work <- work[order(work$pval, -abs(work$lfc)), , drop = FALSE]
      total_n <- max(1, as.integer(value_or(input$heatmap_total_genes, 20)))
      if (isTRUE(input$heatmap_equal_split)) {
        up_n <- floor(total_n / 2)
        down_n <- floor(total_n / 2)
        if (total_n %% 2 == 1) up_n <- up_n + 1
      } else {
        up_n <- total_n
        down_n <- total_n
      }
      up_genes <- head(work$gene[work$lfc > 0], up_n)
      down_genes <- head(work$gene[work$lfc < 0], down_n)
      if (!isTRUE(input$heatmap_equal_split)) {
        combined <- unique(c(up_genes, down_genes))
        return(head(combined, total_n))
      }
      unique(c(up_genes, down_genes))
    })

    if (DT_AVAILABLE) {
      output$deg_table <- DT::renderDT({
        df <- deg_df()
        req(!is.null(df))
        df <- sort_deg_table(df, sort_mode = value_or(input$deg_sort_mode, "up"), p_col = input$deg_p_col)
        q <- input$deg_gene_query
        if (is.null(q) || !nzchar(trimws(q))) {
          show_df <- df
        } else {
          hits <- df[tolower(df$gene_label) == tolower(trimws(q)), , drop = FALSE]
          show_df <- if (nrow(hits) == 0) df else hits
        }
        simple_dt(show_df, page_length = 50)
      }, server = FALSE)
    } else {
      output$deg_table <- renderTable({
        df <- deg_df()
        req(!is.null(df))
        df <- sort_deg_table(df, sort_mode = value_or(input$deg_sort_mode, "up"), p_col = input$deg_p_col)
        q <- input$deg_gene_query
        if (is.null(q) || !nzchar(trimws(q))) {
          show_df <- format_deg_table(df)
        } else {
          hits <- df[tolower(df$gene_label) == tolower(trimws(q)), , drop = FALSE]
          show_df <- if (nrow(hits) == 0) format_deg_table(df) else format_deg_table(hits)
        }
        rownames(show_df) <- make.unique(as.character(show_df$gene_label))
        show_df$gene_label <- NULL
        show_df
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
    }

    output$deg_all_pca_ui <- renderUI({
      if (!deseq2_available) return(tags$p("PCA has not been generated yet."))
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      pca_path <- pca_all_plot_path(input$deg_compare_col)
      if (!file.exists(pca_path)) {
        return(tags$p("All-sample PCA has not been generated yet. Re-run DESeq2 once to create it."))
      }
      tags$img(src = file.path("deseq2_results", basename(pca_path)), style = "max-width: 100%; border: 1px solid #ddd;")
    })

    output$deg_pca_ui <- renderUI({
      if (!deseq2_available) return(tags$p("PCA has not been generated yet."))
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      req(input$deg_treatment, input$deg_control)
      pca_path <- pca_plot_path(input$deg_compare_col, input$deg_treatment, input$deg_control)
      if (!file.exists(pca_path)) {
        return(tags$p("PCA has not been generated for this comparison."))
      }
      tags$img(src = file.path("deseq2_results", basename(pca_path)), style = "max-width: 100%; border: 1px solid #ddd;")
    })

    heatmap_plot_state <- reactive({
      df <- heatmap_source_df()
      req(!is.null(df))
      genes <- heatmap_selected_genes()
      validate(need(length(genes) > 0, "Select one or more genes to plot."))
      hits <- df[df$gene_label %in% genes, , drop = FALSE]
      validate(need(nrow(hits) > 0, "None of the selected genes were found in the normalized counts table."))
      rownames(hits) <- hits$gene_label
      mat <- as.matrix(hits[, setdiff(colnames(hits), "gene_label"), drop = FALSE])
      storage.mode(mat) <- "numeric"
      ann_col <- NULL
      ann_choices <- value_or(input$heatmap_annotation_cols, character(0))
      if (length(ann_choices) > 0) {
        ann_col <- design_df[match(colnames(mat), design_df$sample), ann_choices, drop = FALSE]
        rownames(ann_col) <- colnames(mat)
      }
      if (!isTRUE(input$heatmap_cluster_cols)) {
        order_col <- value_or(input$heatmap_order_col, "")
        if (nzchar(order_col) && order_col %in% colnames(design_df)) {
          order_vals <- design_df[match(colnames(mat), design_df$sample), order_col, drop = TRUE]
          ord <- order(order_vals, colnames(mat))
          mat <- mat[, ord, drop = FALSE]
          if (!is.null(ann_col)) ann_col <- ann_col[ord, , drop = FALSE]
        }
      }
      if (isTRUE(input$heatmap_scale_rows)) {
        mat <- t(scale(t(mat)))
        mat[is.na(mat)] <- 0
      }
      list(mat = mat, ann_col = ann_col)
    })

    output$deg_heatmap_plot <- renderPlot({
      hp <- heatmap_plot_state()
      heatmap_theme <- value_or(input$heatmap_theme, "publication")
      hm_colors <- heatmap_palette_colors(input$heatmap_palette, heatmap_theme)
      pheatmap::pheatmap(
        hp$mat,
        annotation_col = hp$ann_col,
        cluster_rows = isTRUE(input$heatmap_cluster_rows),
        cluster_cols = isTRUE(input$heatmap_cluster_cols),
        scale = "none",
        fontsize = as.numeric(value_or(input$heatmap_col_font_size, 10)),
        fontsize_row = as.numeric(value_or(input$heatmap_row_font_size, 9)),
        fontsize_col = as.numeric(value_or(input$heatmap_col_font_size, 10)),
        fontfamily = heatmap_font_family(input$heatmap_font_family),
        angle_col = value_or(input$heatmap_col_angle, "45"),
        show_rownames = isTRUE(input$heatmap_show_row_names),
        show_colnames = isTRUE(input$heatmap_show_col_names),
        annotation_names_col = isTRUE(input$heatmap_show_annotation_names),
        clustering_distance_rows = value_or(input$heatmap_distance, "euclidean"),
        clustering_distance_cols = value_or(input$heatmap_distance, "euclidean"),
        clustering_method = value_or(input$heatmap_cluster_method, "complete"),
        border_color = heatmap_border_color(input$heatmap_border_style, heatmap_theme),
        color = hm_colors,
        breaks = heatmap_color_breaks(hp$mat, hm_colors)
      )
    }, width = function() {
      as.numeric(value_or(input$heatmap_plot_width, 900))
    }, height = function() {
      as.numeric(value_or(input$heatmap_plot_height, 700))
    })

    observeEvent(input$save_heatmap_btn, {
      hp <- heatmap_plot_state()
      label_suffix <- if (identical(value_or(input$plot_deseq_label_mode, "gene_id"), "gene_name")) "gene_name" else "gene_id"
      default_name <- sprintf(
        "heatmap_%s_%s_vs_%s_%s.png",
        value_or(input$plot_compare_col, "comparison"),
        value_or(input$plot_treatment, "treatment"),
        value_or(input$plot_control, "control"),
        label_suffix
      )
      filename <- sanitize_filename(input$heatmap_filename)
      if (!nzchar(filename)) filename <- default_name
      if (!grepl("\\.png$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".png")
      out_path <- file.path(deseq2_result_dir(value_or(input$plot_deseq_label_mode, "gene_id")), filename)
      display_width <- as.numeric(value_or(input$heatmap_plot_width, 900))
      display_height <- as.numeric(value_or(input$heatmap_plot_height, 700))
      save_mode <- value_or(input$heatmap_save_mode, "double")
      save_scale <- if (identical(save_mode, "quintuple")) 5 else if (identical(save_mode, "triple")) 3 else if (identical(save_mode, "double")) 2 else 1
      save_width <- display_width * save_scale
      save_height <- display_height * save_scale
      save_res <- 96 * save_scale
      if (identical(save_mode, "custom")) {
        save_width <- as.numeric(value_or(input$heatmap_save_width, 1800))
        save_height <- as.numeric(value_or(input$heatmap_save_height, 1400))
        save_res <- as.numeric(value_or(input$heatmap_save_res, 96))
      }
      png(
        out_path,
        width = save_width,
        height = save_height,
        units = "px",
        res = save_res
      )
      heatmap_theme <- value_or(input$heatmap_theme, "publication")
      hm_colors <- heatmap_palette_colors(input$heatmap_palette, heatmap_theme)
      pheatmap::pheatmap(
        hp$mat,
        annotation_col = hp$ann_col,
        cluster_rows = isTRUE(input$heatmap_cluster_rows),
        cluster_cols = isTRUE(input$heatmap_cluster_cols),
        scale = "none",
        fontsize = as.numeric(value_or(input$heatmap_col_font_size, 10)),
        fontsize_row = as.numeric(value_or(input$heatmap_row_font_size, 9)),
        fontsize_col = as.numeric(value_or(input$heatmap_col_font_size, 10)),
        fontfamily = heatmap_font_family(input$heatmap_font_family),
        angle_col = value_or(input$heatmap_col_angle, "45"),
        show_rownames = isTRUE(input$heatmap_show_row_names),
        show_colnames = isTRUE(input$heatmap_show_col_names),
        annotation_names_col = isTRUE(input$heatmap_show_annotation_names),
        clustering_distance_rows = value_or(input$heatmap_distance, "euclidean"),
        clustering_distance_cols = value_or(input$heatmap_distance, "euclidean"),
        clustering_method = value_or(input$heatmap_cluster_method, "complete"),
        border_color = heatmap_border_color(input$heatmap_border_style, heatmap_theme),
        color = hm_colors,
        breaks = heatmap_color_breaks(hp$mat, hm_colors)
      )
      dev.off()
      showNotification(sprintf("Saved heatmap to %s", out_path), type = "message", duration = 6)
    })
  } else {
    empty_compare_input <- function(id, label) {
      renderUI(selectInput(id, label, choices = character(0), selected = character(0), selectize = FALSE))
    }
    output$deseq_treatment_ui <- empty_compare_input("deseq_treatment", "Treatment")
    output$deseq_control_ui <- empty_compare_input("deseq_control", "Control")
    output$deseq_status_ui <- renderUI({
      status_box("DESeq2 normalized counts have not been generated yet.", "warning")
    })
    output$deseq_counts_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })
    output$download_deseq_counts_table <- downloadHandler(
      filename = function() download_filename("counts_deseq2_normalized_unavailable"),
      content = function(file) write_csv_download(safe_download_df(NULL), file)
    )

    output$deg_treatment_ui <- empty_compare_input("deg_treatment", "Treatment")
    output$deg_control_ui <- empty_compare_input("deg_control", "Control")
    output$plot_treatment_ui <- empty_compare_input("plot_treatment", "Treatment")
    output$plot_control_ui <- empty_compare_input("plot_control", "Control")
    output$deg_status_ui <- renderUI({
      status_box("Differential expression results have not been generated yet.", "warning")
    })
    output$plot_status_ui <- renderUI({
      status_box("Plots are available after differential expression outputs are generated. PCA can also use the count matrix when available.", "warning")
    })
    output$deg_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })
    output$pca_include_values_ui <- renderUI({
      col <- value_or(input$pca_color_col, "__none__")
      if (identical(col, "__none__") || !col %in% colnames(design_df)) return(NULL)
      vals <- sort(unique(as.character(design_df[[col]])))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      selectizeInput("pca_include_values", "Include values", choices = vals, selected = vals, multiple = TRUE, options = list(dropdownParent = "body"))
    })
    output$pca_plot <- renderPlot({
      df <- count_matrix_nonzero_df
      if (is.null(df) || !nrow(df)) {
        plot.new()
        text(0.5, 0.5, "PCA is not available until a count matrix or normalized-counts file exists.")
      } else {
        colnames(df)[1] <- "gene_label"
        mat <- expression_matrix_from_df(df, samples)
        color_col <- value_or(input$pca_color_col, "__none__")
        if (identical(color_col, "__none__")) color_col <- ""
        draw_pca_plot(mat, design_df, color_col, value_or(input$pca_include_values, character(0)), isTRUE(input$pca_label_samples))
      }
    }, width = function() {
      as.numeric(value_or(input$pca_plot_width, 900))
    }, height = function() {
      as.numeric(value_or(input$pca_plot_height, 700))
    })
    output$deg_all_pca_ui <- renderUI({
      tags$p("PCA has not been generated yet.")
    })
    output$deg_pca_ui <- renderUI({
      tags$p("PCA has not been generated yet.")
    })
    output$deg_volcano_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "Volcano plot not available until differential expression outputs are generated.")
    })
    output$deg_heatmap_plot <- renderPlot({
      plot.new()
      text(0.5, 0.5, "Heatmap not available until differential expression outputs are generated.")
    })
    observe({
      updateTextAreaInput(session, "deg_up_genes", value = "")
      updateTextAreaInput(session, "deg_down_genes", value = "")
      updateSelectizeInput(session, "deseq_gene_query", choices = character(0), selected = character(0), server = TRUE)
      updateSelectizeInput(session, "deg_gene_query", choices = character(0), selected = character(0), server = TRUE)
      updateSelectizeInput(session, "volcano_label_genes", choices = character(0), selected = character(0), server = TRUE)
      updateSelectizeInput(session, "heatmap_genes", choices = character(0), selected = character(0), server = TRUE)
    })
  }

  {
    get_gsea_values <- reactive({
      if (is.null(input$gsea_compare_col) || identical(input$gsea_compare_col, "NA") || !nzchar(input$gsea_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$gsea_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$gsea_treatment_ui <- renderUI({
      vals <- get_gsea_values()
      selectInput("gsea_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals, selectize = FALSE)
    })

    output$gsea_control_ui <- renderUI({
      vals <- get_gsea_values()
      selectInput("gsea_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals, selectize = FALSE)
    })

    gsea_comp_dir <- reactive({
      req(input$gsea_compare_col)
      if (identical(input$gsea_compare_col, "NA")) return(NULL)
      req(input$gsea_treatment, input$gsea_control)
      gseapy_comparison_dir(input$gsea_treatment, input$gsea_control)
    })

    output$gsea_collection_ui <- renderUI({
      req(identical(input$main_tabs, "GSEA"))
      comp_dir <- gsea_comp_dir()
      collections <- list_gseapy_collections(comp_dir)
      selectInput("gsea_collection", "Pathway collection", choices = collections, selected = if (length(collections)) collections[1] else character(0), selectize = FALSE)
    })

    gsea_report_df <- reactive({
      comp_dir <- gsea_comp_dir()
      req(!is.null(comp_dir), input$gsea_collection)
      read_gseapy_report(comp_dir, input$gsea_collection)
    })

    output$gsea_pathway_ui <- renderUI({
      req(identical(input$main_tabs, "GSEA"))
      df <- gsea_report_df()
      pathways <- if (is.null(df)) character(0) else df$Term
      selectizeInput("gsea_pathway", "Pathway", choices = pathways, selected = if (length(pathways)) pathways[1] else character(0), multiple = FALSE, options = list(dropdownParent = "body"))
    })

    output$gsea_status_ui <- renderUI({
      if (!gseapy_results_available()) {
        return(status_box("GSEA results have not been generated yet.", "warning"))
      }
      if (!nrow(design_df) || !length(comparison_columns)) {
        return(status_box("The design matrix is not available yet, so comparisons cannot be selected.", "warning"))
      }
      req(input$gsea_compare_col)
      if (identical(input$gsea_compare_col, "NA")) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fff7e6; border: 1px solid #f0c36d; border-radius: 6px;",
          "Select a comparison column to view GSEA results."
        ))
      }
      req(input$gsea_treatment, input$gsea_control)
      comp_dir <- gsea_comp_dir()
      if (is.null(comp_dir) || !dir.exists(comp_dir)) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;",
          sprintf("This comparison has not been tested for GSEA: %s vs %s.", input$gsea_treatment, input$gsea_control)
        ))
      }
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf("Showing GSEA for %s vs %s using %s.", input$gsea_treatment, input$gsea_control, input$gsea_compare_col)
      )
    })

    if (DT_AVAILABLE) {
      output$gsea_pathway_table <- DT::renderDT({
        df <- gsea_report_df()
        req(!is.null(df), input$gsea_pathway)
        hit <- df[df$Term == input$gsea_pathway, , drop = FALSE]
        if (!nrow(hit)) return(NULL)
        simple_dt(hit, page_length = 5, scroll_y = "180px", dom = "t")
      }, server = FALSE)

      output$gsea_all_pathways_table <- DT::renderDT({
        df <- gsea_report_df()
        req(!is.null(df))
        simple_dt(df, page_length = 50)
      }, server = FALSE)
    } else {
      output$gsea_pathway_table <- renderTable({
        df <- gsea_report_df()
        req(!is.null(df), input$gsea_pathway)
        hit <- df[df$Term == input$gsea_pathway, , drop = FALSE]
        if (!nrow(hit)) return(NULL)
        format_gsea_table(hit)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

      output$gsea_all_pathways_table <- renderTable({
        df <- gsea_report_df()
        req(!is.null(df))
        format_gsea_table(df)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)
    }

    output$gsea_summary_plots_ui <- renderUI({
      req(input$gsea_treatment, input$gsea_control, input$gsea_collection)
      comp_rel <- gseapy_comparison_rel(input$gsea_treatment, input$gsea_control)
      dotplot_rel <- file.path("gseapy_results", comp_rel, sprintf("DotPlot_Top10.%s.png", input$gsea_collection))
      enrich_rel <- file.path("gseapy_results", comp_rel, sprintf("EnrichmentPlot_Top10.%s.png", input$gsea_collection))
      tags$div(
        if (file.exists(file.path(gseapy_dir, comp_rel, sprintf("DotPlot_Top10.%s.png", input$gsea_collection))))
          tags$img(src = dotplot_rel, style = "max-width: 100%; border: 1px solid #ddd; margin-bottom: 12px;"),
        if (file.exists(file.path(gseapy_dir, comp_rel, sprintf("EnrichmentPlot_Top10.%s.png", input$gsea_collection))))
          tags$img(src = enrich_rel, style = "max-width: 100%; border: 1px solid #ddd;")
      )
    })

    output$gsea_pathway_plots_ui <- renderUI({
      req(input$gsea_treatment, input$gsea_control, input$gsea_pathway)
      pdf_rel <- pathway_pdf_relpath(input$gsea_treatment, input$gsea_control, input$gsea_pathway, heatmap = FALSE)
      heatmap_rel <- pathway_pdf_relpath(input$gsea_treatment, input$gsea_control, input$gsea_pathway, heatmap = TRUE)
      if (is.null(pdf_rel) && is.null(heatmap_rel)) {
        return(tags$p("No pathway-specific plot files found for this pathway."))
      }
      plot_block <- function(title, rel) {
        if (is.null(rel)) return(NULL)
        rendered_rel <- pdf_render_relpath(rel)
        if (!is.null(rendered_rel)) {
          return(tags$div(
            tags$h5(title),
            tags$img(src = rendered_rel, style = "width: 100%; max-width: 100%; border: 1px solid #ddd; margin-bottom: 12px;"),
            tags$p(tags$a(href = rel, target = "_blank", "Open original PDF"))
          ))
        }
        pdf_canvas_block(title, rel)
      }
      tags$div(
        plot_block("Enrichment Plot", pdf_rel),
        plot_block("Pathway Heatmap", heatmap_rel)
      )
    })
  }

  if (kallisto_available) {
    kallisto_matrix_cache <- reactiveVal(NULL)
    kallisto_matrix_cache_path <- reactiveVal(NULL)
    kallisto_conversion_status <- reactiveVal(NULL)
    kallisto_table_loaded <- reactiveVal(FALSE)

    observeEvent(input$kallisto_load_table, {
      kallisto_table_loaded(TRUE)
    }, ignoreInit = TRUE)

    observeEvent(list(input$kallisto_view_mode, input$kallisto_sample, input$kallisto_sample_filter_col, input$kallisto_filter_col, input$kallisto_matrix_metric, input$kallisto_label_mode), {
      kallisto_table_loaded(FALSE)
      kallisto_conversion_status(NULL)
    }, ignoreInit = TRUE)

    get_kallisto_matrix <- function() {
      req(input$kallisto_matrix_metric)
      label_mode <- value_or(input$kallisto_label_mode, "target_id")
      if (identical(label_mode, "gene_name")) {
        base_path <- kallisto_saved_matrix_path(input$kallisto_matrix_metric, "target_id")
        converted_path <- gene_name_matrix_path(base_path)
        if (!file.exists(converted_path)) {
          res <- save_gene_name_matrix(base_path, "kallisto")
          kallisto_conversion_status(res$message)
        }
      }
      matrix_path <- kallisto_saved_matrix_path(input$kallisto_matrix_metric, label_mode)
      cached <- kallisto_matrix_cache()
      if (!is.null(cached) && identical(kallisto_matrix_cache_path(), matrix_path)) return(cached)
      if (file.exists(matrix_path)) {
        cached <- read_saved_count_matrix(matrix_path)
        kallisto_matrix_cache(cached)
        kallisto_matrix_cache_path(matrix_path)
        return(cached)
      }
      NULL
    }

    output$kallisto_convert_ui <- renderUI({
      req(input$kallisto_matrix_metric)
      base_path <- kallisto_saved_matrix_path(input$kallisto_matrix_metric, "target_id")
      converted_path <- gene_name_matrix_path(base_path)
      if (!file.exists(base_path)) {
        return(tags$span(class = "tiny-note", "Saved Kallisto matrix not found yet. Run the Kallisto step first."))
      }
      if (identical(value_or(input$kallisto_label_mode, "target_id"), "gene_name")) {
        if (!file.exists(converted_path)) save_gene_name_matrix(base_path, "kallisto")
        return(tags$span(class = "tiny-note", if (file.exists(converted_path)) "Using saved Kallisto gene-name matrix." else "Kallisto gene-name matrix could not be saved."))
      }
      NULL
    })

    output$kallisto_status_ui <- renderUI({
      if (!kallisto_available) {
        return(status_box("Kallisto transcript abundance has not been generated yet.", "warning"))
      }
      if (!isTRUE(kallisto_table_loaded())) {
        return(status_box("Choose Kallisto settings, then click 'Load Kallisto table' to render the table.", "info"))
      }
      if (identical(input$kallisto_view_mode, "matrix")) {
        matrix_path <- kallisto_saved_matrix_path(value_or(input$kallisto_matrix_metric, "tpm"), value_or(input$kallisto_label_mode, "target_id"))
        if (file.exists(matrix_path)) {
          msg <- kallisto_conversion_status()
          if (is.null(msg) || !nzchar(msg)) msg <- sprintf("Using saved Kallisto matrix: %s", matrix_path)
          return(status_box(msg, "info"))
        }
        return(status_box("Saved Kallisto matrix has not been generated yet. Run the Kallisto step to build it in the counts folder.", "warning"))
      }
      NULL
    })

    observe({
      req(identical(input$main_tabs, "Counts"))
      req(identical(input$counts_subtab, "Kallisto"))
      req(identical(input$kallisto_view_mode, "sample"))
      req(isTRUE(kallisto_table_loaded()))
      req(input$kallisto_sample_filter_col)
      if (identical(input$kallisto_sample_filter_col, "all")) {
        updateSelectizeInput(
          session,
          "kallisto_sample_filter_value",
          choices = character(0),
          selected = character(0),
          server = TRUE
        )
      } else {
        df <- read_kallisto_abundance(input$kallisto_sample)
        req(!is.null(df))
        vals <- sort(unique(df[[input$kallisto_sample_filter_col]]))
        vals <- vals[!is.na(vals) & nzchar(vals)]
        updateSelectizeInput(
          session,
          "kallisto_sample_filter_value",
          choices = vals,
          selected = character(0),
          server = TRUE
        )
      }
    })

    observe({
      req(identical(input$main_tabs, "Counts"))
      req(identical(input$counts_subtab, "Kallisto"))
      req(identical(input$kallisto_view_mode, "matrix"))
      req(isTRUE(kallisto_table_loaded()))
      req(input$kallisto_filter_col)
      matrix_df <- get_kallisto_matrix()
      if (identical(input$kallisto_filter_col, "all") || is.null(matrix_df) || !(input$kallisto_filter_col %in% colnames(matrix_df))) {
        vals <- character(0)
      } else {
        vals <- sort(unique(as.character(matrix_df[[input$kallisto_filter_col]])))
        vals <- vals[!is.na(vals) & nzchar(vals)]
      }
      updateSelectizeInput(
        session,
        "kallisto_filter_value",
        choices = vals,
        selected = character(0),
        server = TRUE
      )
    })

    kallisto_table_df <- reactive({
      req(isTRUE(kallisto_table_loaded()))
      req(input$kallisto_view_mode)
      if (identical(input$kallisto_view_mode, "sample")) {
        req(input$kallisto_sample)
        df <- read_kallisto_abundance(input$kallisto_sample)
        req(!is.null(df))
        show_df <- df[, c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype", "tpm"), drop = FALSE]
        colnames(show_df)[colnames(show_df) == "tpm"] <- input$kallisto_sample
        if (!identical(input$kallisto_sample_filter_col, "all") &&
            !is.null(input$kallisto_sample_filter_value) &&
            nzchar(trimws(input$kallisto_sample_filter_value))) {
          show_df <- show_df[
            tolower(show_df[[input$kallisto_sample_filter_col]]) == tolower(trimws(input$kallisto_sample_filter_value)),
            ,
            drop = FALSE
          ]
        }
        if (nrow(show_df) == 0) {
          show_df <- df[, c("gene_symbol", "transcript_name", "transcript_id", "gene_id", "biotype", "tpm"), drop = FALSE]
          colnames(show_df)[colnames(show_df) == "tpm"] <- input$kallisto_sample
        }
      } else {
        kallisto_matrix_df <- get_kallisto_matrix()
        req(!is.null(kallisto_matrix_df))
        if (identical(input$kallisto_filter_col, "all") || is.null(input$kallisto_filter_value) || !nzchar(trimws(input$kallisto_filter_value)) || !(input$kallisto_filter_col %in% colnames(kallisto_matrix_df))) {
          show_df <- kallisto_matrix_df
        } else {
          show_df <- kallisto_matrix_df[
            tolower(kallisto_matrix_df[[input$kallisto_filter_col]]) == tolower(trimws(input$kallisto_filter_value)),
            ,
            drop = FALSE
          ]
          if (nrow(show_df) == 0) {
            show_df <- kallisto_matrix_df
          }
        }
      }
      sort_df_by_col(show_df, colnames(show_df)[1], "asc")
    })

    observe({
      df <- kallisto_table_df()
      req(!is.null(df))
      updateSelectInput(
        session,
        "kallisto_sort_col",
        choices = colnames(df),
        selected = colnames(df)[1]
      )
    })

    if (DT_AVAILABLE) {
      output$kallisto_table <- DT::renderDT({
        if (!isTRUE(kallisto_table_loaded())) {
          return(simple_dt(data.frame(Message = "Click 'Load Kallisto table' to render this table.", stringsAsFactors = FALSE), page_length = 50, scroll_y = "140px"))
        }
        df <- kallisto_table_df()
        simple_dt(df, page_length = 50)
      }, server = FALSE)
    } else {
      output$kallisto_table <- renderTable({
        if (!isTRUE(kallisto_table_loaded())) return(data.frame(Message = "Click 'Load Kallisto table' to render this table.", stringsAsFactors = FALSE))
        format_numeric_commas(kallisto_table_df())
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)
    }
    output$download_kallisto_table <- downloadHandler(
      filename = function() download_filename(paste0("counts_kallisto_", value_or(input$kallisto_view_mode, "matrix"))),
      content = function(file) write_csv_download(safe_download_df(if (isTRUE(kallisto_table_loaded())) kallisto_table_df() else NULL), file)
    )
  } else {
    output$kallisto_convert_ui <- renderUI(NULL)
    output$kallisto_status_ui <- renderUI({
      status_box("Kallisto transcript abundance has not been generated yet.", "warning")
    })
    output$kallisto_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable(NULL)
    output$download_kallisto_table <- downloadHandler(
      filename = function() download_filename("counts_kallisto_unavailable"),
      content = function(file) write_csv_download(safe_download_df(NULL), file)
    )
    observe({
      updateSelectizeInput(session, "kallisto_sample_filter_value", choices = character(0), selected = character(0), server = TRUE)
      updateSelectizeInput(session, "kallisto_filter_value", choices = character(0), selected = character(0), server = TRUE)
    })
  }
}

shinyApp(ui, server)
