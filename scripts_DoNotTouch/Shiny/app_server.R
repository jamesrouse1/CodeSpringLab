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
count_matrix_path <- file.path(data_dir, "counts", "count_matrix.txt")
featurecounts_summary_path <- file.path(data_dir, "counts", "featurecounts_summary.txt")
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
  tryCatch(
    read.delim(path, ...),
    error = function(e) NULL
  )
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
  ids <- sub("_(R1|R2)_001_(fastqc|screen)\\.html$", "", files)
  sort(unique(ids[ids != files]))
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
kallisto_available <- length(kallisto_sample_names) > 0
rsem_available <- dir.exists(rsem_dir) && length(list.files(rsem_dir, pattern = "\\.(genes|isoforms)\\.results$", recursive = TRUE)) > 0
deseq2_available <- dir.exists(deseq2_dir) && length(list.files(deseq2_dir, pattern = "^normalized_counts_.*\\.txt$", recursive = FALSE)) > 0
gseapy_available <- dir.exists(gseapy_dir) && length(list.files(gseapy_dir, pattern = "report\\.gseapy\\..*\\.csv$", recursive = TRUE)) > 0

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

if (dir.exists(fastqc_dir)) addResourcePath("fastqc_results", fastqc_dir)
if (dir.exists(fastqc_cutadapt_dir)) addResourcePath("fastqc_cutadapt_results", fastqc_cutadapt_dir)
if (dir.exists(deseq2_dir)) addResourcePath("deseq2_results", deseq2_dir)
if (dir.exists(gseapy_dir)) addResourcePath("gseapy_results", gseapy_dir)
if (!is.null(logo_csl_path)) addResourcePath("logo_csl_asset", dirname(logo_csl_path))
if (!is.null(logo_path)) addResourcePath("logo_asset", dirname(logo_path))

normalize_sample_token <- function(x) {
  x <- tolower(trimws(as.character(value_or(x, ""))))
  gsub("[^a-z0-9]+", "", x)
}

find_qc_html <- function(base_dir, sample_name, read_name, suffix) {
  if (!dir.exists(base_dir)) return(file.path(base_dir, paste0(sample_name, "_", read_name, "_001_", suffix, ".html")))
  sample_stem <- sample_fastq_stems[[sample_name]]
  candidates <- unique(c(sample_name, sample_stem, gsub("_", "", sample_name), gsub("_", "", sample_stem)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  expected <- file.path(base_dir, paste0(candidates, "_", read_name, "_001_", suffix, ".html"))
  hit <- expected[file.exists(expected)]
  if (length(hit)) return(hit[1])

  files <- list.files(base_dir, pattern = paste0("_", read_name, "_001_", suffix, "\\.html$"), full.names = TRUE)
  if (!length(files)) return(expected[1])
  file_keys <- normalize_sample_token(sub(paste0("_", read_name, "_001_", suffix, "\\.html$"), "", basename(files)))
  cand_keys <- normalize_sample_token(candidates)
  idx <- match(cand_keys, file_keys, nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx)) return(files[idx[1]])
  expected[1]
}

qc_file_map <- function(sample_name, read_name = c("R1", "R2"), trimmed = FALSE) {
  read_name <- match.arg(read_name)
  base_dir <- if (trimmed) fastqc_cutadapt_dir else fastqc_dir
  list(
    fastqc = find_qc_html(base_dir, sample_name, read_name, "fastqc"),
    screen = find_qc_html(base_dir, sample_name, read_name, "screen")
  )
}

iframe_or_message <- function(path, resource_prefix, height = "calc(100vh - 260px)", missing_message = NULL) {
  if (!file.exists(path)) {
    return(tags$p(value_or(missing_message, sprintf("Missing file: %s", basename(path)))))
  }
  rel <- basename(path)
  height_css <- if (is.numeric(height)) paste0(height, "px") else height
  tags$iframe(
    class = "qc-report-frame",
    src = file.path(resource_prefix, rel),
    style = sprintf("width: 100%%; height: %s; min-height: 760px; border: 1px solid #d7e0ea;", height_css)
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

table_widget <- function(id) {
  if (DT_AVAILABLE) {
    DT::DTOutput(id)
  } else {
    tableOutput(id)
  }
}

simple_dt <- function(df, page_length = 25, scroll_y = "520px", dom = "tip") {
  DT::datatable(
    df,
    rownames = FALSE,
    class = "compact stripe hover cell-border",
    options = list(
      scrollX = TRUE,
      scrollY = scroll_y,
      pageLength = page_length,
      dom = dom
    )
  )
}

format_numeric_commas <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  out <- df
  for (nm in colnames(out)) {
    if (is.numeric(out[[nm]])) {
      out[[nm]] <- ifelse(
        is.na(out[[nm]]),
        NA_character_,
        prettyNum(out[[nm]], big.mark = ",", scientific = FALSE, trim = TRUE)
      )
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
  fname <- if (heatmap) paste0(pathway_name, ".heatmap.pdf") else paste0(pathway_name, ".pdf")
  rel <- file.path(gseapy_comparison_rel(treatment_value, control_value), "gsea", fname)
  full <- file.path(gseapy_dir, rel)
  if (!file.exists(full)) return(NULL)
  file.path("gseapy_results", rel)
}

pdf_render_relpath <- function(resource_rel) {
  if (is.null(resource_rel)) return(NULL)
  rel <- sub("^gseapy_results/", "", resource_rel)
  full_pdf <- file.path(gseapy_dir, rel)
  if (!file.exists(full_pdf)) return(NULL)
  png_rel <- paste0(tools::file_path_sans_ext(rel), ".rendered.png")
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
  if (file.exists(full_png)) file.path("gseapy_results", png_rel) else NULL
}

normalized_counts_path <- function(treatment_value, control_value) {
  file.path(deseq2_dir, sprintf("normalized_counts_%s_vs_%s(ref).txt", treatment_value, control_value))
}

deg_table_path <- function(treatment_value, control_value) {
  file.path(deseq2_dir, sprintf("DEG_%s_vs_%s(ref).txt", treatment_value, control_value))
}

pca_plot_path <- function(compare_col, treatment_value, control_value) {
  file.path(deseq2_dir, sprintf("pca_%s_%s_vs_%s(ref).png", compare_col, treatment_value, control_value))
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
  if (identical(style, "minimal")) return(list(palette = "viridis", border = "none"))
  if (identical(style, "compact")) return(list(palette = "magma", border = "none"))
  if (identical(style, "labeled")) return(list(palette = "blue_red", border = "light"))
  if (identical(style, "contrast")) return(list(palette = "green_purple", border = "light"))
  list(palette = "blue_red", border = "none")
}

heatmap_palette_colors <- function(style, theme = "publication") {
  if (identical(value_or(style, "theme"), "theme")) {
    style <- heatmap_theme_options(theme)$palette
  }
  if (identical(style, "viridis")) return(grDevices::hcl.colors(100, "viridis"))
  if (identical(style, "magma")) return(grDevices::hcl.colors(100, "inferno"))
  if (identical(style, "green_purple")) return(colorRampPalette(c("#1b7837", "white", "#762a83"))(100))
  if (identical(style, "navy_orange")) return(colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100))
  if (identical(style, "gray_red")) return(colorRampPalette(c("#3b3b3b", "white", "#b2182b"))(100))
  colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
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
        selectInput("sample", "Sample", choices = samples, selected = first_or_null(samples)),
        checkboxInput("show_trimmed", "Show cutadapt-trimmed QC", value = default_show_trimmed),
        tags$hr(),
        helpText("This tab renders FastQC and FastQ Screen HTML reports for both reads.")
      ),
      mainPanel(
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
        selectInput("star_sample", "STAR sample", choices = samples, selected = first_or_null(samples)),
        selectInput("star_sample_sort_col", "Sort sample columns by", choices = comparison_columns, selected = if ("treatment" %in% comparison_columns) "treatment" else first_or_null(comparison_columns)),
        tags$hr(),
        helpText("Compact alignment metrics pulled from STAR summary output.")
      ),
      mainPanel(
        uiOutput("star_status_ui"),
        h4("Alignment Summary Across Samples"),
        table_widget("star_summary_table"),
        tags$hr(),
        h4("Selected Sample"),
        tableOutput("star_sample_table")
      )
    )
  ),
  tabPanel(
    "FeatureCounts QC",
    sidebarLayout(
      sidebarPanel(
        selectInput("featurecounts_qc_sample_sort_col", "Sort sample columns by", choices = comparison_columns, selected = if ("treatment" %in% comparison_columns) "treatment" else first_or_null(comparison_columns))
      ),
      mainPanel(
        uiOutput("featurecounts_qc_status_ui"),
        h4("FeatureCounts Summary Across Samples"),
        table_widget("featurecounts_summary_table")
      )
    )
  )
)

counts_subtabs <- list(
  tabPanel(
    "Raw Counts",
    sidebarLayout(
      sidebarPanel(
        selectizeInput("gene_query", "Search gene", choices = NULL, selected = NULL, multiple = FALSE),
        uiOutput("featurecounts_convert_ui"),
        tags$hr(),
        helpText("Search the raw featureCounts matrix and optionally flip between gene_id and gene_name using the local GTF.")
      ),
      mainPanel(
        uiOutput("featurecounts_status_ui"),
        table_widget("gene_search_table")
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
          selectInput("rsem_type", "Show", choices = c("Genes" = "genes", "Isoforms" = "isoforms"), selected = "genes"),
          selectInput("rsem_metric", "Metric", choices = c("TPM", "expected_count", "FPKM"), selected = "TPM"),
          selectizeInput("rsem_query", "Select gene/transcript of interest", choices = NULL, selected = NULL, multiple = FALSE),
          uiOutput("rsem_convert_ui"),
          tags$hr(),
          helpText("Shows an RSEM matrix across all samples. Choose genes or isoforms, then optionally flip displayed gene labels using the local GTF.")
        ),
        mainPanel(
          h4("RSEM Matrix"),
          uiOutput("rsem_status_ui"),
          table_widget("rsem_table")
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
            selected = if ("treatment" %in% comparison_columns) "treatment" else "NA"
          ),
          uiOutput("deseq_treatment_ui"),
          uiOutput("deseq_control_ui"),
          selectizeInput("deseq_gene_query", "Select gene", choices = NULL, selected = NULL, multiple = FALSE),
          uiOutput("deseq_convert_ui"),
          tags$hr(),
          helpText("Shows DESeq2 normalized counts for an available treatment vs control comparison.")
        ),
        mainPanel(
          uiOutput("deseq_status_ui"),
          table_widget("deseq_counts_table")
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
            selectInput("kallisto_sample", "Kallisto sample", choices = kallisto_sample_names, selected = first_or_null(kallisto_sample_names)),
            selectInput("kallisto_sample_filter_col", "Filter by", choices = c("All transcripts" = "all", kallisto_filter_columns), selected = "all"),
            selectizeInput("kallisto_sample_filter_value", "Select value", choices = NULL, selected = NULL, multiple = FALSE)
          ),
          conditionalPanel(
            "input.kallisto_view_mode == 'matrix'",
            selectInput("kallisto_filter_col", "Filter by", choices = kallisto_filter_columns, selected = "gene_symbol"),
            selectizeInput("kallisto_filter_value", "Select value", choices = NULL, selected = NULL, multiple = FALSE),
            uiOutput("kallisto_build_ui")
          ),
          tags$hr(),
          helpText("Single sample matrix shows transcript-level abundance for one sample. Transcript matrix view shows matching transcripts across all samples.")
        ),
        mainPanel(
          h4("Kallisto Transcript Abundance"),
          uiOutput("kallisto_status_ui"),
          table_widget("kallisto_table")
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
              selected = if ("treatment" %in% comparison_columns) "treatment" else "NA"
            ),
            uiOutput("deg_treatment_ui"),
            uiOutput("deg_control_ui"),
            conditionalPanel(
              "input.deg_subtab == 'DEGs'",
              selectizeInput("deg_gene_query", "Select gene", choices = NULL, selected = NULL, multiple = FALSE),
              uiOutput("deg_convert_ui"),
              selectInput("deg_p_col", "P-value column", choices = c("padj", "pvalue"), selected = "padj"),
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
                selected = "up"
              )
            ),
            conditionalPanel(
              "input.deg_subtab == 'Heatmap'",
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
                selected = "deseq"
              ),
              conditionalPanel(
                "input.heatmap_gene_mode == 'manual'",
                selectizeInput("heatmap_genes", "Genes to plot", choices = NULL, selected = NULL, multiple = TRUE)
              ),
              conditionalPanel(
                "input.heatmap_gene_mode == 'topdeg'",
                selectInput("heatmap_rank_p_col", "Rank by p-value column", choices = c("padj", "pvalue"), selected = "padj"),
                numericInput("heatmap_total_genes", "Total genes", value = 20, min = 2, step = 2),
                checkboxInput("heatmap_equal_split", "Equal numbers of up and down", value = TRUE)
              ),
              selectizeInput(
                "heatmap_annotation_cols",
                "Column annotations",
                choices = comparison_columns,
                selected = if ("treatment" %in% comparison_columns) "treatment" else character(0),
                multiple = TRUE
              ),
              checkboxInput("heatmap_cluster_rows", "Cluster rows", value = TRUE),
              checkboxInput("heatmap_cluster_cols", "Cluster columns", value = TRUE),
              conditionalPanel(
                "!input.heatmap_cluster_cols",
                selectInput(
                  "heatmap_order_col",
                  "Order samples by",
                  choices = comparison_columns,
                  selected = if ("treatment" %in% comparison_columns) "treatment" else comparison_columns[1]
                )
              ),
              checkboxInput("heatmap_scale_rows", "Scale heatmap rows", value = TRUE),
              selectInput(
                "heatmap_theme",
                "Heatmap style",
                choices = c(
                  "Clean publication" = "publication",
                  "Minimal" = "minimal",
                  "Compact" = "compact",
                  "Fully labeled" = "labeled",
                  "High contrast" = "contrast"
                ),
                selected = "publication"
              ),
              selectInput(
                "heatmap_palette",
                "Heatmap colors",
                choices = c(
                  "Use style default" = "theme",
                  "Blue-white-red" = "blue_red",
                  "Viridis" = "viridis",
                  "Magma" = "magma",
                  "Green-white-purple" = "green_purple",
                  "Blue-white-red vivid" = "navy_orange",
                  "Gray-white-red" = "gray_red"
                ),
                selected = "theme"
              ),
              selectInput(
                "heatmap_border_style",
                "Cell borders",
                choices = c("Use style default" = "theme", "None" = "none", "Light grid" = "light"),
                selected = "theme"
              ),
              selectInput(
                "heatmap_font_family",
                "Font family",
                choices = c("Sans" = "sans", "Serif" = "serif", "Mono" = "mono"),
                selected = "sans"
              ),
              selectInput(
                "heatmap_row_font_size",
                "Gene label size",
                choices = c("Small" = 7, "Medium" = 9, "Large" = 11, "XL" = 13),
                selected = 9
              ),
              selectInput(
                "heatmap_col_font_size",
                "Sample label size",
                choices = c("Small" = 8, "Medium" = 10, "Large" = 12, "XL" = 14),
                selected = 10
              ),
              selectInput(
                "heatmap_col_angle",
                "Sample label angle",
                choices = c("0" = "0", "45" = "45", "90" = "90", "315" = "315"),
                selected = "45"
              ),
              checkboxInput("heatmap_show_row_names", "Show gene labels", value = TRUE),
              checkboxInput("heatmap_show_col_names", "Show sample labels", value = TRUE),
              checkboxInput("heatmap_show_annotation_names", "Show annotation track labels", value = FALSE)
            ),
            tags$hr(),
            helpText("Shows DEGs, PCA, and a custom heatmap for an available treatment vs control comparison.")
          ),
          mainPanel(
            uiOutput("deg_status_ui"),
            tabsetPanel(
              id = "deg_subtab",
              tabPanel(
                "DEGs",
                tags$div(
                  style = "max-height: 520px; overflow-y: auto; overflow-x: auto;",
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
              ),
              tabPanel(
                "PCA",
                uiOutput("deg_pca_ui")
              ),
              tabPanel(
                "Heatmap",
                textInput("heatmap_filename", "Filename", value = ""),
                actionButton("save_heatmap_btn", "Save heatmap"),
                plotOutput("deg_heatmap_plot", height = "700px")
              )
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
              selected = if ("treatment" %in% comparison_columns) "treatment" else "NA"
            ),
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
      table.dataTable tbody td, table.dataTable thead th {
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        max-width: 210px;
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
      }
      .dataTables_wrapper {
        width: 100%;
        overflow-x: auto;
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

  if (DT_AVAILABLE) {
    output$star_summary_table <- DT::renderDT({
      req(star_available)
      df <- order_sample_columns(star_summary_df, annotation_col = value_or(input$star_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("metric"))
      simple_dt(df, page_length = 25)
    })
  } else {
    output$star_summary_table <- renderTable({
      req(star_available)
      df <- order_sample_columns(star_summary_df, annotation_col = value_or(input$star_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("metric"))
      rownames(df) <- df$metric
      df$metric <- NULL
      format_numeric_commas(df)
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
  }

  output$star_sample_table <- renderTable({
    req(star_available)
    req(input$star_sample)
    sample_col <- input$star_sample
    req(sample_col %in% colnames(star_summary_df))
    format_numeric_commas(data.frame(
      Metric = star_summary_df$metric,
      Value = star_summary_df[[sample_col]],
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$featurecounts_qc_status_ui <- renderUI({
    if (!featurecounts_available) {
      return(status_box("FeatureCounts summary has not been generated yet.", "warning"))
    }
    NULL
  })

  if (DT_AVAILABLE) {
    output$featurecounts_summary_table <- DT::renderDT({
      req(featurecounts_available)
      df <- featurecounts_summary_df
      df <- order_sample_columns(df, annotation_col = value_or(input$featurecounts_qc_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("Status"))
      simple_dt(df, page_length = 25)
    })
  } else {
    output$featurecounts_summary_table <- renderTable({
      req(featurecounts_available)
      df <- featurecounts_summary_df
      df <- order_sample_columns(df, annotation_col = value_or(input$featurecounts_qc_sample_sort_col, "__alpha__"), design_df = design_df, id_cols = c("Status"))
      rownames(df) <- df$Status
      df$Status <- NULL
      format_numeric_commas(df)
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
  }

  output$featurecounts_convert_ui <- renderUI({
    if (count_matrix_available && "Geneid" %in% colnames(count_matrix_nonzero_df) && looks_like_gene_id(count_matrix_nonzero_df$Geneid)) {
      actionButton("featurecounts_convert_btn", "Convert Ensembl IDs to gene names")
    }
  })

  featurecounts_display_df <- reactive({
    df <- count_matrix_nonzero_df
    if ("Geneid" %in% colnames(df) && looks_like_gene_id(df$Geneid) && (value_or(input$featurecounts_convert_btn, 0)) %% 2 == 1) {
      species <- detect_species_from_ids(df$Geneid)
      map_df <- get_gtf_map(species)
      conv <- convert_gene_labels(df$Geneid, map_df)
      df$Geneid <- conv$values
      df <- aggregate_display_matrix(df, "Geneid")
    }
    df
  })

  output$featurecounts_status_ui <- renderUI({
    if (!count_matrix_available) {
      return(status_box("Raw count matrix has not been generated yet.", "warning"))
    }
    if ("Geneid" %in% colnames(count_matrix_nonzero_df) && looks_like_gene_id(count_matrix_nonzero_df$Geneid) && (value_or(input$featurecounts_convert_btn, 0)) %% 2 == 1) {
      species <- detect_species_from_ids(count_matrix_nonzero_df$Geneid)
      map_df <- get_gtf_map(species)
      conv <- convert_gene_labels(head(count_matrix_nonzero_df$Geneid, 100), map_df)
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf("FeatureCounts display converted from %s to %s using %s GTF.", conv$from, conv$to, value_or(species, "detected"))
      )
    } else {
      NULL
    }
  })

  observe({
    if (!count_matrix_available) {
      updateSelectizeInput(session, "gene_query", choices = character(0), selected = character(0), server = TRUE)
      updateSelectInput(session, "raw_counts_sort_col", choices = c("Geneid"), selected = "Geneid")
      return()
    }
    updateSelectInput(
      session,
      "raw_counts_sort_col",
      choices = colnames(featurecounts_display_df()),
      selected = if ("Geneid" %in% colnames(featurecounts_display_df())) "Geneid" else colnames(featurecounts_display_df())[1]
    )
    updateSelectizeInput(
      session,
      "gene_query",
      choices = featurecounts_display_df()$Geneid,
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
    })
  } else {
    output$gene_search_table <- renderTable({
      df <- raw_counts_table_df()
      rownames(df) <- make.unique(as.character(df$Geneid))
      df$Geneid <- NULL
      format_numeric_commas(df)
    }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
  }

  if (rsem_available) {
    rsem_matrix_cache <- reactiveValues()
    rsem_display_cache <- reactiveValues()

    get_rsem_matrix <- reactive({
      req(input$rsem_type, input$rsem_metric)
      cache_key <- paste(input$rsem_type, input$rsem_metric, sep = "::")
      cached <- rsem_matrix_cache[[cache_key]]
      if (is.null(cached)) {
        cached <- build_rsem_matrix(samples, input$rsem_type, input$rsem_metric)
        rsem_matrix_cache[[cache_key]] <- cached
      }
      cached
    })

    output$rsem_convert_ui <- renderUI({
      raw_df <- get_rsem_matrix()
      if (!is.null(raw_df) && "gene_id" %in% colnames(raw_df) && looks_like_gene_id(raw_df$gene_id)) {
        actionButton("rsem_convert_btn", "Convert Ensembl IDs to gene names")
      }
    })

    rsem_display_df <- reactive({
      raw_df <- get_rsem_matrix()
      req(!is.null(raw_df))
      cache_key <- paste(
        input$rsem_type,
        input$rsem_metric,
        if ("gene_id" %in% colnames(raw_df) && looks_like_gene_id(raw_df$gene_id)) value_or(input$rsem_convert_btn, 0) %% 2 else 0,
        sep = "::"
      )
      cached <- rsem_display_cache[[cache_key]]
      if (!is.null(cached)) {
        return(cached)
      }
      df <- raw_df
      if ("gene_id" %in% colnames(df) && looks_like_gene_id(df$gene_id) && (value_or(input$rsem_convert_btn, 0)) %% 2 == 1) {
        species <- detect_species_from_ids(df$gene_id)
        map_df <- get_gtf_map(species)
        conv <- convert_gene_labels(df$gene_id, map_df)
        df$gene_id <- conv$values
        if (identical(input$rsem_type, "genes")) {
          df <- aggregate_display_matrix(df, "gene_id")
        }
      }
      rsem_display_cache[[cache_key]] <- df
      df
    })

    output$rsem_status_ui <- renderUI({
      raw_df <- get_rsem_matrix()
      if (is.null(raw_df)) {
        return(status_box("RSEM results have not been generated yet.", "warning"))
      }
      if ("gene_id" %in% colnames(raw_df) && looks_like_gene_id(raw_df$gene_id) && (value_or(input$rsem_convert_btn, 0)) %% 2 == 1) {
        species <- detect_species_from_ids(raw_df$gene_id)
        map_df <- get_gtf_map(species)
        conv <- convert_gene_labels(head(raw_df$gene_id, 100), map_df)
        tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
          if (isTRUE(conv$changed)) {
            sprintf("RSEM display converted from %s to %s using %s GTF (%s IDs mapped).", conv$from, conv$to, value_or(species, "detected"), conv$mapped)
          } else {
            sprintf("RSEM conversion requested, but no %s labels were found in the %s GTF map.", conv$from, value_or(species, "detected"))
          }
        )
      } else {
        NULL
      }
    })

    observe({
      req(input$rsem_type, input$rsem_metric)
      df <- rsem_display_df()
      req(!is.null(df))
      id_col <- if (identical(input$rsem_type, "genes")) "gene_id" else "transcript_id"
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
      id_col <- if (identical(input$rsem_type, "genes")) "gene_id" else "transcript_id"
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
        df <- rsem_table_df()
        simple_dt(df, page_length = 50)
      })
    } else {
      output$rsem_table <- renderTable({
        df <- rsem_table_df()
        id_col <- if (identical(input$rsem_type, "genes")) "gene_id" else "transcript_id"
        rownames(df) <- make.unique(as.character(df[[id_col]]))
        df[[id_col]] <- NULL
        format_numeric_commas(df)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
    }
  } else {
    output$rsem_status_ui <- renderUI({
      status_box("RSEM results have not been generated yet.", "warning")
    })
    output$rsem_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable(NULL)
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
      selectInput("deseq_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals)
    })

    output$deseq_control_ui <- renderUI({
      vals <- get_deseq_values()
      selectInput("deseq_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals)
    })

    deseq_counts_raw_df <- reactive({
      req(input$deseq_compare_col)
      if (identical(input$deseq_compare_col, "NA")) return(NULL)
      req(input$deseq_treatment, input$deseq_control)
      path <- normalized_counts_path(input$deseq_treatment, input$deseq_control)
      df <- read_normalized_counts(path)
      if (is.null(df)) return(NULL)
      keep_samples <- design_df$sample[design_df[[input$deseq_compare_col]] %in% c(input$deseq_treatment, input$deseq_control)]
      keep_cols <- c("gene_label", intersect(keep_samples, colnames(df)))
      df[, keep_cols, drop = FALSE]
    })

    output$deseq_convert_ui <- renderUI({
      df <- deseq_counts_raw_df()
      if (!is.null(df) && "gene_label" %in% colnames(df) && looks_like_gene_id(df$gene_label)) {
        actionButton("deseq_convert_btn", "Convert Ensembl IDs to gene names")
      }
    })

    deseq_counts_df <- reactive({
      df <- deseq_counts_raw_df()
      if (is.null(df)) return(NULL)
      if ("gene_label" %in% colnames(df) && looks_like_gene_id(df$gene_label) && (value_or(input$deseq_convert_btn, 0)) %% 2 == 1) {
        species <- detect_species_from_ids(df$gene_label)
        map_df <- get_gtf_map(species)
        conv <- convert_gene_labels(df$gene_label, map_df)
        df$gene_label <- conv$values
        df <- aggregate_display_matrix(df, "gene_label")
      }
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
      path <- normalized_counts_path(input$deseq_treatment, input$deseq_control)
      if (!file.exists(path)) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;",
          sprintf("This comparison has not been tested: %s vs %s.", input$deseq_treatment, input$deseq_control)
        ))
      }
      raw_df <- read_normalized_counts(path)
      shown_df <- deseq_counts_df()
      raw_n <- if (is.null(raw_df)) 0 else nrow(raw_df)
      shown_n <- if (is.null(shown_df)) 0 else nrow(shown_df)
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf("Showing normalized counts for %s vs %s using %s. File has %s genes; display has %s genes. DESeq2.R writes normalized counts after keep <- rowSums(counts(dds)) >= 10.", input$deseq_treatment, input$deseq_control, input$deseq_compare_col, raw_n, shown_n)
      )
    })

    observe({
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
        simple_dt(deseq_counts_table_df(), page_length = 50)
      })
    } else {
      output$deseq_counts_table <- renderTable({
        show_df <- deseq_counts_table_df()
        rownames(show_df) <- make.unique(as.character(show_df$gene_label))
        show_df$gene_label <- NULL
        format_numeric_commas(show_df)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = TRUE)
    }

    get_deg_values <- reactive({
      if (is.null(input$deg_compare_col) || identical(input$deg_compare_col, "NA") || !nzchar(input$deg_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$deg_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$deg_treatment_ui <- renderUI({
      vals <- get_deg_values()
      selectInput("deg_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals)
    })

    output$deg_control_ui <- renderUI({
      vals <- get_deg_values()
      selectInput("deg_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals)
    })

    deg_raw_df <- reactive({
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      req(input$deg_treatment, input$deg_control)
      read_deg_table(deg_table_path(input$deg_treatment, input$deg_control))
    })

    output$deg_convert_ui <- renderUI({
      df <- deg_raw_df()
      if (!is.null(df) && "gene_label" %in% colnames(df) && looks_like_gene_id(df$gene_label)) {
        actionButton("deg_convert_btn", "Convert Ensembl IDs to gene names")
      }
    })

    deg_df <- reactive({
      df <- deg_raw_df()
      if (is.null(df)) return(NULL)
      if ("gene_label" %in% colnames(df) && looks_like_gene_id(df$gene_label) && (value_or(input$deg_convert_btn, 0)) %% 2 == 1) {
        species <- detect_species_from_ids(df$gene_label)
        map_df <- get_gtf_map(species)
        conv <- convert_gene_labels(df$gene_label, map_df)
        df$gene_label <- conv$values
      }
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
      path <- deg_table_path(input$deg_treatment, input$deg_control)
      if (!file.exists(path)) {
        return(tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fdeeee; border: 1px solid #efb0b0; border-radius: 6px;",
          sprintf("This comparison has not been tested: %s vs %s.", input$deg_treatment, input$deg_control)
        ))
      }
      tags$div(
        style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
        sprintf("Showing DEG results for %s vs %s using %s.", input$deg_treatment, input$deg_control, input$deg_compare_col)
      )
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

    heatmap_rsem_df <- reactive({
      req(rsem_available)
      req(input$deg_compare_col)
      if (identical(input$deg_compare_col, "NA")) return(NULL)
      req(input$deg_treatment, input$deg_control)
      df <- build_rsem_matrix(samples, data_type = "genes", metric = "TPM")
      if (is.null(df)) return(NULL)
      keep_samples <- design_df$sample[design_df[[input$deg_compare_col]] %in% c(input$deg_treatment, input$deg_control)]
      keep_cols <- c("gene_id", intersect(keep_samples, colnames(df)))
      df <- df[, keep_cols, drop = FALSE]
      colnames(df)[colnames(df) == "gene_id"] <- "gene_label"
      df
    })

    heatmap_source_df <- reactive({
      if (identical(value_or(input$heatmap_source, "deseq"), "rsem")) {
        heatmap_rsem_df()
      } else {
        deg_norm_counts_df()
      }
    })

    observe({
      df <- deg_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "deg_gene_query", choices = character(0), selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session, "deg_gene_query", choices = df$gene_label, selected = character(0), server = TRUE)
      }
    })

    observe({
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
      df <- heatmap_source_df()
      if (is.null(df)) {
        updateSelectizeInput(session, "heatmap_genes", choices = character(0), selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session, "heatmap_genes", choices = df$gene_label, selected = isolate(input$heatmap_genes), server = TRUE)
      }
    })

    heatmap_selected_genes <- reactive({
      if (identical(input$heatmap_gene_mode, "manual")) {
        return(value_or(input$heatmap_genes, character(0)))
      }
      df <- deg_df()
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
      })
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
        border_color = heatmap_border_color(input$heatmap_border_style, heatmap_theme),
        color = heatmap_palette_colors(input$heatmap_palette, heatmap_theme)
      )
    })

    observeEvent(input$save_heatmap_btn, {
      hp <- heatmap_plot_state()
      default_name <- sprintf(
        "heatmap_%s_%s_vs_%s.png",
        value_or(input$deg_compare_col, "comparison"),
        value_or(input$deg_treatment, "treatment"),
        value_or(input$deg_control, "control")
      )
      filename <- sanitize_filename(input$heatmap_filename)
      if (!nzchar(filename)) filename <- default_name
      if (!grepl("\\.png$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".png")
      out_path <- file.path(deseq2_dir, filename)
      png(out_path, width = 2400, height = 1900, res = 220)
      heatmap_theme <- value_or(input$heatmap_theme, "publication")
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
        border_color = heatmap_border_color(input$heatmap_border_style, heatmap_theme),
        color = heatmap_palette_colors(input$heatmap_palette, heatmap_theme)
      )
      dev.off()
      showNotification(sprintf("Saved heatmap to %s", out_path), type = "message", duration = 6)
    })
  } else {
    empty_compare_input <- function(id, label) {
      renderUI(selectInput(id, label, choices = character(0), selected = character(0)))
    }
    output$deseq_treatment_ui <- empty_compare_input("deseq_treatment", "Treatment")
    output$deseq_control_ui <- empty_compare_input("deseq_control", "Control")
    output$deseq_status_ui <- renderUI({
      status_box("DESeq2 normalized counts have not been generated yet.", "warning")
    })
    output$deseq_counts_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })

    output$deg_treatment_ui <- empty_compare_input("deg_treatment", "Treatment")
    output$deg_control_ui <- empty_compare_input("deg_control", "Control")
    output$deg_status_ui <- renderUI({
      status_box("Differential expression results have not been generated yet.", "warning")
    })
    output$deg_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })
    output$deg_pca_ui <- renderUI({
      tags$p("PCA has not been generated yet.")
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
      updateSelectizeInput(session, "heatmap_genes", choices = character(0), selected = character(0), server = TRUE)
    })
  }

  if (gseapy_available) {
    get_gsea_values <- reactive({
      if (is.null(input$gsea_compare_col) || identical(input$gsea_compare_col, "NA") || !nzchar(input$gsea_compare_col)) {
        return(character(0))
      }
      vals <- unique(as.character(design_df[[input$gsea_compare_col]]))
      vals[!is.na(vals) & nzchar(vals)]
    })

    output$gsea_treatment_ui <- renderUI({
      vals <- get_gsea_values()
      selectInput("gsea_treatment", "Treatment", choices = vals, selected = if (length(vals) > 1) vals[2] else vals)
    })

    output$gsea_control_ui <- renderUI({
      vals <- get_gsea_values()
      selectInput("gsea_control", "Control", choices = vals, selected = if (length(vals) > 0) vals[1] else vals)
    })

    gsea_comp_dir <- reactive({
      req(input$gsea_compare_col)
      if (identical(input$gsea_compare_col, "NA")) return(NULL)
      req(input$gsea_treatment, input$gsea_control)
      gseapy_comparison_dir(input$gsea_treatment, input$gsea_control)
    })

    output$gsea_collection_ui <- renderUI({
      comp_dir <- gsea_comp_dir()
      collections <- list_gseapy_collections(comp_dir)
      selectInput("gsea_collection", "Pathway collection", choices = collections, selected = if (length(collections)) collections[1] else character(0))
    })

    gsea_report_df <- reactive({
      comp_dir <- gsea_comp_dir()
      req(!is.null(comp_dir), input$gsea_collection)
      read_gseapy_report(comp_dir, input$gsea_collection)
    })

    output$gsea_pathway_ui <- renderUI({
      df <- gsea_report_df()
      pathways <- if (is.null(df)) character(0) else df$Term
      selectizeInput("gsea_pathway", "Pathway", choices = pathways, selected = if (length(pathways)) pathways[1] else character(0), multiple = FALSE)
    })

    output$gsea_status_ui <- renderUI({
      if (!gseapy_available) {
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
      })

      output$gsea_all_pathways_table <- DT::renderDT({
        df <- gsea_report_df()
        req(!is.null(df))
        simple_dt(df, page_length = 25)
      })
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
        tags$div(
          tags$h5(title),
          if (!is.null(rendered_rel)) {
            tags$img(src = rendered_rel, style = "width: 100%; max-width: 100%; border: 1px solid #ddd; margin-bottom: 12px;")
          } else {
            status_box("This server could not render the PDF to an image. Open the PDF link below.", "warning")
          },
          tags$p(tags$a(href = rel, target = "_blank", "Open original PDF"))
        )
      }
      tags$div(
        plot_block("Enrichment Plot", pdf_rel),
        plot_block("Pathway Heatmap", heatmap_rel)
      )
    })
  } else {
    output$gsea_treatment_ui <- renderUI(selectInput("gsea_treatment", "Treatment", choices = character(0), selected = character(0)))
    output$gsea_control_ui <- renderUI(selectInput("gsea_control", "Control", choices = character(0), selected = character(0)))
    output$gsea_collection_ui <- renderUI(selectInput("gsea_collection", "Pathway collection", choices = character(0), selected = character(0)))
    output$gsea_pathway_ui <- renderUI(selectizeInput("gsea_pathway", "Pathway", choices = character(0), selected = character(0), multiple = FALSE))
    output$gsea_status_ui <- renderUI({
      status_box("GSEA results have not been generated yet.", "warning")
    })
    output$gsea_pathway_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })
    output$gsea_all_pathways_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable({ NULL })
    output$gsea_summary_plots_ui <- renderUI({ tags$p("GSEA summary plots are not available yet.") })
    output$gsea_pathway_plots_ui <- renderUI({ tags$p("No pathway-specific GSEA plots are available yet.") })
  }

  if (kallisto_available) {
    kallisto_matrix_cache <- reactiveVal(read_kallisto_matrix())
    kallisto_filter_cache <- reactiveValues()
    kallisto_building <- reactiveVal(FALSE)

    output$kallisto_build_ui <- renderUI({
      if (is.null(kallisto_matrix_cache())) {
        actionButton("build_kallisto_matrix_btn", "Build transcript matrix")
      } else {
        tags$span("Using saved transcript matrix.")
      }
    })

    observeEvent(input$build_kallisto_matrix_btn, {
      if (isTRUE(kallisto_building()) || !is.null(kallisto_matrix_cache())) return()
      kallisto_building(TRUE)
      showModal(modalDialog("Building transcript matrix...", footer = NULL, easyClose = FALSE))
      built <- build_kallisto_matrix(kallisto_sample_names)
      if (!is.null(built)) write_kallisto_matrix(built)
      kallisto_matrix_cache(built)
      kallisto_building(FALSE)
      removeModal()
    })

    output$kallisto_status_ui <- renderUI({
      if (!kallisto_available) {
        return(status_box("Kallisto transcript abundance has not been generated yet.", "warning"))
      }
      if (identical(input$kallisto_view_mode, "matrix") && !is.null(kallisto_matrix_cache())) {
        return(status_box(sprintf("Using Kallisto transcript matrix: %s", kallisto_matrix_path), "info"))
      }
      if (identical(input$kallisto_view_mode, "matrix") && is.null(kallisto_matrix_cache())) {
        tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #fff7e6; border: 1px solid #f0c36d; border-radius: 6px;",
          "Click 'Build transcript matrix' to load the joined transcript table across all samples."
        )
      } else if (isTRUE(kallisto_building()) && identical(input$kallisto_view_mode, "matrix")) {
        tags$div(
          style = "margin-bottom: 12px; padding: 10px 12px; background: #eef5ff; border: 1px solid #b9d0f5; border-radius: 6px;",
          "Building transcript matrix..."
        )
      } else {
        NULL
      }
    })

    observe({
      req(identical(input$main_tabs, "Counts"))
      req(identical(input$counts_subtab, "Kallisto"))
      req(identical(input$kallisto_view_mode, "sample"))
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
      req(input$kallisto_filter_col)
      col_key <- input$kallisto_filter_col
      vals <- kallisto_filter_cache[[col_key]]
      if (is.null(vals)) {
        vals <- collect_kallisto_filter_values(samples, col_key)
        kallisto_filter_cache[[col_key]] <- vals
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
        kallisto_matrix_df <- kallisto_matrix_cache()
        req(!is.null(kallisto_matrix_df))
        if (is.null(input$kallisto_filter_value) || !nzchar(trimws(input$kallisto_filter_value))) {
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
        simple_dt(kallisto_table_df(), page_length = 50)
      })
    } else {
      output$kallisto_table <- renderTable({
        format_numeric_commas(kallisto_table_df())
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)
    }
  } else {
    output$kallisto_status_ui <- renderUI({
      status_box("Kallisto transcript abundance has not been generated yet.", "warning")
    })
    output$kallisto_table <- if (DT_AVAILABLE) DT::renderDT(NULL) else renderTable(NULL)
    observe({
      updateSelectizeInput(session, "kallisto_sample_filter_value", choices = character(0), selected = character(0), server = TRUE)
      updateSelectizeInput(session, "kallisto_filter_value", choices = character(0), selected = character(0), server = TRUE)
    })
  }
}

shinyApp(ui, server)
