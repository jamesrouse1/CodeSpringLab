args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
root <- tempfile("csl-native-comparison-test-")
data_dir <- file.path(root, "data")
dir.create(file.path(data_dir, "deseq2"), recursive = TRUE)
dir.create(file.path(data_dir, "gseapy", "Treated_vs_Control", "Hallmark"), recursive = TRUE)
dir.create(file.path(data_dir, "manifest"), recursive = TRUE)
dir.create(file.path(data_dir, "star", "S1"), recursive = TRUE)

design_path <- file.path(data_dir, "manifest", "design_matrix.txt")
write.table(
  data.frame(
    sample = c("S1", "S2"),
    treatment = c("Control", "Treated"),
    filename = c("S1", "S2")
  ),
  design_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
writeLines(
  "DESCRIPTION\tS1\tS2",
  file.path(data_dir, "deseq2", "normalized_counts_Treated_vs_Control(ref).txt")
)
writeLines(
  "gene\tbaseMean\tlog2FoldChange\tpvalue\tpadj",
  file.path(data_dir, "deseq2", "DEG_Treated_vs_Control(ref).txt")
)
writeLines(
  "Term,NES,NOM p-val,FDR q-val",
  file.path(data_dir, "gseapy", "Treated_vs_Control", "Hallmark", "report.gseapy.Hallmark.csv")
)
writeLines("synthetic STAR log", file.path(data_dir, "star", "S1", "S1_Log.final.out"))

config_path <- file.path(root, "config.R")
writeLines(c(
  sprintf("project_name <- %s", deparse("test")),
  sprintf("results_root <- %s", deparse(root)),
  sprintf("data_dir <- %s", deparse(data_dir)),
  sprintf("design_matrix_path <- %s", deparse(design_path)),
  sprintf("app_dir <- %s", deparse(file.path(repo_root, "scripts_DoNotTouch", "Shiny")))
), config_path)

old_config <- Sys.getenv("RNASEQ_SHINY_CONFIG", unset = NA_character_)
on.exit({
  if (is.na(old_config)) Sys.unsetenv("RNASEQ_SHINY_CONFIG") else Sys.setenv(RNASEQ_SHINY_CONFIG = old_config)
  unlink(root, recursive = TRUE, force = TRUE)
}, add = TRUE)
Sys.setenv(RNASEQ_SHINY_CONFIG = config_path)

app_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, "scripts_DoNotTouch", "Shiny", "app_server.R"), envir = app_env)

stopifnot(
  NROW(app_env$completed_deseq_catalog) == 1L,
  identical(app_env$completed_deseq_catalog$compare_col[[1]], "treatment"),
  identical(app_env$completed_deseq_catalog$treatment[[1]], "Treated"),
  identical(app_env$completed_deseq_catalog$control[[1]], "Control"),
  NROW(app_env$completed_gsea_catalog) == 1L,
  identical(app_env$completed_gsea_catalog$compare_col[[1]], "treatment"),
  identical(app_env$completed_gsea_catalog$treatment[[1]], "Treated"),
  identical(app_env$completed_gsea_catalog$control[[1]], "Control"),
  all(c("Tool", "Sample", "File", "Size", "Modified", "Absolute path", "Copy path") %in% names(app_env$rna_result_files)),
  any(app_env$rna_result_files$Tool == "STAR" & app_env$rna_result_files$Sample == "S1" & app_env$rna_result_files$File == "S1_Log.final.out"),
  !any(grepl("/", app_env$rna_result_files$File, fixed = TRUE))
)

ui_text <- paste(as.character(app_env$ui), collapse = "\n")
server_text <- paste(deparse(body(app_env$server)), collapse = "\n")
stopifnot(
  grepl("download_selected_rna_file", ui_text, fixed = TRUE),
  grepl("rna-trigger-file-download", ui_text, fixed = TRUE),
  grepl("validated_rna_result_file", server_text, fixed = TRUE),
  grepl("confirm_delete_rna_file", server_text, fixed = TRUE),
  grepl("This cannot be undone", server_text, fixed = TRUE)
)

cat("Native completed RNA-seq comparison fixture passed.\n")
