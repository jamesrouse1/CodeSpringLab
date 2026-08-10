#!/usr/bin/env Rscript

# CodeSpringLab single-cell RNA-seq workflow (Seurat engine)
#
# Usage:
#   Rscript scrna_pipeline_seurat.R <samples.tsv> <out_dir> <params.tsv> [inspect|qc|preprocess|cluster|annotate|score|differential|all]
#
# samples.tsv (required columns): sample_id, input_path
# Optional sample columns (for example condition and batch) are retained as
# per-cell metadata.  input_path may be a filtered 10x matrix directory or a
# .rds Seurat object. AnnData .h5ad inputs are handled natively by Scanpy.
# params.tsv contains key/value pairs; comments and unknown keys are allowed.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) stop("Usage: scrna_pipeline_seurat.R <samples.tsv> <out_dir> <params.tsv> [inspect|qc|preprocess|cluster|annotate|score|differential|all]")

samples_path <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
params_path <- normalizePath(args[[3]], mustWork = TRUE)
stage <- tolower(if (length(args) >= 4L) args[[4]] else "all")
if (!stage %in% c("inspect", "qc", "preprocess", "cluster", "annotate", "score", "differential", "all")) stop("Unknown scRNA stage: ", stage)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is not installed: ", pkg)
}
for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork")) need_pkg(pkg)
suppressPackageStartupMessages(library(Seurat))
# A submitted job already has the requested resources. Keeping Seurat's future
# plan sequential here prevents large pre-existing RDS objects from being
# serialized to nested workers, which otherwise triggers future.globals.maxSize
# failures and can duplicate the full expression matrix in memory.
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
  # SCTransform still validates globals under a sequential plan. Large RDS
  # inputs can legitimately exceed the interactive 500 MiB default without
  # transferring data to another worker, so raise only this job-local ceiling.
  options(future.globals.maxSize = 8 * 1024^3)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x

# The jpplot gradient is reserved for continuous expression and score values.
jpplot_colors <- c("#90C3DD", "#C2E4EF", "#ECF7E1", "#FEF4AF", "#FDD484", "#FBA25B", "#F0653F", "#D42D26", "#A50026")
categorical_palette <- function(values) {
  levels <- sort(unique(as.character(values)))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  if (!length(levels)) return(character(0))
  colors <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#6A3D9A",
    "#56B4E9", "#B15928", "#E31A1C", "#1B9E77", "#7570B3", "#E7298A",
    "#66A61E", "#A6761D", "#1F78B4", "#FF7F00", "#33A02C", "#984EA3",
    "#A65628", "#F781BF", "#17BECF", "#BCBD22", "#8C564B", "#4D4D4D"
  )
  if (length(levels) > length(colors)) colors <- c(colors, grDevices::hcl.colors(length(levels) - length(colors), "Dynamic"))
  colors <- colors[seq_along(levels)]
  stats::setNames(colors, levels)
}
categorical_dimplot <- function(obj, group.by, ...) {
  Seurat::DimPlot(obj, group.by = group.by, cols = categorical_palette(obj[[group.by]][, 1]), ...)
}
jpplot_point_layer <- function(plot, color = jpplot_colors[[7]]) {
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, "GeomPoint")) plot$layers[[i]]$aes_params$colour <- color
  }
  plot
}

read_delim_safe <- function(path) {
  tryCatch(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = ""),
           error = function(e) stop("Could not read ", path, ": ", conditionMessage(e)))
}

read_params <- function(path) {
  x <- read_delim_safe(path)
  if (!all(c("key", "value") %in% names(x))) stop("Parameter file must contain key and value columns.")
  values <- as.list(as.character(x$value)); names(values) <- as.character(x$key)
  get <- function(key, default = "") trimws(as.character(values[[key]] %||% default))
  get_bool <- function(key, default = FALSE) tolower(get(key, if (default) "true" else "false")) %in% c("1", "true", "t", "yes", "y")
  list(
    normalization = tolower(get("normalization", "sct")),
    integration = tolower(get("integration", "auto")),
    batch_column = get("batch_column", "batch"),
    cluster_resolution = suppressWarnings(as.numeric(get("cluster_resolution", "0.6"))),
    n_pcs = suppressWarnings(as.integer(get("n_pcs", "30"))),
    min_features = suppressWarnings(as.integer(get("min_features", "200"))),
    min_counts = suppressWarnings(as.integer(get("min_counts", "0"))),
    max_features = suppressWarnings(as.integer(get("max_features", "0"))),
    max_percent_mt = suppressWarnings(as.numeric(get("max_percent_mt", "20"))),
    min_cells_per_gene = suppressWarnings(as.integer(get("min_cells_per_gene", "3"))),
    doublet_method = tolower(get("doublet_method", "none")),
    doublet_rate = suppressWarnings(as.numeric(get("doublet_rate", "0.05"))),
    remove_doublets = get_bool("remove_doublets", TRUE),
    marker_file = get("marker_file", ""),
    celltype_file = get("celltype_file", ""),
    annotation_name = get("annotation_name", "cell_type"),
    signature_file = get("signature_file", ""),
    de_group_column = get("de_group_column", "condition"),
    de_reference = get("de_reference", ""),
    de_comparison = get("de_comparison", ""),
    de_annotation_column = get("de_annotation_column", ""),
    de_annotation_values = Filter(nzchar, strsplit(get("de_annotation_values", get("de_annotation_value", "all")), "||", fixed = TRUE)[[1]]),
    de_method = tolower(get("de_method", "both")),
    de_covariates = Filter(nzchar, trimws(strsplit(get("de_covariates", ""), ",", fixed = TRUE)[[1]])),
    harmony_theta = suppressWarnings(as.numeric(get("harmony_theta", "2"))),
    harmony_lambda = suppressWarnings(as.numeric(get("harmony_lambda", "1"))),
    harmony_max_iter = suppressWarnings(as.integer(get("harmony_max_iter", "20"))),
    seed = suppressWarnings(as.integer(get("seed", "1234")))
  )
}

params <- read_params(params_path)
if (!params$normalization %in% c("sct", "lognormalize")) stop("normalization must be SCT or LogNormalize")
if (!params$integration %in% c("auto", "none", "rpca", "cca", "harmony")) stop("integration must be auto, none, rpca, cca, or harmony")
if (!is.finite(params$cluster_resolution) || params$cluster_resolution <= 0) params$cluster_resolution <- 0.6
if (!is.finite(params$n_pcs) || params$n_pcs < 5) params$n_pcs <- 30L
if (!is.finite(params$min_features) || params$min_features < 0) params$min_features <- 200L
if (!is.finite(params$min_counts) || params$min_counts < 0) params$min_counts <- 0L
if (!is.finite(params$max_features) || params$max_features < 0) params$max_features <- 0L
if (!is.finite(params$max_percent_mt) || params$max_percent_mt < 0) params$max_percent_mt <- 20
if (!is.finite(params$min_cells_per_gene) || params$min_cells_per_gene < 1) params$min_cells_per_gene <- 3L
if (!params$doublet_method %in% c("auto", "none", "scdblfinder")) stop("doublet_method must be auto, none, or scDblFinder")
if (!is.finite(params$doublet_rate) || params$doublet_rate < 0 || params$doublet_rate >= 1) params$doublet_rate <- 0.05
if (!is.finite(params$harmony_theta) || params$harmony_theta < 0) params$harmony_theta <- 2
if (!is.finite(params$harmony_lambda) || params$harmony_lambda <= 0) params$harmony_lambda <- 1
if (!is.finite(params$harmony_max_iter) || params$harmony_max_iter < 1) params$harmony_max_iter <- 20L
if (!is.finite(params$seed)) params$seed <- 1234L
set.seed(params$seed)

samples <- read_delim_safe(samples_path)
required <- c("sample_id", "input_path")
if (!all(required %in% names(samples))) stop("Sample manifest requires columns: ", paste(required, collapse = ", "))
samples$sample_id <- trimws(as.character(samples$sample_id))
samples$input_path <- trimws(as.character(samples$input_path))
if (!NROW(samples)) stop("No sample rows were supplied.")
if (any(!nzchar(samples$sample_id))) stop("Every sample manifest row needs a sample_id.")
if (any(!nzchar(samples$input_path))) stop("Every sample manifest row needs an input_path.")
if (anyDuplicated(samples$sample_id)) stop("sample_id values must be unique in the input manifest.")
samples$input_path <- normalizePath(path.expand(as.character(samples$input_path)), winslash = "/", mustWork = FALSE)
missing <- samples$input_path[!file.exists(samples$input_path)]
if (length(missing)) stop("Missing unreadable input path(s): ", paste(missing, collapse = ", "))

tables_dir <- file.path(out_dir, "tables"); figures_dir <- file.path(out_dir, "figures"); objects_dir <- file.path(out_dir, "objects"); checkpoints_dir <- file.path(out_dir, "checkpoints")
for (d in c(tables_dir, figures_dir, objects_dir, checkpoints_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
utils::write.table(samples, file.path(tables_dir, "sample_manifest_used.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
stage_marker <- function(name) writeLines("complete", file.path(out_dir, paste0("_STAGE_", toupper(name), "_COMPLETE")))
checkpoint_path <- function(name) file.path(checkpoints_dir, paste0(name, "_seurat.rds"))
require_checkpoint <- function(name, prior) {
  path <- checkpoint_path(name)
  if (!file.exists(path)) stop("The ", prior, " stage has not completed. Run ", prior, " before this stage.")
  readRDS(path)
}

input_kind <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.rds$", lower)) return("seurat_rds")
  if (grepl("\\.h5ad$", lower)) return("scanpy_h5ad")
  if (dir.exists(path)) return("filtered_10x_matrix")
  stop("Unsupported input type: ", path, ". Use a Seurat .rds or filtered 10x matrix directory.")
}

read_10x_matrix <- function(path) {
  value <- tryCatch(Seurat::Read10X(path), error = function(e) NULL)
  if (!is.null(value)) return(value)
  pick <- function(stems) {
    hits <- unlist(lapply(stems, function(stem) c(file.path(path, stem), file.path(path, paste0(stem, ".gz")))), use.names = FALSE)
    hit <- hits[file.exists(hits)]
    if (!length(hit)) return("")
    hit[[1]]
  }
  matrix_file <- pick(c("matrix.mtx")); feature_file <- pick(c("features.tsv", "genes.tsv")); barcode_file <- pick(c("barcodes.tsv"))
  if (!nzchar(matrix_file) || !nzchar(feature_file) || !nzchar(barcode_file)) {
    stop("A filtered 10x matrix directory must contain matrix.mtx(.gz), features.tsv/genes.tsv(.gz), and barcodes.tsv(.gz): ", path)
  }
  Seurat::ReadMtx(mtx = matrix_file, features = feature_file, cells = barcode_file, feature.column = 2)
}

assay_layers_safe <- function(obj, assay) {
  if (!assay %in% names(obj@assays)) return(character(0))
  tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
}

input_status_one <- function(obj, sample_id, kind) {
  rna_layers <- assay_layers_safe(obj, "RNA")
  annotation_columns <- c("cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype")
  annotation_columns <- annotation_columns[annotation_columns %in% colnames(obj@meta.data)]
  data.frame(
    sample_id = sample_id,
    input_kind = kind,
    cells_input = ncol(obj),
    features_input = nrow(obj),
    assays_detected = paste(names(obj@assays), collapse = "; "),
    rna_count_layers_detected = paste(grep("^counts($|\\.)", rna_layers, value = TRUE), collapse = "; "),
    rna_normalized_layers_detected = paste(grep("^data($|\\.)", rna_layers, value = TRUE), collapse = "; "),
    sct_detected = "SCT" %in% names(obj@assays),
    integrated_detected = "integrated" %in% names(obj@assays),
    reductions_detected = paste(names(obj@reductions), collapse = "; "),
    pca_detected = any(tolower(names(obj@reductions)) == "pca"),
    umap_detected = any(tolower(names(obj@reductions)) == "umap"),
    clusters_detected = any(c("seurat_clusters", "cluster") %in% colnames(obj@meta.data)),
    annotation_columns_detected = paste(annotation_columns, collapse = "; "),
    workflow_action = "Raw RNA counts retained; CodeSpring reruns selected QC and downstream processing reproducibly.",
    stringsAsFactors = FALSE
  )
}

save_input_umap <- function(obj, sample_id) {
  umap_name <- names(obj@reductions)[tolower(names(obj@reductions)) == "umap"]
  if (!length(umap_name)) return(invisible(NULL))
  preferred <- c("cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype", "condition", "Condition", "treatment", "Treatment", "group", "Group", "orig.ident", "batch", "sample_id")
  candidates <- unique(c(intersect(preferred, colnames(obj@meta.data)), colnames(obj@meta.data)))
  candidates <- candidates[vapply(candidates, function(col) {
    raw_values <- obj[[col]][, 1]
    if (is.numeric(raw_values) || is.integer(raw_values)) return(FALSE)
    values <- unique(as.character(raw_values))
    length(values[nzchar(values)]) > 1L && length(values[nzchar(values)]) <= 40L
  }, logical(1))]
  group_by <- if (length(candidates)) candidates[[1]] else NULL
  plot <- categorical_dimplot(obj, reduction = umap_name[[1]], group.by = group_by, shuffle = TRUE) + ggplot2::ggtitle("UMAP supplied with input object")
  file_name <- paste0("00_input_umap_", gsub("[^A-Za-z0-9_.-]", "_", sample_id), ".png")
  ggplot2::ggsave(file.path(figures_dir, file_name), plot = plot, width = 8, height = 6, dpi = 160)
  invisible(file_name)
}

read_one <- function(row) {
  path <- row[["input_path"]]; sample_id <- row[["sample_id"]]
  kind <- input_kind(path)
  obj <- switch(kind,
    seurat_rds = readRDS(path),
    scanpy_h5ad = stop("AnnData .h5ad input requires the Scanpy engine. CodeSpringApp chooses it automatically."),
    filtered_10x_matrix = {
      counts <- read_10x_matrix(path)
      if (is.list(counts)) counts <- counts[[if ("Gene Expression" %in% names(counts)) "Gene Expression" else 1]]
      Seurat::CreateSeuratObject(counts = counts, project = sample_id, min.cells = params$min_cells_per_gene)
    }
  )
  if (!inherits(obj, "Seurat")) stop("Input for ", sample_id, " did not yield a Seurat object.")
  # `cell` is reserved for the final processed barcode in CodeSpring output.
  # Preserve a same-named input metadata field before Seurat merge/QC steps,
  # which otherwise may silently drop it when objects are combined.
  if ("cell" %in% colnames(obj@meta.data)) names(obj@meta.data)[names(obj@meta.data) == "cell"] <- "input_cell"
  if (!"RNA" %in% names(obj@assays)) {
    assay <- DefaultAssay(obj)
    if (!nzchar(assay)) stop("No usable expression assay was found for ", sample_id)
    if (!length(grep("^counts($|\\.)", assay_layers_safe(obj, assay), value = TRUE))) {
      stop("Seurat input ", sample_id, " has no RNA assay with raw count layers. A full QC/normalization workflow must start from raw counts.")
    }
    obj[["RNA"]] <- obj[[assay]]
  }
  input_status <- input_status_one(obj, sample_id, kind)
  rna_count_layers <- grep("^counts($|\\.)", assay_layers_safe(obj, "RNA"), value = TRUE)
  if (!length(rna_count_layers)) stop("Seurat input ", sample_id, " has no raw RNA count layer. A full QC/normalization workflow must start from raw counts.")
  # Existing Seurat v5 objects may keep one count layer per sample. Join only
  # those raw layers before QC; normalized/integrated layers remain preserved
  # in the source object and are recorded in the provenance table.
  if (length(rna_count_layers) > 1L) obj <- SeuratObject::JoinLayers(obj, assay = "RNA", layers = "counts", new = "counts")
  DefaultAssay(obj) <- "RNA"
  # Preserve the unmodified barcode for cell-level annotation files. This is
  # safer than trying to strip a sample prefix later (sample IDs can contain
  # underscores themselves).
  obj$source_barcode <- colnames(obj)
  obj <- Seurat::RenameCells(obj, add.cell.id = sample_id)
  # Add all manifest fields as cell-level metadata.
  for (field in names(row)) obj[[field]] <- rep(as.character(row[[field]]), ncol(obj))
  obj[["input_kind"]] <- rep(kind, ncol(obj))
  obj@misc$codespring_input_status <- input_status
  if (identical(kind, "seurat_rds")) save_input_umap(obj, sample_id)
  obj
}

# Used both when QC figures are created and when a later annotation stage is
# resumed from the clustered checkpoint in a new SLURM process.
save_plot <- function(plot, file, width = 9, height = 6) {
  ggplot2::ggsave(file.path(figures_dir, file), plot = plot, width = width, height = height, dpi = 160)
}

if (stage %in% c("inspect", "qc", "all")) {
  objects <- lapply(seq_len(NROW(samples)), function(i) read_one(as.list(samples[i, , drop = FALSE])))
  names(objects) <- samples$sample_id
  cells_before_qc <- vapply(objects, ncol, numeric(1))
  input_processing_status <- do.call(rbind, lapply(objects, function(obj) obj@misc$codespring_input_status))
  utils::write.table(input_processing_status, file.path(tables_dir, "input_processing_detected.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  # Show raw distributions before a user commits to filtering thresholds.
  # These metrics are calculated from the retained raw RNA count layers and
  # become the evidence presented beside the QC controls in the app.
  objects <- lapply(objects, function(obj) {
    mt <- grep("^(MT-|mt-)", rownames(obj), value = TRUE)
    obj[["percent.mt"]] <- if (length(mt)) Seurat::PercentageFeatureSet(obj, features = mt) else rep(0, ncol(obj))
    obj
  })
  pre_qc_cells <- do.call(rbind, lapply(objects, function(obj) data.frame(
    sample_id = as.character(obj$sample_id), nCount_RNA = obj$nCount_RNA,
    nFeature_RNA = obj$nFeature_RNA, percent.mt = obj$percent.mt,
    stringsAsFactors = FALSE
  )))
  utils::write.table(pre_qc_cells, file.path(tables_dir, "qc_pre_filter_cell_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  pre_qc_merged <- Reduce(function(a, b) merge(a, y = b), objects)
  Seurat::DefaultAssay(pre_qc_merged) <- "RNA"
  save_plot(Seurat::VlnPlot(pre_qc_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", cols = categorical_palette(pre_qc_merged$sample_id), ncol = 3, pt.size = 0, layer = "counts"), "00_qc_pre_filter_violin.png", 14, 5)
  save_plot(jpplot_point_layer(Seurat::FeatureScatter(pre_qc_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")) + jpplot_point_layer(Seurat::FeatureScatter(pre_qc_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")), "00_qc_pre_filter_scatter.png", 12, 5)
  saveRDS(list(objects = objects, cells_before_qc = cells_before_qc, samples = samples), checkpoint_path("01_input"))
  stage_marker("inspect")
  if (identical(stage, "inspect")) quit(save = "no", status = 0L)
} else if (identical(stage, "qc")) {
  input_state <- require_checkpoint("01_input", "input inspection")
  objects <- input_state$objects
  cells_before_qc <- input_state$cells_before_qc
  samples <- input_state$samples
}

if (stage %in% c("qc", "all")) {
qc_one <- function(obj) {
  genes <- rownames(obj)
  mt <- grep("^(MT-|mt-)", genes, value = TRUE)
  obj[["percent.mt"]] <- if (length(mt)) Seurat::PercentageFeatureSet(obj, features = mt) else rep(0, ncol(obj))
  keep <- obj$nFeature_RNA >= params$min_features & obj$nCount_RNA >= params$min_counts & obj$percent.mt <= params$max_percent_mt
  if (params$max_features > 0) keep <- keep & obj$nFeature_RNA <= params$max_features
  if (!any(keep)) stop("QC thresholds removed every cell in sample ", unique(obj$sample_id)[1], ".")
  subset(obj, cells = colnames(obj)[keep])
}
objects <- lapply(objects, qc_one)

filter_genes_one <- function(obj) {
  counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  detected_cells <- Matrix::rowSums(counts > 0)
  keep <- detected_cells >= params$min_cells_per_gene
  if (!any(keep)) stop("Minimum cells per retained gene removed every gene in sample ", unique(obj$sample_id)[1], ".")
  summary <- data.frame(
    sample_id = unique(as.character(obj$sample_id))[[1]],
    genes_before_filtering = nrow(obj),
    min_cells_per_gene = params$min_cells_per_gene,
    genes_retained = sum(keep),
    stringsAsFactors = FALSE
  )
  list(object = subset(obj, features = rownames(obj)[keep]), summary = summary)
}

feature_filter_runs <- lapply(objects, filter_genes_one)
objects <- lapply(feature_filter_runs, `[[`, "object")
feature_filter_summary <- do.call(rbind, lapply(feature_filter_runs, `[[`, "summary"))
utils::write.table(feature_filter_summary, file.path(tables_dir, "feature_filtering_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sample_doublet_rate <- function(obj) {
  if (!"expected_doublet_rate" %in% colnames(obj@meta.data)) return(params$doublet_rate)
  raw <- unique(trimws(as.character(obj$expected_doublet_rate)))
  raw <- raw[nzchar(raw) & !is.na(raw)]
  if (!length(raw)) return(params$doublet_rate)
  if (length(raw) != 1L) stop("expected_doublet_rate must have one value per sample.")
  value <- suppressWarnings(as.numeric(raw[[1]]))
  if (!is.finite(value) || value <= 0 || value >= 1) stop("expected_doublet_rate must be greater than 0 and less than 1 for each sample.")
  value
}

run_doublet_one <- function(obj) {
  sample_id <- unique(as.character(obj$sample_id))[[1]]
  expected_rate <- if (identical(params$doublet_method, "none")) params$doublet_rate else sample_doublet_rate(obj)
  obj$doublet_score <- NA_real_
  obj$predicted_doublet <- FALSE
  if (identical(params$doublet_method, "none")) {
    return(list(
      object = obj,
      calls = data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = NA_real_, predicted_doublet = FALSE, removed_as_doublet = FALSE),
      summary = data.frame(sample_id = sample_id, method = "none", expected_doublet_rate = expected_rate, cells_before = ncol(obj), predicted_doublets = 0L, removed_doublets = 0L, note = "Doublet removal disabled")
    ))
  }
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    if (identical(params$doublet_method, "scdblfinder")) {
      stop("Doublet removal with the Seurat engine requires the Bioconductor package scDblFinder. Install it, choose the Scanpy/Scrublet engine, or explicitly select no doublet removal.")
    }
    return(list(
      object = obj,
      calls = data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = NA_real_, predicted_doublet = FALSE, removed_as_doublet = FALSE),
      summary = data.frame(sample_id = sample_id, method = "scDblFinder", expected_doublet_rate = expected_rate, cells_before = ncol(obj), predicted_doublets = 0L, removed_doublets = 0L, note = "Automatic scDblFinder skipped: package unavailable")
    ))
  }
  if (ncol(obj) < 100L || nrow(obj) < 20L) {
    return(list(
      object = obj,
      calls = data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = NA_real_, predicted_doublet = FALSE, removed_as_doublet = FALSE),
      summary = data.frame(sample_id = sample_id, method = "scDblFinder", expected_doublet_rate = expected_rate, cells_before = ncol(obj), predicted_doublets = 0L, removed_doublets = 0L, note = "Skipped: fewer than 100 cells or 20 genes")
    ))
  }
  sce <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder::scDblFinder(sce, dbr = expected_rate)
  score <- as.numeric(SummarizedExperiment::colData(sce)[["scDblFinder.score"]])
  call <- as.character(SummarizedExperiment::colData(sce)[["scDblFinder.class"]]) == "doublet"
  obj$doublet_score <- score
  obj$predicted_doublet <- call
  calls <- data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = score, predicted_doublet = call,
                      removed_as_doublet = call & params$remove_doublets, stringsAsFactors = FALSE)
  summary <- data.frame(sample_id = sample_id, method = "scDblFinder", expected_doublet_rate = expected_rate, cells_before = ncol(obj),
                        predicted_doublets = sum(call), removed_doublets = if (params$remove_doublets) sum(call) else 0L, note = "", stringsAsFactors = FALSE)
  if (params$remove_doublets) obj <- subset(obj, cells = colnames(obj)[!call])
  if (!ncol(obj)) stop("Doublet removal removed every cell in sample ", sample_id, ".")
  list(object = obj, calls = calls, summary = summary)
}

# Doublet detection is deliberately placed after initial per-sample QC and
# before normalization, feature selection, integration, and clustering.
doublet_runs <- lapply(objects, run_doublet_one)
objects <- lapply(doublet_runs, `[[`, "object")
doublet_calls <- do.call(rbind, lapply(doublet_runs, `[[`, "calls"))
doublet_summary <- do.call(rbind, lapply(doublet_runs, `[[`, "summary"))
utils::write.table(doublet_calls, file.path(tables_dir, "doublet_calls.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(doublet_summary, file.path(tables_dir, "doublet_summary_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

qc_cells <- do.call(rbind, lapply(objects, function(obj) {
  condition_value <- if ("condition" %in% colnames(obj@meta.data)) as.character(obj[["condition"]][, 1]) else rep("", ncol(obj))
  data.frame(cell = colnames(obj), sample_id = as.character(obj$sample_id), condition = condition_value,
             nCount_RNA = obj$nCount_RNA, nFeature_RNA = obj$nFeature_RNA, percent.mt = obj$percent.mt,
             doublet_score = obj$doublet_score, predicted_doublet = obj$predicted_doublet, stringsAsFactors = FALSE)
}))
utils::write.table(qc_cells, file.path(tables_dir, "qc_cell_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
qc_by_sample <- do.call(rbind, lapply(split(qc_cells, qc_cells$sample_id), function(x) data.frame(
  sample_id = x$sample_id[1], cells_input = unname(cells_before_qc[x$sample_id[1]]), cells_after_qc_and_doublets = NROW(x), median_umis = stats::median(x$nCount_RNA),
  median_genes = stats::median(x$nFeature_RNA), median_percent_mt = stats::median(x$percent.mt), stringsAsFactors = FALSE
)))
utils::write.table(qc_by_sample, file.path(tables_dir, "qc_summary_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

save_plot <- function(plot, file, width = 9, height = 6) ggplot2::ggsave(file.path(figures_dir, file), plot = plot, width = width, height = height, dpi = 160)
qc_merged <- Reduce(function(a, b) merge(a, y = b), objects)
DefaultAssay(qc_merged) <- "RNA"
save_plot(Seurat::VlnPlot(qc_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", cols = categorical_palette(qc_merged$sample_id), ncol = 3, pt.size = 0, layer = "counts"), "01_qc_violin.png", 14, 5)
save_plot(jpplot_point_layer(Seurat::FeatureScatter(qc_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")) + jpplot_point_layer(Seurat::FeatureScatter(qc_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")), "02_qc_scatter.png", 12, 5)
if (any(is.finite(doublet_calls$doublet_score))) {
  doublet_plot <- ggplot2::ggplot(doublet_calls[is.finite(doublet_calls$doublet_score), , drop = FALSE], ggplot2::aes(x = doublet_score)) +
    ggplot2::geom_histogram(bins = 40, fill = jpplot_colors[[7]], color = "white") +
    ggplot2::labs(title = "Doublet-score distribution before removal", x = "Doublet score", y = "Cells") +
    ggplot2::theme_classic()
  save_plot(doublet_plot, "02b_doublet_scores.png", 7, 4.5)
}

saveRDS(list(objects = objects, cells_before_qc = cells_before_qc, samples = samples, doublet_summary = doublet_summary), checkpoint_path("02_qc"))
stage_marker("qc")
if (identical(stage, "qc")) quit(save = "no", status = 0L)
} else if (identical(stage, "preprocess")) {
  qc_state <- require_checkpoint("02_qc", "QC and doublet handling")
  objects <- qc_state$objects
  cells_before_qc <- qc_state$cells_before_qc
  samples <- qc_state$samples
  doublet_summary <- qc_state$doublet_summary
}

integration <- params$integration
batch_values <- unlist(lapply(objects, function(obj) {
  if (params$batch_column %in% colnames(obj@meta.data)) as.character(obj[[params$batch_column]][, 1]) else character(0)
}), use.names = FALSE)
if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "rpca" else "none"
if (integration %in% c("rpca", "cca", "harmony") && length(unique(batch_values[nzchar(batch_values)])) < 2L) {
  stop(integration, " integration requires at least two values in the selected batch column (", params$batch_column, "). Choose none or supply the appropriate technical batch column.")
}
if (stage %in% c("preprocess", "all")) {
  # Anchor-based integration treats each object as an integration unit. Group
  # inputs by the explicitly selected technical batch first, so two biological
  # samples from the same library/run are not incorrectly corrected apart.
  if (integration %in% c("rpca", "cca")) {
    batch_for_object <- vapply(objects, function(obj) {
      values <- unique(trimws(as.character(obj[[params$batch_column]][, 1])))
      values <- values[nzchar(values) & !is.na(values)]
      if (length(values) != 1L) stop("Each input must have exactly one value for anchor-integration batch column ", params$batch_column, ".")
      values[[1]]
    }, character(1))
    grouped <- split(names(objects), batch_for_object)
    objects <- lapply(grouped, function(ids) {
      grouped_object <- Reduce(function(a, b) merge(a, y = b), objects[ids])
      # Seurat v5 merge retains one raw count layer per source sample. Within
      # one technical-batch integration unit these are not separate datasets,
      # so join them before normalization and anchor construction.
      if (length(grep("^counts($|\\.)", assay_layers_safe(grouped_object, "RNA"), value = TRUE)) > 1L) {
        grouped_object <- SeuratObject::JoinLayers(grouped_object, assay = "RNA", layers = "counts", new = "counts")
      }
      grouped_object
    })
  }
  if (identical(params$normalization, "sct")) {
    objects <- lapply(objects, function(obj) Seurat::SCTransform(obj, assay = "RNA", vst.flavor = "v2", verbose = FALSE))
  } else {
    objects <- lapply(objects, function(obj) {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
      Seurat::ScaleData(obj, verbose = FALSE)
    })
  }
  # Diagnostic embedding before any batch correction. Raw RNA counts remain
  # untouched and are retained for marker and differential-expression work.
  features <- Seurat::SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
  unintegrated <- Reduce(function(a, b) merge(a, y = b), objects)
  DefaultAssay(unintegrated) <- if (identical(params$normalization, "sct")) "SCT" else "RNA"
  Seurat::VariableFeatures(unintegrated) <- intersect(features, rownames(unintegrated))
  if (!identical(params$normalization, "sct")) unintegrated <- Seurat::ScaleData(unintegrated, features = Seurat::VariableFeatures(unintegrated), verbose = FALSE)
  pre_npcs <- max(2L, min(params$n_pcs, ncol(unintegrated) - 1L, length(Seurat::VariableFeatures(unintegrated)) - 1L))
  unintegrated <- Seurat::RunPCA(unintegrated, features = Seurat::VariableFeatures(unintegrated), npcs = pre_npcs, verbose = FALSE)
  pre_dims <- seq_len(min(params$n_pcs, ncol(Seurat::Embeddings(unintegrated, "pca"))))
  unintegrated <- Seurat::FindNeighbors(unintegrated, reduction = "pca", dims = pre_dims, verbose = FALSE)
  unintegrated <- Seurat::RunUMAP(unintegrated, reduction = "pca", dims = pre_dims, reduction.name = "umap.unintegrated", seed.use = params$seed, verbose = FALSE)
  save_plot(Seurat::DimPlot(unintegrated, reduction = "umap.unintegrated", group.by = "sample_id", shuffle = TRUE), "02_preintegration_umap_sample.png", 8, 6)
  if (nzchar(params$batch_column) && params$batch_column %in% colnames(unintegrated@meta.data) && length(unique(as.character(unintegrated[[params$batch_column]][, 1]))) > 1L && !identical(params$batch_column, "sample_id")) {
    save_plot(Seurat::DimPlot(unintegrated, reduction = "umap.unintegrated", group.by = params$batch_column, shuffle = TRUE), "02_preintegration_umap_batch.png", 8, 6)
  }
  pre_coords <- as.data.frame(Seurat::Embeddings(unintegrated, "umap.unintegrated"))
  names(pre_coords)[1:2] <- c("UMAP_1", "UMAP_2"); pre_coords$cell <- rownames(pre_coords); pre_coords$sample_id <- as.character(unintegrated$sample_id)
  if (nzchar(params$batch_column) && params$batch_column %in% colnames(unintegrated@meta.data)) pre_coords[[params$batch_column]] <- as.character(unintegrated[[params$batch_column]][, 1])
  utils::write.table(pre_coords[, c("cell", setdiff(names(pre_coords), "cell")), drop = FALSE], file.path(tables_dir, "preintegration_umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  rm(unintegrated); invisible(gc())
  saveRDS(list(objects = objects, cells_before_qc = cells_before_qc, samples = samples, doublet_summary = doublet_summary), checkpoint_path("03_preprocessed"))
  stage_marker("preprocess")
  if (identical(stage, "preprocess")) quit(save = "no", status = 0L)
} else if (identical(stage, "cluster")) {
  pre_state <- require_checkpoint("03_preprocessed", "normalization and PCA")
  objects <- pre_state$objects
  cells_before_qc <- pre_state$cells_before_qc
  samples <- pre_state$samples
  doublet_summary <- pre_state$doublet_summary
}

# Determine integration only after loading the preprocessed checkpoint, so a
# resumed stage uses the exact same retained cells and technical metadata.
if (stage %in% c("cluster", "all")) {
integration <- params$integration
batch_values <- unlist(lapply(objects, function(obj) {
  if (params$batch_column %in% colnames(obj@meta.data)) as.character(obj[[params$batch_column]][, 1]) else character(0)
}), use.names = FALSE)
if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "rpca" else "none"
if (integration %in% c("rpca", "cca", "harmony") && length(unique(batch_values[nzchar(batch_values)])) < 2L) {
  stop(integration, " integration requires at least two values in the selected batch column (", params$batch_column, "). Choose none or supply the appropriate technical batch column.")
}

  integration_k_weight <- max(5L, min(50L, floor(min(vapply(objects, ncol, numeric(1))) / 2)))
  reduction_for_graph <- "pca"
  if (identical(integration, "harmony")) {
    if (!requireNamespace("harmony", quietly = TRUE)) stop("Seurat Harmony integration requires the R package harmony.")
    features <- Seurat::SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
    obj <- Reduce(function(a, b) merge(a, y = b), objects)
    DefaultAssay(obj) <- if (identical(params$normalization, "sct")) "SCT" else "RNA"
    Seurat::VariableFeatures(obj) <- intersect(features, rownames(obj))
    if (!identical(params$normalization, "sct")) obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
    obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = params$n_pcs, verbose = FALSE)
    harmony_dims <- seq_len(min(params$n_pcs, ncol(Seurat::Embeddings(obj, "pca"))))
    obj <- harmony::RunHarmony(obj, group.by.vars = params$batch_column, reduction.use = "pca", dims.use = harmony_dims, theta = params$harmony_theta, lambda = params$harmony_lambda, max_iter = params$harmony_max_iter, reduction.save = "harmony", verbose = FALSE)
    reduction_for_graph <- "harmony"
  } else if (identical(params$normalization, "sct")) {
    if (integration %in% c("rpca", "cca") && length(objects) > 1L) {
      features <- Seurat::SelectIntegrationFeatures(objects, nfeatures = 3000)
      objects <- Seurat::PrepSCTIntegration(objects, anchor.features = features, verbose = FALSE)
      if (identical(integration, "rpca")) objects <- lapply(objects, function(obj) Seurat::RunPCA(obj, features = features, npcs = params$n_pcs, verbose = FALSE))
      anchors <- Seurat::FindIntegrationAnchors(object.list = objects, normalization.method = "SCT", anchor.features = features,
                                                 reduction = integration, dims = seq_len(params$n_pcs), verbose = FALSE)
      obj <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = seq_len(params$n_pcs), k.weight = integration_k_weight, verbose = FALSE)
      DefaultAssay(obj) <- "integrated"; obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
    } else {
      obj <- Reduce(function(a, b) merge(a, y = b), objects)
      DefaultAssay(obj) <- "SCT"; obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
    }
  } else if (integration %in% c("rpca", "cca") && length(objects) > 1L) {
    features <- Seurat::SelectIntegrationFeatures(objects, nfeatures = 3000)
    if (identical(integration, "rpca")) objects <- lapply(objects, function(obj) Seurat::RunPCA(obj, features = features, npcs = params$n_pcs, verbose = FALSE))
    anchors <- Seurat::FindIntegrationAnchors(object.list = objects, anchor.features = features, reduction = integration,
                                               dims = seq_len(params$n_pcs), verbose = FALSE)
    obj <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(params$n_pcs), k.weight = integration_k_weight, verbose = FALSE)
    DefaultAssay(obj) <- "integrated"; obj <- Seurat::ScaleData(obj, verbose = FALSE); obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  } else {
    obj <- Reduce(function(a, b) merge(a, y = b), objects)
    DefaultAssay(obj) <- "RNA"; obj <- Seurat::NormalizeData(obj, verbose = FALSE); obj <- Seurat::FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE); obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  }
  n_available <- ncol(Seurat::Embeddings(obj, reduction_for_graph)); dims <- seq_len(min(params$n_pcs, n_available))
  pca_stdev <- Seurat::Stdev(obj, reduction = "pca"); pca_variance <- (pca_stdev ^ 2) / sum(pca_stdev ^ 2)
  utils::write.table(data.frame(PC = seq_along(pca_variance), variance_explained = pca_variance, percent_variance_explained = 100 * pca_variance), file.path(tables_dir, "pca_variance_explained.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  variable_genes <- Seurat::VariableFeatures(obj)
  utils::write.table(data.frame(gene = variable_genes, highly_variable = TRUE), file.path(tables_dir, "highly_variable_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_for_graph, dims = dims, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = params$cluster_resolution, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = reduction_for_graph, dims = dims, seed.use = params$seed, verbose = FALSE)
  obj$cluster <- as.character(Seurat::Idents(obj))
  saveRDS(list(object = obj, integration = integration, samples = samples, doublet_summary = doublet_summary), checkpoint_path("04_clustered"))
  stage_marker("cluster")
  if (identical(stage, "cluster")) quit(save = "no", status = 0L)
} else if (identical(stage, "annotate")) {
  processed_path <- file.path(objects_dir, "processed_seurat.rds")
  clustered_path <- checkpoint_path("04_clustered")
  use_processed <- file.exists(processed_path) && (!file.exists(clustered_path) || file.info(processed_path)$mtime >= file.info(clustered_path)$mtime)
  if (use_processed) {
    obj <- readRDS(processed_path)
    integration <- attr(obj, "codespring_integration") %||% params$integration
    doublet_summary <- data.frame(removed_doublets = attr(obj, "codespring_doublets_removed") %||% 0)
  } else {
    cluster_state <- require_checkpoint("04_clustered", "integration and clustering")
    obj <- cluster_state$object
    integration <- cluster_state$integration
    samples <- cluster_state$samples
    doublet_summary <- cluster_state$doublet_summary
  }
} else if (stage %in% c("score", "differential")) {
  processed_path <- file.path(objects_dir, "processed_seurat.rds")
  if (!file.exists(processed_path)) stop("The annotation stage has not completed. Run annotation before ", stage, ".")
  obj <- readRDS(processed_path)
  integration <- attr(obj, "codespring_integration") %||% params$integration
  doublet_summary <- data.frame(removed_doublets = 0)
}

safe_metadata_name <- function(value, default = "cell_type") {
  value <- gsub("[^A-Za-z0-9_]+", "_", trimws(as.character(value)))
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) value <- default
  if (grepl("^[0-9]", value)) value <- paste0("annotation_", value)
  value
}

apply_celltype_mapping <- function(obj, path, annotation_name) {
  if (!nzchar(path) || identical(tolower(path), "none")) return(obj)
  map <- read_delim_safe(path)
  cell_col <- intersect(c("cell", "cell_id", "barcode", "Cell", "CellID"), names(map))
  type_col <- intersect(c("cell_type", "celltype", "CellType", "annotation"), names(map))
  if (!length(cell_col) || !length(type_col)) stop("Cell-type mapping needs a cell/barcode column and a cell_type column.")
  labels <- setNames(as.character(map[[type_col[[1]]]]), as.character(map[[cell_col[[1]]]]))
  direct <- unname(labels[colnames(obj)])
  raw_barcode <- as.character(obj$source_barcode)
  direct[is.na(direct)] <- unname(labels[raw_barcode[is.na(direct)]])
  obj[[annotation_name]] <- ifelse(is.na(direct) | !nzchar(direct), "Unassigned", direct)
  obj[[paste0("annotation_source__", annotation_name)]] <- "provided cell metadata"
  obj$annotation_source <- paste0("provided cell metadata (", annotation_name, ")")
  obj
}

apply_marker_annotation <- function(obj, path, annotation_name) {
  if (!nzchar(path) || identical(tolower(path), "none")) return(obj)
  markers <- read_delim_safe(path)
  if (!all(c("cell_type", "gene") %in% names(markers))) stop("Marker list needs cell_type and gene columns.")
  marker_list <- split(unlist(strsplit(as.character(markers$gene), "[,;[:space:]]+")), as.character(markers$cell_type))
  marker_list <- lapply(marker_list, function(x) intersect(unique(x[nzchar(x)]), rownames(obj)))
  marker_list <- marker_list[lengths(marker_list) > 0]
  if (!length(marker_list)) stop("None of the supplied marker genes were present in the expression object.")
  # Average normalized expression per cluster is stable and transparent for
  # cluster-level marker-list annotation.
  DefaultAssay(obj) <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  scores <- sapply(marker_list, function(genes) Matrix::colMeans(Seurat::GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")[genes, , drop = FALSE]))
  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1, dimnames = list(colnames(obj), names(marker_list)))
  cluster <- as.character(obj$cluster)
  means <- rowsum(scores, group = cluster) / as.numeric(table(cluster)[rownames(rowsum(scores, group = cluster))])
  assignment <- colnames(means)[max.col(means, ties.method = "first")]
  cluster_label <- setNames(assignment, rownames(means))
  obj[[annotation_name]] <- unname(cluster_label[cluster])
  obj[[paste0("annotation_source__", annotation_name)]] <- "marker-list cluster scoring"
  obj$annotation_source <- paste0("marker-list cluster scoring (", annotation_name, ")")
  score_table <- data.frame(cluster = rownames(means), means, check.names = FALSE)
  utils::write.table(score_table, file.path(tables_dir, "marker_annotation_cluster_scores.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(score_table, file.path(tables_dir, paste0("marker_annotation_cluster_scores__", annotation_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  obj
}

apply_existing_annotation <- function(obj) {
  candidates <- c("cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype")
  candidates <- candidates[candidates %in% colnames(obj@meta.data)]
  if (!length(candidates)) return(NULL)
  for (column in candidates) {
    values <- as.character(obj[[column]][, 1])
    if (any(nzchar(values) & !is.na(values))) {
      obj$cell_type <- ifelse(is.na(values) | !nzchar(values), "Unassigned", values)
      obj$annotation_source <- paste0("existing input metadata: ", column)
      return(obj)
    }
  }
  NULL
}

if (stage %in% c("annotate", "all")) {
annotation_name <- safe_metadata_name(params$annotation_name)
# Rejoin Seurat v5 RNA layers before any annotation scoring or differential
# expression. Values are unchanged; the final object becomes portable and its
# raw counts remain directly accessible.
rna_layers <- assay_layers_safe(obj, "RNA")
if (any(grepl("^(counts|data)\\.", rna_layers))) obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
if (nzchar(params$celltype_file) && !identical(tolower(params$celltype_file), "none")) {
  obj <- apply_celltype_mapping(obj, params$celltype_file, annotation_name)
} else if (nzchar(params$marker_file) && !identical(tolower(params$marker_file), "none")) {
  obj <- apply_marker_annotation(obj, params$marker_file, annotation_name)
} else {
  existing_annotation <- apply_existing_annotation(obj)
  if (is.null(existing_annotation)) {
    obj[[annotation_name]] <- obj$cluster
    obj[[paste0("annotation_source__", annotation_name)]] <- "cluster ID (no annotation supplied)"
    obj$annotation_source <- paste0("cluster ID (", annotation_name, "; no annotation supplied)")
  } else {
    obj <- existing_annotation
    if (!identical(annotation_name, "cell_type")) {
      obj[[annotation_name]] <- obj$cell_type
      obj[[paste0("annotation_source__", annotation_name)]] <- obj$annotation_source
    }
  }
}

save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "sample_id", shuffle = TRUE), "03_umap_sample.png", 8, 6)
save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "cluster", label = TRUE, repel = TRUE), "04_umap_clusters.png", 8, 6)
save_plot(categorical_dimplot(obj, reduction = "umap", group.by = annotation_name, label = TRUE, repel = TRUE), paste0("05_umap_", annotation_name, ".png"), 9, 6)
if ("condition" %in% colnames(obj@meta.data) && length(unique(obj$condition)) > 1L) save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "condition", shuffle = TRUE), "06_umap_condition.png", 8, 6)

# Cluster markers should use the normalized SCT representation when that is
# the selected workflow.
DefaultAssay(obj) <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
if (ncol(obj) >= 20 && length(unique(obj$cluster)) > 1L) {
  # An SCT-integrated object can retain one model per input sample. Bring
  # those models to a common sequencing-depth scale before marker testing;
  # otherwise Seurat may silently return an empty marker table.
  if (identical(DefaultAssay(obj), "SCT") && length(obj[["SCT"]]@SCTModel.list) > 1L) {
    obj <- Seurat::PrepSCTFindMarkers(obj, verbose = FALSE)
  }
  markers <- tryCatch(Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE), error = function(e) data.frame(error = conditionMessage(e)))
  utils::write.table(markers, file.path(tables_dir, "cluster_markers.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!"error" %in% names(markers) && NROW(markers)) {
    gene_col <- intersect(c("gene", "features"), names(markers))
    if (length(gene_col)) {
      top <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) utils::head(x[order(x$p_val_adj), , drop = FALSE], 10)))
      utils::write.table(top, file.path(tables_dir, "top10_markers_per_cluster.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

# Match the Scanpy output schema: retain an input metadata column named
# `cell` as `input_cell`, while always reserving `cell` for the exact final
# barcode/cell identifier.  This prevents duplicate column names in tables
# consumed by the interactive results explorer.
cell_metadata <- obj@meta.data
if ("cell" %in% names(cell_metadata)) names(cell_metadata)[names(cell_metadata) == "cell"] <- "input_cell"
cell_metadata <- data.frame(cell = colnames(obj), cell_metadata, check.names = FALSE)
utils::write.table(cell_metadata, file.path(tables_dir, "cell_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# The interactive dashboard requests one normalized gene at a time from the
# processed RDS. Keep the complete symbol list small and engine-consistent.
utils::write.table(data.frame(gene = rownames(obj)), file.path(tables_dir, "dashboard_all_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
umap_coordinates <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"), check.names = FALSE)
colnames(umap_coordinates)[seq_len(min(2L, NCOL(umap_coordinates)))] <- c("UMAP_1", "UMAP_2")[seq_len(min(2L, NCOL(umap_coordinates)))]
if (!all(c("UMAP_1", "UMAP_2") %in% names(umap_coordinates))) stop("The final UMAP does not contain two coordinates.")
umap_metadata <- obj@meta.data
if ("cell" %in% names(umap_metadata)) names(umap_metadata)[names(umap_metadata) == "cell"] <- "input_cell"
umap_table <- data.frame(cell = rownames(umap_coordinates), umap_coordinates[, c("UMAP_1", "UMAP_2"), drop = FALSE], umap_metadata, check.names = FALSE)
utils::write.table(umap_table, file.path(tables_dir, "umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cluster_sizes <- as.data.frame(table(cluster = obj$cluster, annotation_label = obj[[annotation_name]][, 1]), stringsAsFactors = FALSE)
names(cluster_sizes)[names(cluster_sizes) == "Freq"] <- "cells"
cluster_sizes <- cluster_sizes[cluster_sizes$cells > 0, , drop = FALSE]
cluster_sizes$annotation_field <- annotation_name
utils::write.table(cluster_sizes[, c("cluster", "annotation_field", "annotation_label", "cells")], file.path(tables_dir, "cluster_annotation_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
legacy_cluster <- cluster_sizes[, c("cluster", "annotation_label", "cells")]; names(legacy_cluster)[2] <- "cell_type"
utils::write.table(legacy_cluster, file.path(tables_dir, "cluster_cell_type_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cell_type_by_sample <- as.data.frame(table(sample_id = obj$sample_id, annotation_label = obj[[annotation_name]][, 1]), stringsAsFactors = FALSE)
names(cell_type_by_sample)[names(cell_type_by_sample) == "Freq"] <- "cells"
cell_type_by_sample <- cell_type_by_sample[cell_type_by_sample$cells > 0, , drop = FALSE]
cell_type_by_sample$proportion_within_sample <- cell_type_by_sample$cells / ave(cell_type_by_sample$cells, cell_type_by_sample$sample_id, FUN = sum)
cell_type_by_sample$annotation_field <- annotation_name
composition_path <- file.path(tables_dir, "composition_by_sample.tsv")
if (file.exists(composition_path)) {
  previous <- read_delim_safe(composition_path)
  if ("annotation_field" %in% names(previous)) previous <- previous[as.character(previous$annotation_field) != annotation_name, , drop = FALSE]
  cell_type_by_sample <- rbind(previous[, names(cell_type_by_sample), drop = FALSE], cell_type_by_sample)
}
utils::write.table(cell_type_by_sample[, c("sample_id", "annotation_field", "annotation_label", "cells", "proportion_within_sample")], composition_path, sep = "\t", row.names = FALSE, quote = FALSE)
legacy_composition <- cell_type_by_sample[as.character(cell_type_by_sample$annotation_field) == annotation_name, c("sample_id", "annotation_label", "cells", "proportion_within_sample")]; names(legacy_composition)[2] <- "cell_type"
utils::write.table(legacy_composition, file.path(tables_dir, "cell_type_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

attr(obj, "codespring_active_annotation") <- annotation_name
attr(obj, "codespring_integration") <- integration
attr(obj, "codespring_doublets_removed") <- sum(doublet_summary$removed_doublets)
saveRDS(obj, file.path(objects_dir, "processed_seurat.rds"))
summary_lines <- c(
  paste("engine: seurat"), paste("normalization:", params$normalization), paste("integration:", integration),
  paste("doublet_method:", params$doublet_method), paste("doublets_removed:", sum(doublet_summary$removed_doublets)),
  paste("input_processing_inventory:", file.path("tables", "input_processing_detected.tsv")),
  paste("input_samples:", NROW(samples)), paste("cells_after_qc:", ncol(obj)), paste("clusters:", length(unique(obj$cluster))), paste("active_annotation:", annotation_name),
  paste("annotation_source:", unique(obj$annotation_source)[1]), paste("generated:", as.character(Sys.time()))
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
stage_marker("annotate")
writeLines(as.character(Sys.time()), file.path(out_dir, "_COMPLETE"))
}

if (identical(stage, "score")) {
  if (!nzchar(params$signature_file) || identical(tolower(params$signature_file), "none")) stop("Choose a signature TSV with signature and gene columns.")
  signatures <- read_delim_safe(params$signature_file)
  if (!all(c("signature", "gene") %in% names(signatures))) stop("Signature list needs signature and gene columns.")
  annotation_name <- safe_metadata_name(attr(obj, "codespring_active_annotation") %||% params$annotation_name)
  DefaultAssay(obj) <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  signature_names <- unique(as.character(signatures$signature))
  coverage <- list(); score_columns <- character(0)
  for (signature in signature_names) {
    requested <- unique(unlist(strsplit(as.character(signatures$gene[as.character(signatures$signature) == signature]), "[,;[:space:]]+")))
    requested <- requested[nzchar(requested)]
    present <- intersect(requested, rownames(obj))
    column <- paste0("signature__", safe_metadata_name(signature, "score"))
    coverage[[length(coverage) + 1L]] <- data.frame(signature = signature, metadata_column = if (length(present)) column else "", genes_requested = length(requested), genes_found = length(present), coverage = length(present) / max(1L, length(requested)))
    if (!length(present)) next
    nbin <- min(24L, max(2L, floor(nrow(obj) / 5L)))
    ctrl <- min(100L, max(1L, floor((nrow(obj) - length(present)) / nbin)))
    obj <- Seurat::AddModuleScore(obj, features = list(present), name = paste0(column, "__tmp"), assay = DefaultAssay(obj), nbin = nbin, ctrl = ctrl, seed = params$seed)
    generated <- grep(paste0("^", column, "__tmp"), colnames(obj@meta.data), value = TRUE)[1]
    obj[[column]] <- obj[[generated]][, 1]
    obj[[generated]] <- NULL
    score_columns <- c(score_columns, column)
  }
  utils::write.table(do.call(rbind, coverage), file.path(tables_dir, "signature_gene_coverage.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!length(score_columns)) stop("None of the supplied signature genes were found in the normalized expression object.")
  score_metadata <- obj@meta.data[, unique(c("sample_id", annotation_name, score_columns)), drop = FALSE]
  score_metadata <- data.frame(cell = rownames(score_metadata), score_metadata, check.names = FALSE)
  utils::write.table(score_metadata, file.path(tables_dir, "signature_scores_per_cell.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  grouping <- interaction(obj@meta.data[, intersect(c("sample_id", annotation_name), colnames(obj@meta.data)), drop = FALSE], drop = TRUE, sep = " | ")
  summary_rows <- do.call(rbind, lapply(split(seq_len(ncol(obj)), grouping), function(idx) {
    base <- obj@meta.data[idx[1], intersect(c("sample_id", annotation_name), colnames(obj@meta.data)), drop = FALSE]
    stats <- unlist(lapply(score_columns, function(column) c(mean = mean(obj@meta.data[idx, column]), median = stats::median(obj@meta.data[idx, column]), count = length(idx))))
    names(stats) <- unlist(lapply(score_columns, function(column) paste0(column, "__", c("mean", "median", "count"))))
    data.frame(base, as.list(stats), check.names = FALSE)
  }))
  utils::write.table(summary_rows, file.path(tables_dir, "signature_scores_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cell_metadata <- obj@meta.data; if ("cell" %in% names(cell_metadata)) names(cell_metadata)[names(cell_metadata) == "cell"] <- "input_cell"
  utils::write.table(data.frame(cell = colnames(obj), cell_metadata, check.names = FALSE), file.path(tables_dir, "cell_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  umap_coordinates <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap")); names(umap_coordinates)[1:2] <- c("UMAP_1", "UMAP_2")
  utils::write.table(data.frame(cell = rownames(umap_coordinates), umap_coordinates[, 1:2], cell_metadata, check.names = FALSE), file.path(tables_dir, "umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(obj, file.path(objects_dir, "processed_seurat.rds"))
  stage_marker("score")
}

if (identical(stage, "differential")) {
  group_column <- params$de_group_column
  reference <- params$de_reference; comparison <- params$de_comparison
  if (!group_column %in% colnames(obj@meta.data)) stop("Differential-expression group field is absent from cell metadata: ", group_column)
  if (!nzchar(reference) || !nzchar(comparison) || identical(reference, comparison)) stop("Choose two different reference and comparison groups for differential expression.")
  targets <- unique(params$de_annotation_values %||% "all")
  if (any(!tolower(targets) %in% c("all", "all cells")) && (!nzchar(params$de_annotation_column) || !params$de_annotation_column %in% colnames(obj@meta.data))) stop("Choose an annotation metadata field before requesting individual-population differential expression.")
  contrast_slug <- paste0(safe_metadata_name(comparison, "comparison"), "_vs_", safe_metadata_name(reference, "reference"))
  root <- file.path(out_dir, "differential_expression", contrast_slug); dir.create(root, recursive = TRUE, showWarnings = FALSE)
  manifest_rows <- list(); completed <- 0L
  for (annotation_value in targets) {
    is_global <- tolower(annotation_value) %in% c("all", "all cells")
    population_slug <- if (is_global) "global" else paste0(safe_metadata_name(params$de_annotation_column, "annotation"), "__", safe_metadata_name(annotation_value, "population"))
    job_dir <- file.path(root, population_slug); dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    file_slug <- paste0(contrast_slug, "__", population_slug)
    row <- data.frame(job_id = file_slug, comparison = comparison, reference = reference, annotation_column = if (is_global) "" else params$de_annotation_column, annotation_value = if (is_global) "all" else annotation_value, output_dir = normalizePath(job_dir, winslash = "/", mustWork = FALSE), file_slug = file_slug, status = "prepared", note = "", count_file = "", design_file = "", cell_result_file = "", stringsAsFactors = FALSE)
    error <- tryCatch({
      keep <- as.character(obj@meta.data[[group_column]]) %in% c(reference, comparison)
      if (!is_global) keep <- keep & as.character(obj@meta.data[[params$de_annotation_column]]) == annotation_value
      subset_obj <- subset(obj, cells = colnames(obj)[keep])
      if (!ncol(subset_obj)) stop("No cells remain for this population.")
      if (!all(c(reference, comparison) %in% unique(as.character(subset_obj@meta.data[[group_column]])))) stop("Both comparison groups must contain cells in this population.")
      if (params$de_method %in% c("both", "pseudobulk")) {
        count_layers <- grep("^counts($|\\.)", assay_layers_safe(subset_obj, "RNA"), value = TRUE)
        if (length(count_layers) > 1L) subset_obj <- SeuratObject::JoinLayers(subset_obj, assay = "RNA", layers = "counts", new = "counts")
        counts <- Seurat::GetAssayData(subset_obj, assay = "RNA", layer = "counts")
        count_columns <- list(); design_rows <- list()
        for (sample_id in unique(as.character(subset_obj$sample_id))) {
          cells <- which(as.character(subset_obj$sample_id) == sample_id)
          groups <- unique(as.character(subset_obj@meta.data[[group_column]][cells]))
          if (length(groups) != 1L) stop("Biological replicate ", sample_id, " contains multiple comparison-group values.")
          count_columns[[sample_id]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
          design_row <- data.frame(sample_id = sample_id, group = groups[[1]], cells = length(cells), check.names = FALSE)
          for (covariate in params$de_covariates) {
            if (!covariate %in% colnames(subset_obj@meta.data)) stop("Requested covariate is absent: ", covariate)
            values <- unique(as.character(subset_obj@meta.data[[covariate]][cells]))
            if (length(values) != 1L) stop("Covariate ", covariate, " is not constant within replicate ", sample_id, ".")
            design_row[[covariate]] <- values[[1]]
          }
          design_rows[[sample_id]] <- design_row
        }
        pseudobulk <- do.call(cbind, count_columns)
        count_file <- file.path(job_dir, paste0("pseudobulk_counts__", file_slug, ".tsv")); design_file <- file.path(job_dir, paste0("pseudobulk_design__", file_slug, ".tsv"))
        utils::write.table(data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE), count_file, sep = "\t", row.names = FALSE, quote = FALSE)
        utils::write.table(do.call(rbind, design_rows), design_file, sep = "\t", row.names = FALSE, quote = FALSE)
        row$count_file <- count_file; row$design_file <- design_file
      }
      if (params$de_method %in% c("both", "cell", "cell_level")) {
        cell_error <- tryCatch({
          DefaultAssay(subset_obj) <- if ("SCT" %in% names(subset_obj@assays)) "SCT" else "RNA"
          if (identical(DefaultAssay(subset_obj), "SCT") && length(subset_obj[["SCT"]]@SCTModel.list) > 1L) subset_obj <- Seurat::PrepSCTFindMarkers(subset_obj, verbose = FALSE)
          Seurat::Idents(subset_obj) <- as.character(subset_obj@meta.data[[group_column]])
          cell_de <- Seurat::FindMarkers(subset_obj, ident.1 = comparison, ident.2 = reference, test.use = "wilcox", logfc.threshold = 0, min.pct = 0, verbose = FALSE)
          cell_de <- data.frame(gene = rownames(cell_de), comparison = paste0(comparison, "_vs_", reference), population = if (is_global) "global" else annotation_value, cell_de, analysis_level = "cell-level exploratory Wilcoxon; cells are not biological replicates", check.names = FALSE)
          cell_file <- file.path(job_dir, paste0("cell_level_Wilcoxon__", file_slug, ".tsv")); utils::write.table(cell_de, cell_file, sep = "\t", row.names = FALSE, quote = FALSE); row$cell_result_file <- cell_file
          if (is_global) utils::write.table(cell_de, file.path(tables_dir, "cell_level_differential_expression.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
          ""
        }, error = function(e) conditionMessage(e))
        if (nzchar(cell_error)) row$note <- paste0("Cell-level test skipped: ", cell_error)
      }
      completed <- completed + 1L
      ""
    }, error = function(e) conditionMessage(e))
    if (nzchar(error)) { row$status <- "skipped"; row$note <- error }
    manifest_rows[[length(manifest_rows) + 1L]] <- row
  }
  utils::write.table(do.call(rbind, manifest_rows), file.path(tables_dir, "differential_jobs.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!completed) stop("No requested differential-expression population could be prepared. Review tables/differential_jobs.tsv.")
}
