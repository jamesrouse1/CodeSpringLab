#!/usr/bin/env Rscript

# CodeSpringLab single-cell RNA-seq workflow (Seurat engine)
#
# Usage:
#   Rscript scrna_pipeline_seurat.R <samples.tsv> <out_dir> <params.tsv>
#
# samples.tsv (required columns): sample_id, input_path
# Optional sample columns (for example condition and batch) are retained as
# per-cell metadata.  input_path may be a filtered 10x matrix directory or a
# .rds Seurat object. AnnData .h5ad inputs are handled natively by Scanpy.
# params.tsv contains key/value pairs; comments and unknown keys are allowed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: scrna_pipeline_seurat.R <samples.tsv> <out_dir> <params.tsv>")

samples_path <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
params_path <- normalizePath(args[[3]], mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is not installed: ", pkg)
}
for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork")) need_pkg(pkg)
suppressPackageStartupMessages(library(Seurat))

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x

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
    seed = suppressWarnings(as.integer(get("seed", "1234")))
  )
}

params <- read_params(params_path)
if (!params$normalization %in% c("sct", "lognormalize")) stop("normalization must be SCT or LogNormalize")
if (!params$integration %in% c("auto", "none", "rpca", "cca")) stop("integration must be auto, none, rpca, or cca")
if (!is.finite(params$cluster_resolution) || params$cluster_resolution <= 0) params$cluster_resolution <- 0.6
if (!is.finite(params$n_pcs) || params$n_pcs < 5) params$n_pcs <- 30L
if (!is.finite(params$min_features) || params$min_features < 0) params$min_features <- 200L
if (!is.finite(params$min_counts) || params$min_counts < 0) params$min_counts <- 0L
if (!is.finite(params$max_features) || params$max_features < 0) params$max_features <- 0L
if (!is.finite(params$max_percent_mt) || params$max_percent_mt < 0) params$max_percent_mt <- 20
if (!is.finite(params$min_cells_per_gene) || params$min_cells_per_gene < 1) params$min_cells_per_gene <- 3L
if (!params$doublet_method %in% c("auto", "none", "scdblfinder")) stop("doublet_method must be auto, none, or scDblFinder")
if (!is.finite(params$doublet_rate) || params$doublet_rate < 0 || params$doublet_rate >= 1) params$doublet_rate <- 0.05
if (!is.finite(params$seed)) params$seed <- 1234L
set.seed(params$seed)

samples <- read_delim_safe(samples_path)
required <- c("sample_id", "input_path")
if (!all(required %in% names(samples))) stop("Sample manifest requires columns: ", paste(required, collapse = ", "))
samples <- samples[nzchar(trimws(as.character(samples$sample_id))) & nzchar(trimws(as.character(samples$input_path))), , drop = FALSE]
if (!NROW(samples)) stop("No valid sample rows were supplied.")
samples$sample_id <- make.unique(make.names(as.character(samples$sample_id)))
samples$input_path <- normalizePath(path.expand(as.character(samples$input_path)), winslash = "/", mustWork = FALSE)
missing <- samples$input_path[!file.exists(samples$input_path)]
if (length(missing)) stop("Missing unreadable input path(s): ", paste(missing, collapse = ", "))

tables_dir <- file.path(out_dir, "tables"); figures_dir <- file.path(out_dir, "figures"); objects_dir <- file.path(out_dir, "objects")
for (d in c(tables_dir, figures_dir, objects_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
utils::write.table(samples, file.path(tables_dir, "sample_manifest_used.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

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
  if (!"RNA" %in% names(obj@assays)) {
    assay <- DefaultAssay(obj)
    if (!nzchar(assay)) stop("No usable expression assay was found for ", sample_id)
    obj[["RNA"]] <- obj[[assay]]
  }
  DefaultAssay(obj) <- "RNA"
  obj <- Seurat::RenameCells(obj, add.cell.id = sample_id)
  # Add all manifest fields as cell-level metadata.
  for (field in names(row)) obj[[field]] <- rep(as.character(row[[field]]), ncol(obj))
  obj[["input_kind"]] <- rep(kind, ncol(obj))
  obj
}

objects <- lapply(seq_len(NROW(samples)), function(i) read_one(as.list(samples[i, , drop = FALSE])))
names(objects) <- samples$sample_id
cells_before_qc <- vapply(objects, ncol, numeric(1))

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

run_doublet_one <- function(obj) {
  sample_id <- unique(as.character(obj$sample_id))[[1]]
  obj$doublet_score <- NA_real_
  obj$predicted_doublet <- FALSE
  if (identical(params$doublet_method, "none")) {
    return(list(
      object = obj,
      calls = data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = NA_real_, predicted_doublet = FALSE, removed_as_doublet = FALSE),
      summary = data.frame(sample_id = sample_id, method = "none", cells_before = ncol(obj), predicted_doublets = 0L, removed_doublets = 0L, note = "Doublet removal disabled")
    ))
  }
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("Doublet removal with the Seurat engine requires the Bioconductor package scDblFinder. Install it, choose the Scanpy/Scrublet engine, or explicitly select no doublet removal.")
  }
  if (ncol(obj) < 100L || nrow(obj) < 20L) {
    return(list(
      object = obj,
      calls = data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = NA_real_, predicted_doublet = FALSE, removed_as_doublet = FALSE),
      summary = data.frame(sample_id = sample_id, method = "scDblFinder", cells_before = ncol(obj), predicted_doublets = 0L, removed_doublets = 0L, note = "Skipped: fewer than 100 cells or 20 genes")
    ))
  }
  sce <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder::scDblFinder(sce, dbr = params$doublet_rate)
  score <- as.numeric(SummarizedExperiment::colData(sce)[["scDblFinder.score"]])
  call <- as.character(SummarizedExperiment::colData(sce)[["scDblFinder.class"]]) == "doublet"
  obj$doublet_score <- score
  obj$predicted_doublet <- call
  calls <- data.frame(cell = colnames(obj), sample_id = sample_id, doublet_score = score, predicted_doublet = call,
                      removed_as_doublet = call & params$remove_doublets, stringsAsFactors = FALSE)
  summary <- data.frame(sample_id = sample_id, method = "scDblFinder", cells_before = ncol(obj),
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
  data.frame(cell = colnames(obj), sample_id = as.character(obj$sample_id), condition = as.character(obj$condition %||% ""),
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
save_plot(Seurat::VlnPlot(qc_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", ncol = 3, pt.size = 0), "01_qc_violin.png", 14, 5)
save_plot(Seurat::FeatureScatter(qc_merged, feature1 = "nCount_RNA", feature2 = "percent.mt") + Seurat::FeatureScatter(qc_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), "02_qc_scatter.png", 12, 5)
if (any(is.finite(doublet_calls$doublet_score))) {
  doublet_plot <- ggplot2::ggplot(doublet_calls[is.finite(doublet_calls$doublet_score), , drop = FALSE], ggplot2::aes(x = doublet_score)) +
    ggplot2::geom_histogram(bins = 40, fill = "#5b7db1", color = "white") +
    ggplot2::labs(title = "Doublet-score distribution before removal", x = "Doublet score", y = "Cells") +
    ggplot2::theme_classic()
  save_plot(doublet_plot, "02b_doublet_scores.png", 7, 4.5)
}

integration <- params$integration
batch_values <- unlist(lapply(objects, function(obj) {
  if (params$batch_column %in% colnames(obj@meta.data)) as.character(obj[[params$batch_column]][, 1]) else character(0)
}), use.names = FALSE)
if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "rpca" else "none"
if (integration %in% c("rpca", "cca") && length(unique(batch_values[nzchar(batch_values)])) < 2L) {
  stop(integration, " integration requires at least two values in the selected batch column (", params$batch_column, "). Choose none or supply the appropriate technical batch column.")
}
used_reduction <- "pca"
integration_k_weight <- max(5L, min(50L, floor(min(vapply(objects, ncol, numeric(1))) / 2)))
if (identical(params$normalization, "sct")) {
  objects <- lapply(objects, function(obj) Seurat::SCTransform(obj, assay = "RNA", vst.flavor = "v2", verbose = FALSE))
  if (integration %in% c("rpca", "cca") && length(objects) > 1L) {
    features <- Seurat::SelectIntegrationFeatures(objects, nfeatures = 3000)
    objects <- Seurat::PrepSCTIntegration(objects, anchor.features = features, verbose = FALSE)
    if (identical(integration, "rpca")) objects <- lapply(objects, function(obj) Seurat::RunPCA(obj, features = features, npcs = params$n_pcs, verbose = FALSE))
    anchors <- Seurat::FindIntegrationAnchors(object.list = objects, normalization.method = "SCT", anchor.features = features,
                                               reduction = integration, dims = seq_len(params$n_pcs), verbose = FALSE)
    obj <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = seq_len(params$n_pcs), k.weight = integration_k_weight, verbose = FALSE)
    DefaultAssay(obj) <- "integrated"
    obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  } else {
    obj <- Reduce(function(a, b) merge(a, y = b), objects)
    DefaultAssay(obj) <- "SCT"
    obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  }
} else {
  objects <- lapply(objects, function(obj) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    Seurat::ScaleData(obj, verbose = FALSE)
  })
  if (integration %in% c("rpca", "cca") && length(objects) > 1L) {
    features <- Seurat::SelectIntegrationFeatures(objects, nfeatures = 3000)
    if (identical(integration, "rpca")) objects <- lapply(objects, function(obj) Seurat::RunPCA(obj, features = features, npcs = params$n_pcs, verbose = FALSE))
    anchors <- Seurat::FindIntegrationAnchors(object.list = objects, anchor.features = features, reduction = integration,
                                               dims = seq_len(params$n_pcs), verbose = FALSE)
    obj <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(params$n_pcs), k.weight = integration_k_weight, verbose = FALSE)
    DefaultAssay(obj) <- "integrated"; obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  } else {
    obj <- Reduce(function(a, b) merge(a, y = b), objects)
    DefaultAssay(obj) <- "RNA"
    obj <- Seurat::NormalizeData(obj, verbose = FALSE); obj <- Seurat::FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE); obj <- Seurat::RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  }
}

n_available <- ncol(Seurat::Embeddings(obj, "pca")); dims <- seq_len(min(params$n_pcs, n_available))
pca_stdev <- Seurat::Stdev(obj, reduction = "pca")
pca_variance <- (pca_stdev ^ 2) / sum(pca_stdev ^ 2)
utils::write.table(data.frame(PC = seq_along(pca_variance), variance_explained = pca_variance,
                              percent_variance_explained = 100 * pca_variance),
                   file.path(tables_dir, "pca_variance_explained.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
variable_genes <- Seurat::VariableFeatures(obj)
utils::write.table(data.frame(gene = variable_genes, highly_variable = TRUE), file.path(tables_dir, "highly_variable_genes.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = dims, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = params$cluster_resolution, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = dims, seed.use = params$seed, verbose = FALSE)
obj$cluster <- as.character(Seurat::Idents(obj))

apply_celltype_mapping <- function(obj, path) {
  if (!nzchar(path) || identical(tolower(path), "none")) return(obj)
  map <- read_delim_safe(path)
  cell_col <- intersect(c("cell", "cell_id", "barcode", "Cell", "CellID"), names(map))
  type_col <- intersect(c("cell_type", "celltype", "CellType", "annotation"), names(map))
  if (!length(cell_col) || !length(type_col)) stop("Cell-type mapping needs a cell/barcode column and a cell_type column.")
  labels <- setNames(as.character(map[[type_col[[1]]]]), as.character(map[[cell_col[[1]]]]))
  direct <- unname(labels[colnames(obj)])
  raw_barcode <- sub("^[^_]+_", "", colnames(obj))
  direct[is.na(direct)] <- unname(labels[raw_barcode[is.na(direct)]])
  obj$cell_type <- ifelse(is.na(direct) | !nzchar(direct), "Unassigned", direct)
  obj$annotation_source <- "provided cell metadata"
  obj
}

apply_marker_annotation <- function(obj, path) {
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
  obj$cell_type <- unname(cluster_label[cluster])
  obj$annotation_source <- "marker-list cluster scoring"
  score_table <- data.frame(cluster = rownames(means), means, check.names = FALSE)
  utils::write.table(score_table, file.path(tables_dir, "marker_annotation_cluster_scores.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  obj
}

if (nzchar(params$celltype_file) && !identical(tolower(params$celltype_file), "none")) {
  obj <- apply_celltype_mapping(obj, params$celltype_file)
} else if (nzchar(params$marker_file) && !identical(tolower(params$marker_file), "none")) {
  obj <- apply_marker_annotation(obj, params$marker_file)
} else {
  obj$cell_type <- obj$cluster
  obj$annotation_source <- "cluster ID (no annotation supplied)"
}

save_plot(Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_id", shuffle = TRUE), "03_umap_sample.png", 8, 6)
save_plot(Seurat::DimPlot(obj, reduction = "umap", group.by = "cluster", label = TRUE, repel = TRUE), "04_umap_clusters.png", 8, 6)
save_plot(Seurat::DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE), "05_umap_cell_types.png", 9, 6)
if ("condition" %in% colnames(obj@meta.data) && length(unique(obj$condition)) > 1L) save_plot(Seurat::DimPlot(obj, reduction = "umap", group.by = "condition", shuffle = TRUE), "06_umap_condition.png", 8, 6)

DefaultAssay(obj) <- "RNA"
if (ncol(obj) >= 20 && length(unique(obj$cluster)) > 1L) {
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

cell_metadata <- data.frame(cell = colnames(obj), obj@meta.data, check.names = FALSE)
utils::write.table(cell_metadata, file.path(tables_dir, "cell_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cluster_sizes <- as.data.frame(table(cluster = obj$cluster, cell_type = obj$cell_type), stringsAsFactors = FALSE)
utils::write.table(cluster_sizes, file.path(tables_dir, "cluster_cell_type_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(obj, file.path(objects_dir, "processed_seurat.rds"))
summary_lines <- c(
  paste("engine: seurat"), paste("normalization:", params$normalization), paste("integration:", integration),
  paste("doublet_method:", params$doublet_method), paste("doublets_removed:", sum(doublet_summary$removed_doublets)),
  paste("input_samples:", NROW(samples)), paste("cells_after_qc:", ncol(obj)), paste("clusters:", length(unique(obj$cluster))),
  paste("annotation_source:", unique(obj$annotation_source)[1]), paste("generated:", as.character(Sys.time()))
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
writeLines(as.character(Sys.time()), file.path(out_dir, "_COMPLETE"))
