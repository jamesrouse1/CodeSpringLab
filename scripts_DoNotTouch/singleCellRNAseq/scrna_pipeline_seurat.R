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
  # Seurat's integration helpers still inventory future globals under a
  # sequential plan. For large inputs, the same object can be counted once as
  # object.list and again inside FUN even though it is not exported to another
  # process. Disable that validation ceiling only for this sequential job; this
  # does not create workers or duplicate the expression matrix in worker RAM.
  options(future.globals.maxSize = Inf)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x

allocated_cpu_count <- function(cap = Inf) {
  value <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
  if (is.na(value) || value < 1L) value <- 1L
  max(1L, min(value, as.integer(cap)))
}

# Seurat's uwot backend derives its thread count from future::nbrOfWorkers().
# Keep the global plan sequential for memory-heavy integration/transfer calls,
# but expose up to eight allocated cores only while UMAP runs. On Linux,
# multicore workers use forked memory and the UMAP call itself consumes the
# worker count as native threads; the full Seurat object is not exported to a
# persistent multisession cluster.
run_umap_with_allocated_threads <- function(object, ...) {
  workers <- allocated_cpu_count(8L)
  if (workers <= 1L || !requireNamespace("future", quietly = TRUE) || !isTRUE(future::supportsMulticore())) {
    return(Seurat::RunUMAP(object, ...))
  }
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::multicore, workers = workers)
  message("Running Seurat UMAP with ", workers, " allocated threads.")
  Seurat::RunUMAP(object, ...)
}

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
  object <- obj
  grouping <- group.by
  dots <- list(...)
  values <- if (!is.null(grouping) && nzchar(grouping) && grouping %in% colnames(object@meta.data)) object[[grouping]][, 1] else Seurat::Idents(object)
  colors <- categorical_palette(values)
  args <- c(list(object = object, group.by = grouping), if (length(colors)) list(cols = colors) else list(), dots)
  do.call(Seurat::DimPlot, args)
}
jpplot_point_layer <- function(plot, color = jpplot_colors[[7]]) {
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, "GeomPoint")) plot$layers[[i]]$aes_params$colour <- color
  }
  plot
}

read_delim_safe <- function(path) {
  reader <- if (grepl("\\.csv$", path, ignore.case = TRUE)) utils::read.csv else utils::read.delim
  tryCatch(reader(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = ""),
           error = function(e) stop("Could not read ", path, ": ", conditionMessage(e)))
}

normalized_column_name <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
pick_ortholog_column <- function(tab, candidates) {
  keys <- stats::setNames(names(tab), normalized_column_name(names(tab)))
  hit <- keys[normalized_column_name(candidates)]
  unname(hit[!is.na(hit)][1] %||% "")
}
read_ortholog_symbols <- function(path) {
  tab <- read_delim_safe(path)
  mouse_col <- pick_ortholog_column(tab, c("mouse_gene_symbol", "mouse_symbol", "mgi_symbol", "marker_symbol", "external_gene_name", "mouse_gene_name"))
  human_col <- pick_ortholog_column(tab, c("human_gene_symbol", "human_symbol", "hgnc_symbol", "human_gene_name", "hsapiens_homolog_associated_gene_name"))
  if (!nzchar(mouse_col) || !nzchar(human_col)) stop("Ortholog table needs recognizable mouse and human gene-symbol columns (for example mouse_gene_symbol and human_gene_symbol).")
  out <- unique(data.frame(mouse = trimws(as.character(tab[[mouse_col]])), human = trimws(as.character(tab[[human_col]])), stringsAsFactors = FALSE))
  out <- out[nzchar(out$mouse) & nzchar(out$human) & !is.na(out$mouse) & !is.na(out$human), , drop = FALSE]
  out
}
map_cross_species_genes <- function(genes, source_species, ortholog_path, expression_genes, audit_path, set_names = NULL) {
  genes <- trimws(as.character(genes)); keep <- nzchar(genes); genes <- genes[keep]
  if (!is.null(set_names)) set_names <- as.character(set_names)[keep]
  source_species <- tolower(source_species %||% "same")
  if (identical(source_species, "same")) return(genes)
  orth <- read_ortholog_symbols(ortholog_path)
  mouse_overlap <- sum(expression_genes %in% orth$mouse)
  human_overlap <- sum(expression_genes %in% orth$human)
  target_species <- if (mouse_overlap > human_overlap) "mouse" else if (human_overlap > mouse_overlap) "human" else "unknown"
  if (identical(target_species, "unknown")) stop("Could not infer whether expression features are mouse or human from the ortholog table. Choose 'Same as the expression dataset' if conversion is unnecessary.")
  if (identical(source_species, target_species)) {
    audit <- data.frame(set = set_names %||% "", original_gene = genes, mapped_gene = genes, status = "same_species", source_species = source_species, target_species = target_species)
    utils::write.table(audit, audit_path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(genes)
  }
  source_col <- source_species; target_col <- target_species
  source_counts <- table(orth[[source_col]])
  usable <- orth[source_counts[orth[[source_col]]] == 1L, , drop = FALSE]
  lookup <- stats::setNames(usable[[target_col]], usable[[source_col]])
  mapped <- unname(lookup[genes]); status <- ifelse(is.na(mapped) | !nzchar(mapped), "unmapped_or_ambiguous", "mapped")
  audit <- data.frame(set = set_names %||% "", original_gene = genes, mapped_gene = ifelse(is.na(mapped), "", mapped), status = status, source_species = source_species, target_species = target_species)
  utils::write.table(audit, audit_path, sep = "\t", row.names = FALSE, quote = FALSE)
  unique(mapped[status == "mapped"])
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
    batch_column = get("batch_column", "sample_id"),
    cluster_resolution = suppressWarnings(as.numeric(get("cluster_resolution", "0.6"))),
    n_pcs = suppressWarnings(as.integer(get("n_pcs", "30"))),
    n_neighbors = suppressWarnings(as.integer(get("n_neighbors", "15"))),
    umap_min_dist = suppressWarnings(as.numeric(get("umap_min_dist", "0.3"))),
    umap_spread = suppressWarnings(as.numeric(get("umap_spread", "1"))),
    umap_metric = tolower(get("umap_metric", "euclidean")),
    min_features = suppressWarnings(as.integer(get("min_features", "200"))),
    min_counts = suppressWarnings(as.integer(get("min_counts", "0"))),
    max_features = suppressWarnings(as.integer(get("max_features", "0"))),
    max_percent_mt = suppressWarnings(as.numeric(get("max_percent_mt", "20"))),
    min_cells_per_gene = suppressWarnings(as.integer(get("min_cells_per_gene", "3"))),
    doublet_method = tolower(get("doublet_method", "none")),
    # Zero means that scDblFinder should estimate the rate independently for
    # each capture. A positive value is an explicit global override.
    doublet_rate = suppressWarnings(as.numeric(get("doublet_rate", "0"))),
    remove_doublets = get_bool("remove_doublets", TRUE),
    marker_file = get("marker_file", ""),
    celltype_file = get("celltype_file", ""),
    reference_file = get("reference_file", ""),
    reference_label_column = get("reference_label_column", ""),
    reference_ortholog_file = get("reference_ortholog_file", ""),
    find_cluster_markers = tolower(get("find_cluster_markers", "false")) %in% c("true", "1", "yes"),
    marker_species = tolower(get("marker_species", "same")),
    marker_ortholog_file = get("marker_ortholog_file", ""),
    annotation_name = get("annotation_name", "cell_type"),
    signature_file = get("signature_file", ""),
    signature_species = tolower(get("signature_species", "same")),
    signature_ortholog_file = get("signature_ortholog_file", ""),
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
if (!is.finite(params$n_neighbors) || params$n_neighbors < 2L || params$n_neighbors > 200L) params$n_neighbors <- 15L
if (!is.finite(params$umap_min_dist) || params$umap_min_dist < 0 || params$umap_min_dist > 2) params$umap_min_dist <- 0.3
if (!is.finite(params$umap_spread) || params$umap_spread < 0.1 || params$umap_spread > 10) params$umap_spread <- 1
if (!params$umap_metric %in% c("euclidean", "cosine", "manhattan", "correlation")) params$umap_metric <- "euclidean"
if (!is.finite(params$min_features) || params$min_features < 0) params$min_features <- 200L
if (!is.finite(params$min_counts) || params$min_counts < 0) params$min_counts <- 0L
if (!is.finite(params$max_features) || params$max_features < 0) params$max_features <- 0L
if (!is.finite(params$max_percent_mt) || params$max_percent_mt < 0) params$max_percent_mt <- 20
if (!is.finite(params$min_cells_per_gene) || params$min_cells_per_gene < 1) params$min_cells_per_gene <- 3L
if (!params$doublet_method %in% c("auto", "none", "scdblfinder")) stop("doublet_method must be auto, none, or scDblFinder")
if (!is.finite(params$doublet_rate) || params$doublet_rate < 0 || params$doublet_rate >= 1) params$doublet_rate <- 0
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
if (!"capture_id" %in% names(samples)) samples$capture_id <- samples$sample_id
samples$capture_id <- trimws(as.character(samples$capture_id))
samples$capture_id[!nzchar(samples$capture_id)] <- samples$sample_id[!nzchar(samples$capture_id)]
if (any(!nzchar(samples$capture_id))) stop("Every sample needs a capture_id or sample_id.")
if ("expected_doublet_rate" %in% names(samples)) {
  raw_rate <- trimws(as.character(samples$expected_doublet_rate))
  supplied <- nzchar(raw_rate) & !is.na(raw_rate)
  parsed <- suppressWarnings(as.numeric(raw_rate))
  if (any(supplied & (!is.finite(parsed) | parsed <= 0 | parsed >= 1))) stop("expected_doublet_rate must be blank for automatic estimation or greater than 0 and less than 1.")
  rates_by_capture <- split(parsed[supplied], samples$capture_id[supplied])
  inconsistent <- names(Filter(function(x) length(unique(x)) > 1L, rates_by_capture))
  if (length(inconsistent)) stop("expected_doublet_rate must be constant within capture_id: ", inconsistent[[1]])
}
samples$input_path <- normalizePath(path.expand(as.character(samples$input_path)), winslash = "/", mustWork = FALSE)
missing <- samples$input_path[!file.exists(samples$input_path)]
if (length(missing)) stop("Missing unreadable input path(s): ", paste(missing, collapse = ", "))

tables_dir <- file.path(out_dir, "tables"); figures_dir <- file.path(out_dir, "figures"); objects_dir <- file.path(out_dir, "objects"); checkpoints_dir <- file.path(out_dir, "checkpoints")
for (d in c(tables_dir, figures_dir, objects_dir, checkpoints_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
utils::write.table(samples, file.path(tables_dir, "sample_manifest_used.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
stage_marker <- function(name) writeLines("complete", file.path(out_dir, paste0("_STAGE_", toupper(name), "_COMPLETE")))
checkpoint_path <- function(name) file.path(checkpoints_dir, paste0(name, "_seurat.rds"))
atomic_save_rds <- function(object, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (!file.exists(temporary) || file.info(temporary)$size <= 0) stop("Temporary RDS write was empty: ", temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically replace the processed Seurat object: ", path)
  invisible(path)
}
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
    workflow_action = "Raw RNA counts retained; a complete existing UMAP and cluster assignment can be continued directly, while reconstruction remains optional.",
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

# Keep extreme mitochondrial-content droplets from flattening the useful QC
# range. This affects only the displayed axis; exported values and filtering
# always use the original, uncapped percent.mt measurements.
mt_display_limit <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  percentile_99 <- if (length(values)) stats::quantile(values, 0.99, names = FALSE, na.rm = TRUE) else 10
  max(10, min(50, ceiling(percentile_99 / 5) * 5))
}

round_qc_threshold <- function(value, direction = c("down", "up"), minimum = 0) {
  direction <- match.arg(direction)
  value <- as.numeric(value)
  if (!is.finite(value)) return(as.integer(minimum))
  increment <- if (abs(value) >= 100) 100 else 10
  rounded <- if (identical(direction, "down")) floor(value / increment) * increment else ceiling(value / increment) * increment
  as.integer(max(minimum, rounded))
}

qc_recommendations <- function(metrics) {
  robust_tail <- function(values, direction = c("lower", "upper"), log_scale = FALSE) {
    direction <- match.arg(direction)
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (!length(values)) return(NA_real_)
    working <- if (log_scale) log1p(values) else values
    center <- stats::median(working)
    spread <- stats::mad(working, center = center, constant = 1.4826)
    sign <- if (identical(direction, "lower")) -1 else 1
    threshold <- center + sign * 3 * spread
    if (!is.finite(threshold) || spread == 0) {
      threshold <- stats::quantile(working, if (identical(direction, "lower")) 0.01 else 0.99, names = FALSE, na.rm = TRUE)
    }
    if (log_scale) exp(threshold) - 1 else threshold
  }
  rows <- lapply(split(metrics, metrics$sample_id), function(sample) {
    low_genes <- robust_tail(sample$nFeature_RNA, "lower", log_scale = TRUE)
    low_counts <- robust_tail(sample$nCount_RNA, "lower", log_scale = TRUE)
    high_genes <- robust_tail(sample$nFeature_RNA, "upper", log_scale = TRUE)
    high_mt <- robust_tail(sample$percent.mt, "upper", log_scale = FALSE)
    data.frame(
      sample_id = as.character(sample$sample_id[[1]]),
      cells_before_qc = NROW(sample),
      min_features = round_qc_threshold(low_genes, "down", minimum = 200),
      min_counts = round_qc_threshold(low_counts, "down", minimum = 0),
      max_features = round_qc_threshold(high_genes, "up", minimum = 0),
      max_percent_mt = as.integer(max(5, min(50, ceiling(high_mt)))),
      recommendation_basis = "Per-sample 3-MAD outlier thresholds; counts and genes calculated on log1p scale; mitochondrial percentage on original scale; editable starting values",
      stringsAsFactors = FALSE
    )
  })
  recommendations <- do.call(rbind, rows)
  if (is.null(recommendations) || !NROW(recommendations)) return(data.frame())
  global <- data.frame(
    sample_id = "Recommended global",
    cells_before_qc = sum(recommendations$cells_before_qc),
    min_features = min(recommendations$min_features),
    min_counts = min(recommendations$min_counts),
    max_features = max(recommendations$max_features),
    max_percent_mt = max(recommendations$max_percent_mt),
    recommendation_basis = "Conservative shared starting values across samples; inspect per-sample distributions before applying",
    stringsAsFactors = FALSE
  )
  rbind(global, recommendations)
}

qc_cutoff_lines <- function(plot, values) {
  values <- as.numeric(values)
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) return(plot)
  plot + ggplot2::geom_hline(yintercept = values, color = "#C62828", linetype = "dashed", linewidth = 0.55)
}

qc_violin_plot <- function(obj, cutoffs = NULL, state_label = "") {
  plots <- Seurat::VlnPlot(
    obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample_id", cols = categorical_palette(obj$sample_id),
    ncol = 3, pt.size = 0, layer = "counts", combine = FALSE
  )
  mt_cap <- mt_display_limit(obj$percent.mt)
  panel_titles <- c("Genes per cell", "Counts per cell", "Mitochondrial reads (%)")
  for (i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] +
      ggplot2::labs(title = panel_titles[[i]], x = "Sample", y = panel_titles[[i]]) +
      ggplot2::theme(legend.position = "none")
  }
  if (!is.null(cutoffs)) {
    plots[[1]] <- qc_cutoff_lines(plots[[1]], c(cutoffs$min_features, cutoffs$max_features))
    plots[[2]] <- qc_cutoff_lines(plots[[2]], cutoffs$min_counts)
    plots[[3]] <- qc_cutoff_lines(plots[[3]], cutoffs$max_percent_mt)
  }
  plots[[3]] <- plots[[3]] +
    ggplot2::coord_cartesian(ylim = c(0, mt_cap)) +
    ggplot2::labs(caption = paste0("Display capped at ", mt_cap, "% (99th-percentile rule; maximum 50%)."))
  combined <- patchwork::wrap_plots(plots, ncol = 3)
  if (nzchar(state_label)) combined <- combined + patchwork::plot_annotation(title = state_label)
  combined
}

qc_scatter_plot <- function(obj, cutoffs = NULL, state_label = "") {
  mt_cap <- mt_display_limit(obj$percent.mt)
  counts_mt <- jpplot_point_layer(Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")) +
    ggplot2::coord_cartesian(ylim = c(0, mt_cap)) +
    ggplot2::labs(title = "Counts versus mitochondrial reads", x = "Counts per cell", y = "Mitochondrial reads (%)", caption = paste0("Mitochondrial axis capped at ", mt_cap, "% for display only."))
  counts_features <- jpplot_point_layer(Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")) +
    ggplot2::labs(title = "Counts versus detected genes", x = "Counts per cell", y = "Genes per cell")
  if (!is.null(cutoffs)) {
    if (is.finite(as.numeric(cutoffs$min_counts)) && as.numeric(cutoffs$min_counts) > 0) {
      counts_mt <- counts_mt + ggplot2::geom_vline(xintercept = as.numeric(cutoffs$min_counts), color = "#C62828", linetype = "dashed", linewidth = 0.55)
      counts_features <- counts_features + ggplot2::geom_vline(xintercept = as.numeric(cutoffs$min_counts), color = "#C62828", linetype = "dashed", linewidth = 0.55)
    }
    counts_mt <- qc_cutoff_lines(counts_mt, cutoffs$max_percent_mt)
    counts_features <- qc_cutoff_lines(counts_features, c(cutoffs$min_features, cutoffs$max_features))
  }
  combined <- patchwork::wrap_plots(counts_mt, counts_features, ncol = 2)
  if (nzchar(state_label)) combined <- combined + patchwork::plot_annotation(title = state_label)
  combined
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
  recommended_qc <- qc_recommendations(pre_qc_cells)
  utils::write.table(recommended_qc, file.path(tables_dir, "qc_recommended_thresholds.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  global_qc <- recommended_qc[recommended_qc$sample_id == "Recommended global", , drop = FALSE]
  suggested_cutoffs <- if (NROW(global_qc)) as.list(global_qc[1, c("min_features", "min_counts", "max_features", "max_percent_mt")]) else NULL
  pre_qc_merged <- Reduce(function(a, b) merge(a, y = b), objects)
  Seurat::DefaultAssay(pre_qc_merged) <- "RNA"
  save_plot(qc_violin_plot(pre_qc_merged, suggested_cutoffs, "Unfiltered cells — suggested starting cutoffs"), "00_qc_pre_filter_violin.png", 14, 5)
  save_plot(qc_scatter_plot(pre_qc_merged, suggested_cutoffs, "Unfiltered cells — suggested starting cutoffs"), "00_qc_pre_filter_scatter.png", 12, 5)
  # A single processed Seurat object can be continued without rebuilding its
  # existing embedding. Preserve a project-local copy only when both UMAP and
  # clusters were detected; source data remain read-only and unchanged.
  reusable_processed_input <- NROW(samples) == 1L &&
    identical(input_processing_status$input_kind[[1]], "seurat_rds") &&
    isTRUE(input_processing_status$umap_detected[[1]]) &&
    isTRUE(input_processing_status$clusters_detected[[1]])
  if (reusable_processed_input) {
    if (!"cluster" %in% colnames(pre_qc_merged@meta.data)) {
      cluster_source <- intersect(c("seurat_clusters", "cluster"), colnames(pre_qc_merged@meta.data))
      if (length(cluster_source)) pre_qc_merged$cluster <- as.character(pre_qc_merged[[cluster_source[[1]]]][, 1])
    }
    attr(pre_qc_merged, "codespring_integration") <- "existing input object"
    attr(pre_qc_merged, "codespring_doublets_removed") <- 0
    atomic_save_rds(pre_qc_merged, file.path(objects_dir, "processed_seurat.rds"), compress = FALSE)
  }
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
capture_doublet_rate <- function(obj) {
  if (!"expected_doublet_rate" %in% colnames(obj@meta.data)) return(if (params$doublet_rate > 0) params$doublet_rate else NULL)
  raw <- unique(trimws(as.character(obj$expected_doublet_rate)))
  raw <- raw[nzchar(raw) & !is.na(raw)]
  if (!length(raw)) return(if (params$doublet_rate > 0) params$doublet_rate else NULL)
  if (length(raw) != 1L) stop("expected_doublet_rate must have one value per capture.")
  value <- suppressWarnings(as.numeric(raw[[1]]))
  if (!is.finite(value) || value <= 0 || value >= 1) stop("expected_doublet_rate must be greater than 0 and less than 1 for each capture.")
  value
}

run_doublet_one <- function(obj) {
  capture_id <- unique(as.character(obj$capture_id))
  if (length(capture_id) != 1L) stop("A doublet partition contains multiple capture_id values.")
  capture_id <- capture_id[[1]]
  sample_ids <- paste(sort(unique(as.character(obj$sample_id))), collapse = ",")
  expected_rate <- capture_doublet_rate(obj)
  manifest_rate_supplied <- FALSE
  if ("expected_doublet_rate" %in% colnames(obj@meta.data)) {
    manifest_values <- trimws(as.character(obj$expected_doublet_rate))
    manifest_rate_supplied <- any(!is.na(manifest_values) & nzchar(manifest_values))
  }
  rate_source <- if (is.null(expected_rate)) "automatic" else if (manifest_rate_supplied) "manifest" else "global override"
  reported_rate <- if (is.null(expected_rate)) NA_real_ else expected_rate
  obj$doublet_score <- NA_real_
  obj$predicted_doublet <- FALSE
  calls_frame <- function(score, call, removed) data.frame(
    cell = colnames(obj), sample_id = as.character(obj$sample_id), capture_id = capture_id,
    doublet_score = score, predicted_doublet = call, removed_as_doublet = removed,
    stringsAsFactors = FALSE
  )
  summary_frame <- function(method, predicted = 0L, removed = 0L, note = "") data.frame(
    capture_id = capture_id, sample_ids = sample_ids, method = method,
    expected_doublet_rate = reported_rate, rate_source = rate_source,
    cells_before = ncol(obj), predicted_doublets = predicted,
    removed_doublets = removed, note = note, stringsAsFactors = FALSE
  )
  if (identical(params$doublet_method, "none")) {
    return(list(
      object = obj,
      calls = calls_frame(NA_real_, FALSE, FALSE),
      summary = summary_frame("none", note = "Doublet removal disabled")
    ))
  }
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    if (identical(params$doublet_method, "scdblfinder")) {
      stop("Doublet removal with the Seurat engine requires the Bioconductor package scDblFinder. Install it, choose the Scanpy/Scrublet engine, or explicitly select no doublet removal.")
    }
    return(list(
      object = obj,
      calls = calls_frame(NA_real_, FALSE, FALSE),
      summary = summary_frame("scDblFinder", note = "Automatic scDblFinder skipped: package unavailable")
    ))
  }
  if (ncol(obj) < 100L || nrow(obj) < 20L) {
    return(list(
      object = obj,
      calls = calls_frame(NA_real_, FALSE, FALSE),
      summary = summary_frame("scDblFinder", note = "Skipped: fewer than 100 cells or 20 genes")
    ))
  }
  sce <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  sce <- if (is.null(expected_rate)) scDblFinder::scDblFinder(sce) else scDblFinder::scDblFinder(sce, dbr = expected_rate)
  score <- as.numeric(SummarizedExperiment::colData(sce)[["scDblFinder.score"]])
  call <- as.character(SummarizedExperiment::colData(sce)[["scDblFinder.class"]]) == "doublet"
  obj$doublet_score <- score
  obj$predicted_doublet <- call
  calls <- calls_frame(score, call, call & params$remove_doublets)
  summary <- summary_frame("scDblFinder", sum(call), if (params$remove_doublets) sum(call) else 0L)
  if (params$remove_doublets) obj <- subset(obj, cells = colnames(obj)[!call])
  if (!ncol(obj)) stop("Doublet removal removed every cell in capture ", capture_id, ".")
  list(object = obj, calls = calls, summary = summary)
}

# Remove only extremely low-coverage barcodes before doublet detection. Full
# mitochondrial/count/feature QC is deliberately performed afterward, per the
# scDblFinder guidance, so potential low-quality members of true doublets are
# still available to the classifier.
minimal_prefilter_one <- function(obj) {
  keep <- obj$nCount_RNA >= 200
  if (!any(keep)) stop("The minimal 200-count prefilter removed every cell in sample ", unique(obj$sample_id)[1], ".")
  subset(obj, cells = colnames(obj)[keep])
}
objects <- lapply(objects, minimal_prefilter_one)
capture_for_object <- vapply(objects, function(obj) unique(as.character(obj$capture_id))[[1]], character(1))
capture_objects <- lapply(split(objects, capture_for_object), function(parts) {
  obj <- if (length(parts) == 1L) parts[[1]] else Reduce(function(a, b) merge(a, y = b), parts)
  count_layers <- grep("^counts($|\\.)", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
  if (length(count_layers) > 1L) obj <- SeuratObject::JoinLayers(obj, assay = "RNA", layers = "counts", new = "counts")
  obj
})

# Each independent droplet capture is processed separately. Biological samples
# multiplexed in the same capture therefore share one realistic doublet model.
doublet_runs <- lapply(capture_objects, run_doublet_one)
clean_capture_objects <- lapply(doublet_runs, `[[`, "object")
doublet_calls <- do.call(rbind, lapply(doublet_runs, `[[`, "calls"))
doublet_summary <- do.call(rbind, lapply(doublet_runs, `[[`, "summary"))
utils::write.table(doublet_calls, file.path(tables_dir, "doublet_calls.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(doublet_summary, file.path(tables_dir, "doublet_summary_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(doublet_summary, file.path(tables_dir, "doublet_summary_by_capture.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Split multiplexed captures back into biological samples, then apply the full
# user-reviewed QC thresholds independently to every sample.
objects <- lapply(samples$sample_id, function(sample_id) {
  capture_id <- samples$capture_id[match(sample_id, samples$sample_id)]
  obj <- clean_capture_objects[[capture_id]]
  keep <- as.character(obj$sample_id) == sample_id
  if (!any(keep)) stop("No cells remain for sample ", sample_id, " after doublet detection.")
  subset(obj, cells = colnames(obj)[keep])
})
names(objects) <- samples$sample_id
objects <- lapply(objects, qc_one)

filter_genes_one <- function(obj) {
  counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  detected_cells <- Matrix::rowSums(counts > 0)
  keep <- detected_cells >= params$min_cells_per_gene
  if (!any(keep)) stop("Minimum cells per retained gene removed every gene in sample ", unique(obj$sample_id)[1], ".")
  summary <- data.frame(
    sample_id = unique(as.character(obj$sample_id))[[1]], genes_before_filtering = nrow(obj),
    min_cells_per_gene = params$min_cells_per_gene, genes_retained = sum(keep), stringsAsFactors = FALSE
  )
  list(object = subset(obj, features = rownames(obj)[keep]), summary = summary)
}
feature_filter_runs <- lapply(objects, filter_genes_one)
objects <- lapply(feature_filter_runs, `[[`, "object")
feature_filter_summary <- do.call(rbind, lapply(feature_filter_runs, `[[`, "summary"))
utils::write.table(feature_filter_summary, file.path(tables_dir, "feature_filtering_by_sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

qc_cells <- do.call(rbind, lapply(objects, function(obj) {
  condition_value <- if ("condition" %in% colnames(obj@meta.data)) as.character(obj[["condition"]][, 1]) else rep("", ncol(obj))
  data.frame(cell = colnames(obj), sample_id = as.character(obj$sample_id), capture_id = as.character(obj$capture_id), condition = condition_value,
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
applied_cutoffs <- list(min_features = params$min_features, min_counts = params$min_counts, max_features = params$max_features, max_percent_mt = params$max_percent_mt)
# Use the same names as the Scanpy runner so the app can display retained-cell
# QC figures uniformly regardless of the selected analysis engine.
save_plot(qc_violin_plot(qc_merged, applied_cutoffs, "Retained cells — applied QC cutoffs"), "01_qc_post_filter_violin.png", 14, 5)
save_plot(qc_scatter_plot(qc_merged, applied_cutoffs, "Retained cells — applied QC cutoffs"), "01_qc_post_filter_scatter.png", 12, 5)
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

if (stage %in% c("preprocess", "all")) {
integration <- params$integration
batch_values <- unlist(lapply(objects, function(obj) {
  if (params$batch_column %in% colnames(obj@meta.data)) as.character(obj[[params$batch_column]][, 1]) else character(0)
}), use.names = FALSE)
if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "harmony" else "none"
if (integration %in% c("rpca", "cca", "harmony") && length(unique(batch_values[nzchar(batch_values)])) < 2L) {
  stop(integration, " integration requires at least two values in the selected batch column (", params$batch_column, "). Choose none or supply the appropriate technical batch column.")
}
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
  pre_neighbors <- min(params$n_neighbors, max(2L, ncol(unintegrated) - 1L))
  unintegrated <- Seurat::FindNeighbors(unintegrated, reduction = "pca", dims = pre_dims, k.param = pre_neighbors, verbose = FALSE)
  unintegrated <- run_umap_with_allocated_threads(unintegrated, reduction = "pca", dims = pre_dims, reduction.name = "umap.unintegrated",
                                  n.neighbors = pre_neighbors, min.dist = params$umap_min_dist, spread = params$umap_spread,
                                  metric = params$umap_metric,
                                  seed.use = params$seed, verbose = FALSE)
  save_plot(Seurat::DimPlot(unintegrated, reduction = "umap.unintegrated", group.by = "sample_id", shuffle = TRUE), "02_preintegration_umap_sample.png", 8, 6)
  if (nzchar(params$batch_column) && params$batch_column %in% colnames(unintegrated@meta.data) && length(unique(as.character(unintegrated[[params$batch_column]][, 1]))) > 1L && !identical(params$batch_column, "sample_id")) {
    save_plot(Seurat::DimPlot(unintegrated, reduction = "umap.unintegrated", group.by = params$batch_column, shuffle = TRUE), "02_preintegration_umap_batch.png", 8, 6)
  }
  pre_coords <- as.data.frame(Seurat::Embeddings(unintegrated, "umap.unintegrated"))
  names(pre_coords)[1:2] <- c("UMAP_1", "UMAP_2")
  pre_metadata <- unintegrated@meta.data
  if ("cell" %in% names(pre_metadata)) names(pre_metadata)[names(pre_metadata) == "cell"] <- "input_cell"
  pre_table <- data.frame(cell = rownames(pre_coords), pre_coords[, c("UMAP_1", "UMAP_2"), drop = FALSE], pre_metadata, check.names = FALSE)
  utils::write.table(pre_table, file.path(tables_dir, "preintegration_umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
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
if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "harmony" else "none"
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
  graph_neighbors <- min(params$n_neighbors, max(2L, ncol(obj) - 1L))
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_for_graph, dims = dims, k.param = graph_neighbors, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = params$cluster_resolution, verbose = FALSE)
  obj <- run_umap_with_allocated_threads(obj, reduction = reduction_for_graph, dims = dims,
                         n.neighbors = graph_neighbors, min.dist = params$umap_min_dist, spread = params$umap_spread,
                         metric = params$umap_metric,
                         seed.use = params$seed, verbose = FALSE)
  obj$cluster <- as.character(Seurat::Idents(obj))
  # Publish the final/integrated embedding immediately after clustering so it
  # is visible in both Run Pipeline and the interactive Results Explorer. Cell
  # type annotations can be joined later without recomputing the coordinates.
  save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "cluster", label = TRUE, repel = TRUE), "04_umap_clusters_pre_annotation.png", 8, 6)
  if (length(unique(as.character(obj$sample_id))) > 1L) {
    save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "sample_id", shuffle = TRUE), "04_umap_samples_pre_annotation.png", 8, 6)
  }
  cluster_umap <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"), check.names = FALSE)
  colnames(cluster_umap)[seq_len(min(2L, NCOL(cluster_umap)))] <- c("UMAP_1", "UMAP_2")[seq_len(min(2L, NCOL(cluster_umap)))]
  if (!all(c("UMAP_1", "UMAP_2") %in% names(cluster_umap))) stop("The final UMAP does not contain two coordinates.")
  cluster_metadata <- obj@meta.data
  if ("cell" %in% names(cluster_metadata)) names(cluster_metadata)[names(cluster_metadata) == "cell"] <- "input_cell"
  cluster_umap_table <- data.frame(cell = rownames(cluster_umap), cluster_umap[, c("UMAP_1", "UMAP_2"), drop = FALSE], cluster_metadata, check.names = FALSE)
  utils::write.table(cluster_umap_table, file.path(tables_dir, "umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(list(object = obj, integration = integration, samples = samples, doublet_summary = doublet_summary), checkpoint_path("04_clustered"))
  # The interactive dashboard reads this normalized checkpoint directly until
  # annotation or signature scoring creates a modified final object. Keep the
  # small gene list beside it without duplicating a multi-gigabyte RDS.
  utils::write.table(data.frame(gene = rownames(obj)), file.path(tables_dir, "dashboard_all_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  stage_marker("cluster")
  if (identical(stage, "cluster")) quit(save = "no", status = 0L)
} else if (identical(stage, "annotate")) {
  processed_path <- file.path(objects_dir, "processed_seurat.rds")
  clustered_path <- checkpoint_path("04_clustered")
  use_processed <- file.exists(processed_path) && (!file.exists(clustered_path) || file.info(processed_path)$mtime >= file.info(clustered_path)$mtime)
  if (use_processed) {
    processed_result <- tryCatch(readRDS(processed_path), error = function(e) e)
    if (!inherits(processed_result, "error")) {
      obj <- processed_result
      message("Loaded existing processed Seurat object for annotation.")
      integration <- attr(obj, "codespring_integration") %||% params$integration
      doublet_summary <- data.frame(removed_doublets = attr(obj, "codespring_doublets_removed") %||% 0)
    } else if (file.exists(clustered_path)) {
      message("Existing processed Seurat object is unreadable; recovering from the intact clustered checkpoint.")
      cluster_state <- require_checkpoint("04_clustered", "integration and clustering")
      obj <- cluster_state$object
      integration <- cluster_state$integration
      samples <- cluster_state$samples
      doublet_summary <- cluster_state$doublet_summary
    } else if (NROW(samples) == 1L && identical(input_kind(samples$input_path[[1]]), "seurat_rds")) {
      message("Existing processed Seurat object is unreadable; recovering from the original read-only Seurat input.")
      obj <- read_one(samples[1, , drop = FALSE])
      integration <- "existing input object"
      doublet_summary <- data.frame(removed_doublets = 0)
    } else {
      stop("The processed Seurat object is unreadable and no intact clustered checkpoint is available: ", conditionMessage(processed_result))
    }
  } else {
    cluster_state <- require_checkpoint("04_clustered", "integration and clustering")
    obj <- cluster_state$object
    integration <- cluster_state$integration
    samples <- cluster_state$samples
    doublet_summary <- cluster_state$doublet_summary
  }
} else if (stage %in% c("score", "differential")) {
  processed_path <- file.path(objects_dir, "processed_seurat.rds")
  clustered_path <- checkpoint_path("04_clustered")
  use_processed <- file.exists(processed_path) && (!file.exists(clustered_path) || file.info(processed_path)$mtime >= file.info(clustered_path)$mtime)
  if (use_processed) {
    obj <- readRDS(processed_path)
    integration <- attr(obj, "codespring_integration") %||% params$integration
    doublet_summary <- data.frame(removed_doublets = attr(obj, "codespring_doublets_removed") %||% 0)
  } else if (identical(stage, "score") && file.exists(clustered_path)) {
    cluster_state <- require_checkpoint("04_clustered", "UMAP and clustering")
    obj <- cluster_state$object
    integration <- cluster_state$integration
    doublet_summary <- cluster_state$doublet_summary
  } else {
    stop("The annotation stage has not completed. Run annotation before ", stage, ".")
  }
}

# Older uploaded Seurat objects commonly store clusters only as
# `seurat_clusters` or active identities. Keep the pipeline's stable `cluster`
# alias available after either a normal load or recovery from the source file.
if (!"cluster" %in% colnames(obj@meta.data)) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    obj$cluster <- as.character(obj$seurat_clusters)
  } else {
    obj$cluster <- as.character(Seurat::Idents(obj))
  }
}

safe_metadata_name <- function(value, default = "cell_type") {
  value <- gsub("[^A-Za-z0-9_]+", "_", trimws(as.character(value)))
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) value <- default
  if (grepl("^[0-9]", value)) value <- paste0("annotation_", value)
  value
}

write_interactive_metadata_tables <- function(obj) {
  cell_metadata <- obj@meta.data
  if ("cell" %in% names(cell_metadata)) names(cell_metadata)[names(cell_metadata) == "cell"] <- "input_cell"
  cell_metadata <- data.frame(cell = colnames(obj), cell_metadata, check.names = FALSE)
  utils::write.table(cell_metadata, file.path(tables_dir, "cell_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  if (!"umap" %in% names(obj@reductions)) return(invisible(FALSE))
  umap_coordinates <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"), check.names = FALSE)
  colnames(umap_coordinates)[seq_len(min(2L, NCOL(umap_coordinates)))] <- c("UMAP_1", "UMAP_2")[seq_len(min(2L, NCOL(umap_coordinates)))]
  if (!all(c("UMAP_1", "UMAP_2") %in% names(umap_coordinates))) stop("The final UMAP does not contain two coordinates.")
  umap_metadata <- obj@meta.data[rownames(umap_coordinates), , drop = FALSE]
  if ("cell" %in% names(umap_metadata)) names(umap_metadata)[names(umap_metadata) == "cell"] <- "input_cell"
  umap_table <- data.frame(cell = rownames(umap_coordinates), umap_coordinates[, c("UMAP_1", "UMAP_2"), drop = FALSE], umap_metadata, check.names = FALSE)
  utils::write.table(umap_table, file.path(tables_dir, "umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
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
  expanded <- do.call(rbind, lapply(seq_len(NROW(markers)), function(i) {
    genes <- unique(Filter(nzchar, strsplit(as.character(markers$gene[[i]]), "[,;[:space:]]+")[[1]]))
    data.frame(cell_type = as.character(markers$cell_type[[i]]), gene = genes, stringsAsFactors = FALSE)
  }))
  if (!identical(params$marker_species, "same")) {
    mapped <- map_cross_species_genes(expanded$gene, params$marker_species, params$marker_ortholog_file, rownames(obj), file.path(tables_dir, "marker_ortholog_mapping.tsv"), expanded$cell_type)
    audit <- read_delim_safe(file.path(tables_dir, "marker_ortholog_mapping.tsv"))
    expanded <- data.frame(cell_type = audit$set[audit$status %in% c("mapped", "same_species")], gene = audit$mapped_gene[audit$status %in% c("mapped", "same_species")], stringsAsFactors = FALSE)
  }
  marker_list <- split(expanded$gene, expanded$cell_type)
  marker_list <- lapply(marker_list, function(x) intersect(unique(x[nzchar(x)]), rownames(obj)))
  marker_list <- marker_list[lengths(marker_list) > 0]
  if (!length(marker_list)) stop("None of the supplied marker genes were present in the expression object.")
  # Retain the exact usable markers so the final dot plots and heatmaps always
  # match the genes that were actually used for scoring.
  attr(obj, "codespring_marker_gene_sets") <- marker_list
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

marker_list_panels <- function(marker_list, max_cell_types = 6L, max_genes = 32L) {
  marker_list <- marker_list[lengths(marker_list) > 0]
  panels <- list(); current <- list(); current_genes <- 0L
  add_panel <- function() {
    if (length(current)) panels[[length(panels) + 1L]] <<- current
    current <<- list(); current_genes <<- 0L
  }
  for (label in names(marker_list)) {
    genes <- unique(as.character(marker_list[[label]]))
    # Very large single labels are split without mixing them with another
    # label, preserving legible axes rather than producing a compressed plot.
    blocks <- split(genes, ceiling(seq_along(genes) / max_genes))
    for (block_index in seq_along(blocks)) {
      block <- blocks[[block_index]]
      block_label <- if (length(blocks) > 1L) paste0(label, " (part ", block_index, ")") else label
      if (length(current) && (length(current) >= max_cell_types || current_genes + length(block) > max_genes)) add_panel()
      current[[block_label]] <- block
      current_genes <- current_genes + length(block)
    }
  }
  add_panel()
  panels
}

save_marker_annotation_panels <- function(obj, marker_list, annotation_name) {
  if (!length(marker_list)) return(invisible(FALSE))
  assay <- DefaultAssay(obj)
  expression <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
  groupings <- list(cluster = as.character(obj$cluster))
  if (annotation_name %in% colnames(obj@meta.data)) groupings$cell_type <- as.character(obj[[annotation_name]][, 1])
  panels <- marker_list_panels(marker_list)
  for (panel_index in seq_along(panels)) {
    features <- panels[[panel_index]]
    # Seurat::DotPlot uses factor levels for feature labels and rejects a
    # gene repeated under two biologically related labels. Shared genes still
    # contribute to both annotation scores; show each once in the visual.
    seen_genes <- character(0)
    dot_features <- lapply(features, function(set) {
      unique_set <- unique(setdiff(set, seen_genes))
      seen_genes <<- unique(c(seen_genes, set))
      unique_set
    })
    dot_features <- dot_features[lengths(dot_features) > 0]
    marker_rows <- do.call(rbind, lapply(names(features), function(label) data.frame(marker_set = label, gene = as.character(features[[label]]), stringsAsFactors = FALSE)))
    marker_rows <- marker_rows[!duplicated(marker_rows$gene) & marker_rows$gene %in% rownames(expression), , drop = FALSE]
    genes <- marker_rows$gene
    if (!length(genes)) next
    panel_label <- paste0("Marker-list annotation: ", paste(names(features), collapse = "; "))
    for (group_key in names(groupings)) {
      grouping <- groupings[[group_key]]
      group_levels <- if (identical(group_key, "cluster")) unique(grouping[order(suppressWarnings(as.numeric(grouping)), grouping, na.last = TRUE)]) else sort(unique(grouping))
      group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
      if (!length(group_levels)) next
      group_column <- if (identical(group_key, "cluster")) "cluster" else annotation_name
      group_label <- if (identical(group_key, "cluster")) "Cluster" else "Cell type"
      width <- max(8, 3.8 + 0.58 * length(group_levels)); height <- max(5.5, 2.8 + 0.22 * length(genes))
      dot <- Seurat::DotPlot(obj, features = dot_features, group.by = group_column, assay = assay, cols = c("#E9F2FA", "#B2182B")) +
        Seurat::RotatedAxis() +
        ggplot2::labs(title = panel_label, subtitle = paste0("Grouped by ", tolower(group_label), "; shared markers are displayed once"), x = "Marker gene", y = group_label, color = "Average expression", size = "% expressing") +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      save_plot(dot, sprintf("07_marker_annotation_dotplot_by_%s_panel_%02d.png", group_key, panel_index), width, height)

      means <- do.call(cbind, lapply(group_levels, function(group) Matrix::rowMeans(expression[genes, grouping == group, drop = FALSE])))
      row_labels <- make.unique(paste0(marker_rows$marker_set, " | ", genes))
      means <- as.matrix(means); rownames(means) <- row_labels; colnames(means) <- group_levels
      scaled <- t(scale(t(means))); scaled[!is.finite(scaled)] <- 0
      scaled <- matrix(pmax(-2.5, pmin(2.5, scaled)), nrow = NROW(means), ncol = NCOL(means), dimnames = dimnames(means))
      heatmap_path <- file.path(figures_dir, sprintf("07_marker_annotation_heatmap_by_%s_panel_%02d.png", group_key, panel_index))
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        annotation_row <- data.frame(`Marker set` = marker_rows$marker_set, row.names = rownames(scaled), check.names = FALSE)
        marker_colors <- categorical_palette(marker_rows$marker_set)
        group_breaks <- which(marker_rows$marker_set[-1] != marker_rows$marker_set[-NROW(marker_rows)])
        grDevices::png(heatmap_path, width = width * 160, height = height * 160, res = 160)
        pheatmap::pheatmap(scaled, color = grDevices::colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(101), breaks = seq(-2.5, 2.5, length.out = 102), cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = annotation_row, annotation_colors = list(`Marker set` = marker_colors), gaps_row = group_breaks, border_color = NA, fontsize = 11, main = paste0(panel_label, " — grouped by ", tolower(group_label)))
        grDevices::dev.off()
      } else {
        heatmap_data <- data.frame(gene = rep(rownames(scaled), times = NCOL(scaled)), group = rep(colnames(scaled), each = NROW(scaled)), z_score = as.vector(scaled), stringsAsFactors = FALSE)
        heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$group, y = .data$gene, fill = .data$z_score)) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), name = "Row Z-score") + ggplot2::labs(title = panel_label, subtitle = paste0("Mean normalized expression, row-scaled and grouped by ", tolower(group_label)), x = group_label, y = "Marker gene") + ggplot2::theme_classic(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
        save_plot(heatmap, basename(heatmap_path), width, height)
      }
    }
  }
  invisible(TRUE)
}

save_cluster_marker_heatmaps <- function(obj, markers, genes_per_cluster = 10L, clusters_per_panel = 4L) {
  if (!NROW(markers) || "error" %in% names(markers) || !"cluster" %in% names(markers)) return(invisible(FALSE))
  gene_hits <- intersect(c("gene", "features"), names(markers))
  gene_column <- if (length(gene_hits)) gene_hits[[1]] else ""
  if (!nzchar(gene_column)) return(invisible(FALSE))
  top <- do.call(rbind, lapply(split(markers, as.character(markers$cluster)), function(x) {
    if ("p_val_adj" %in% names(x)) x <- x[order(x$p_val_adj, na.last = TRUE), , drop = FALSE]
    utils::head(x, genes_per_cluster)
  }))
  if (!NROW(top)) return(invisible(FALSE))
  cluster_order <- unique(as.character(top$cluster))
  cluster_order <- cluster_order[order(suppressWarnings(as.numeric(cluster_order)), cluster_order, na.last = TRUE)]
  panels <- split(cluster_order, ceiling(seq_along(cluster_order) / clusters_per_panel))
  expression <- Seurat::GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
  all_clusters <- as.character(obj$cluster)
  all_cluster_levels <- unique(all_clusters[order(suppressWarnings(as.numeric(all_clusters)), all_clusters, na.last = TRUE)])
  for (panel_index in seq_along(panels)) {
    panel_clusters <- panels[[panel_index]]
    panel_top <- top[as.character(top$cluster) %in% panel_clusters, , drop = FALSE]
    panel_top <- panel_top[panel_top[[gene_column]] %in% rownames(expression), , drop = FALSE]
    if (!NROW(panel_top)) next
    genes <- as.character(panel_top[[gene_column]])
    row_labels <- make.unique(paste0(genes, "  [cluster ", as.character(panel_top$cluster), "]"))
    means <- do.call(cbind, lapply(all_cluster_levels, function(cluster) Matrix::rowMeans(expression[genes, all_clusters == cluster, drop = FALSE])))
    means <- as.matrix(means); rownames(means) <- row_labels; colnames(means) <- all_cluster_levels
    scaled <- t(scale(t(means))); scaled[!is.finite(scaled)] <- 0
    scaled <- matrix(pmax(-2.5, pmin(2.5, scaled)), nrow = NROW(means), ncol = NCOL(means), dimnames = dimnames(means))
    heatmap_data <- data.frame(gene = rep(rownames(scaled), times = NCOL(scaled)), cluster = rep(colnames(scaled), each = NROW(scaled)), z_score = as.vector(scaled), stringsAsFactors = FALSE)
    plot <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$cluster, y = .data$gene, fill = .data$z_score)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), name = "Row Z-score") +
      ggplot2::labs(title = paste0("Top cluster markers: ", paste(panel_clusters, collapse = ", ")), subtitle = "Top 10 ranked markers per selected cluster; mean normalized expression, row-scaled", x = "Cluster", y = "Marker gene [source cluster]") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
    save_plot(plot, sprintf("08_cluster_marker_heatmap_panel_%02d.png", panel_index), max(8, 3.8 + 0.58 * length(all_cluster_levels)), max(6, 2.8 + 0.22 * NROW(panel_top)))
  }
  invisible(TRUE)
}

load_seurat_reference <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "rds")) {
    reference <- readRDS(path)
    source_name <- basename(path)
  } else if (identical(extension, "rda")) {
    reference_environment <- new.env(parent = emptyenv())
    loaded_names <- load(path, envir = reference_environment)
    seurat_names <- loaded_names[vapply(loaded_names, function(name) inherits(reference_environment[[name]], "Seurat"), logical(1))]
    if (!length(seurat_names)) stop("The .rda file does not contain a Seurat object.")
    if (length(seurat_names) > 1L) message("The reference .rda contains multiple Seurat objects; using: ", seurat_names[[1]])
    source_name <- seurat_names[[1]]
    reference <- reference_environment[[source_name]]
  } else {
    stop("Seurat reference transfer accepts only .rda or .rds files.")
  }
  if (!inherits(reference, "Seurat")) stop("The selected reference is not a Seurat object; detected class: ", paste(class(reference), collapse = ", "))
  list(object = SeuratObject::UpdateSeuratObject(reference), source_name = source_name)
}

ensure_log_normalized_rna <- function(object, label) {
  if (!"RNA" %in% names(object@assays)) stop(label, " does not contain an RNA assay required for reference transfer.")
  # JoinLayers is a Seurat v5 Assay5 operation. Published references such as
  # the Baccin bone-marrow object are valid legacy Seurat Assay objects after
  # UpdateSeuratObject(), but do not implement JoinLayers. Their counts/data
  # slots are already single matrices and must be left intact.
  if (inherits(object[["RNA"]], "Assay5")) {
    rna_layers <- assay_layers_safe(object, "RNA")
    if (any(grepl("^(counts|data)\\.", rna_layers))) {
      object <- SeuratObject::JoinLayers(object, assay = "RNA")
    }
  }
  DefaultAssay(object) <- "RNA"
  layers <- assay_layers_safe(object, "RNA")
  if (!"data" %in% layers) {
    if (!"counts" %in% layers) stop(label, " RNA assay has neither a normalized data layer nor a counts layer.")
    object <- Seurat::NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
  }
  object
}

infer_rna_symbol_species <- function(genes, orthologs) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  mouse_symbols <- unique(orthologs$mouse)
  human_symbols <- unique(orthologs$human)
  # Symbols that are identical in both species provide no directional
  # evidence. Restricting detection to species-specific exact symbols avoids
  # classifying a mixed or Ensembl-ID feature set from capitalization alone.
  mouse_only <- setdiff(mouse_symbols, human_symbols)
  human_only <- setdiff(human_symbols, mouse_symbols)
  mouse_hits <- sum(genes %in% mouse_only)
  human_hits <- sum(genes %in% human_only)
  minimum_hits <- max(10L, min(100L, ceiling(length(genes) * 0.005)))
  species <- if (mouse_hits >= minimum_hits && mouse_hits >= 1.5 * max(1L, human_hits)) {
    "mouse"
  } else if (human_hits >= minimum_hits && human_hits >= 1.5 * max(1L, mouse_hits)) {
    "human"
  } else {
    "unknown"
  }
  list(species = species, mouse_hits = mouse_hits, human_hits = human_hits, minimum_hits = minimum_hits)
}

convert_reference_rna_species <- function(reference, reference_labels, source_species, target_species, orthologs, annotation_name) {
  if (!source_species %in% c("mouse", "human") || !target_species %in% c("mouse", "human") || identical(source_species, target_species)) {
    stop("Reference conversion requires two different recognized species.")
  }
  source_counts <- table(orthologs[[source_species]])
  target_counts <- table(orthologs[[target_species]])
  one_to_one <- orthologs[
    source_counts[orthologs[[source_species]]] == 1L & target_counts[orthologs[[target_species]]] == 1L,
    , drop = FALSE
  ]
  lookup <- stats::setNames(one_to_one[[target_species]], one_to_one[[source_species]])
  original_genes <- rownames(reference)
  mapped_genes <- unname(lookup[original_genes])
  mapped <- !is.na(mapped_genes) & nzchar(mapped_genes)
  mapping_audit <- data.frame(
    original_reference_gene = original_genes,
    converted_reference_gene = ifelse(mapped, mapped_genes, ""),
    status = ifelse(mapped, "mapped_one_to_one", "unmapped_or_ambiguous"),
    reference_species = source_species,
    query_species = target_species,
    stringsAsFactors = FALSE
  )
  mapping_path <- file.path(tables_dir, paste0("reference_transfer_ortholog_mapping__", annotation_name, ".tsv"))
  utils::write.table(mapping_audit, mapping_path, sep = "\t", row.names = FALSE, quote = FALSE)
  if (sum(mapped) < 50L) {
    stop("Only ", sum(mapped), " one-to-one ", source_species, "-to-", target_species, " orthologs were available for the reference; at least 50 are required.")
  }
  counts <- tryCatch(Seurat::GetAssayData(reference, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) {
    stop("Cross-species reference conversion requires the reference RNA counts layer; normalized values alone are not converted as counts.")
  }
  converted_counts <- counts[mapped, , drop = FALSE]
  rownames(converted_counts) <- mapped_genes[mapped]
  converted <- Seurat::CreateSeuratObject(
    counts = converted_counts,
    assay = "RNA",
    project = paste0("reference_", target_species),
    meta.data = reference@meta.data[colnames(converted_counts), , drop = FALSE]
  )
  reference_labels <- reference_labels[colnames(converted)]
  Seurat::Idents(converted) <- factor(reference_labels)
  list(object = converted, labels = reference_labels, mapped_genes = sum(mapped), total_genes = length(original_genes), mapping_path = mapping_path)
}

apply_reference_annotation <- function(obj, path, annotation_name, label_column = "") {
  if (!file.exists(path)) stop("Seurat reference file was not found: ", path)
  loaded <- load_seurat_reference(path)
  reference <- loaded$object
  label_column <- trimws(as.character(label_column %||% ""))
  reference_labels <- if (nzchar(label_column)) {
    if (!label_column %in% colnames(reference@meta.data)) {
      available <- colnames(reference@meta.data)
      stop(
        "Reference label column was not found: ", label_column, ". ",
        "Leave the label selection on Active identities or choose a column reported by Inspect reference. ",
        "Available metadata columns: ", if (length(available)) paste(available, collapse = ", ") else "none"
      )
    }
    as.character(reference[[label_column]][, 1])
  } else {
    as.character(Seurat::Idents(reference))
  }
  names(reference_labels) <- colnames(reference)
  valid_labels <- !is.na(reference_labels) & nzchar(trimws(reference_labels))
  if (!all(valid_labels)) {
    reference <- subset(reference, cells = names(reference_labels)[valid_labels])
    reference_labels <- reference_labels[valid_labels]
  }
  if (ncol(reference) < 20L || length(unique(reference_labels)) < 2L) stop("Reference transfer needs at least 20 labeled cells spanning at least two labels.")
  reference <- ensure_log_normalized_rna(reference, "Reference")
  obj <- ensure_log_normalized_rna(obj, "Query")

  ortholog_path <- trimws(as.character(params$reference_ortholog_file %||% ""))
  if (!nzchar(ortholog_path)) {
    script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else ""
    if (nzchar(script_file)) {
      ortholog_path <- file.path(dirname(dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))), "reference", "mouse_human_orthologs_MGI.tsv")
    }
  }
  if (!nzchar(ortholog_path) || !file.exists(ortholog_path)) {
    stop("Reference transfer requires the bundled readable mouse-human ortholog table: ", ortholog_path)
  }
  orthologs <- read_ortholog_symbols(ortholog_path)
  reference_species_info <- infer_rna_symbol_species(rownames(reference), orthologs)
  query_species_info <- infer_rna_symbol_species(rownames(obj), orthologs)
  reference_species <- reference_species_info$species
  query_species <- query_species_info$species
  shared_before_conversion <- length(intersect(rownames(reference), rownames(obj)))
  converted_reference <- FALSE
  converted_genes <- 0L
  reference_genes_before_conversion <- nrow(reference)
  force_reference_pca <- FALSE
  if (reference_species %in% c("mouse", "human") && query_species %in% c("mouse", "human") && !identical(reference_species, query_species)) {
    message("Detected ", reference_species, " reference and ", query_species, " query; converting the reference with one-to-one MGI orthologs.")
    conversion <- convert_reference_rna_species(reference, reference_labels, reference_species, query_species, orthologs, annotation_name)
    reference <- conversion$object
    reference_labels <- conversion$labels
    converted_reference <- TRUE
    converted_genes <- conversion$mapped_genes
    force_reference_pca <- TRUE
    reference <- ensure_log_normalized_rna(reference, "Converted reference")
  } else if (identical(reference_species, "unknown") || identical(query_species, "unknown")) {
    if (shared_before_conversion < 50L) {
      stop(
        "Could not confidently detect reference/query species from gene symbols, and they share only ", shared_before_conversion, " exact RNA features. ",
        "Use mouse or human gene symbols (not mixed identifiers) so automatic ortholog conversion can be applied."
      )
    }
    message("Species detection was inconclusive, but reference and query share ", shared_before_conversion, " exact RNA features; proceeding without conversion.")
  } else {
    message("Detected matching ", reference_species, " reference and query; no ortholog conversion is needed.")
  }

  reduction_name <- "pca"
  usable_pca <- !force_reference_pca && reduction_name %in% names(reference@reductions) &&
    ncol(Seurat::Embeddings(reference, reduction = reduction_name)) >= 5L &&
    identical(reference[[reduction_name]]@assay.used, "RNA")
  if (!usable_pca) {
    reference <- Seurat::FindVariableFeatures(reference, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    variable_genes <- intersect(Seurat::VariableFeatures(reference, assay = "RNA"), rownames(reference))
    if (length(variable_genes) < 50L) stop("The reference has too few variable RNA genes for label transfer.")
    reference <- Seurat::ScaleData(reference, assay = "RNA", features = variable_genes, verbose = FALSE)
    reference <- Seurat::RunPCA(reference, assay = "RNA", features = variable_genes, npcs = min(50L, params$n_pcs), reduction.name = "codespring_reference_pca", verbose = FALSE)
    reduction_name <- "codespring_reference_pca"
  }
  transfer_features <- intersect(Seurat::VariableFeatures(reference, assay = "RNA"), rownames(obj))
  if (length(transfer_features) < 50L) {
    transfer_features <- intersect(rownames(reference), rownames(obj))
    transfer_features <- utils::head(transfer_features, 3000L)
  }
  if (length(transfer_features) < 50L) stop("Reference and query share too few RNA genes for reliable label transfer (", length(transfer_features), " shared).")
  available_pcs <- ncol(Seurat::Embeddings(reference, reduction = reduction_name))
  dims <- seq_len(min(30L, params$n_pcs, available_pcs))
  if (length(dims) < 5L) stop("Reference PCA has fewer than five usable dimensions.")
  transfer_started <- Sys.time()
  message(
    "Phase 1/3: finding reference-query anchors for ", ncol(obj), " query cells using ",
    length(transfer_features), " shared RNA features and ", length(dims), " PCs."
  )
  flush.console()
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = obj,
    normalization.method = "LogNormalize",
    reference.assay = "RNA",
    query.assay = "RNA",
    reference.reduction = reduction_name,
    features = transfer_features,
    dims = dims,
    verbose = TRUE
  )
  message("Phase 1/3 complete after ", round(as.numeric(difftime(Sys.time(), transfer_started, units = "mins")), 1), " minutes.")
  message("Phase 2/3: transferring ", length(unique(reference_labels)), " reference labels to all query cells.")
  flush.console()
  prediction_started <- Sys.time()
  predictions <- Seurat::TransferData(anchorset = anchors, refdata = reference_labels, dims = dims, verbose = TRUE)
  message("Phase 2/3 complete after ", round(as.numeric(difftime(Sys.time(), prediction_started, units = "mins")), 1), " minutes; preparing metadata outputs.")
  flush.console()
  labels <- as.character(predictions$predicted.id)
  scores <- as.numeric(predictions$prediction.score.max)
  names(labels) <- rownames(predictions); names(scores) <- rownames(predictions)
  labels <- labels[colnames(obj)]; scores <- scores[colnames(obj)]
  if (any(is.na(labels) | !nzchar(labels))) stop("Reference transfer did not return a label for every query cell.")
  score_column <- paste0(annotation_name, "_prediction_score")
  obj[[annotation_name]] <- unname(labels)
  obj[[score_column]] <- unname(scores)
  obj[[paste0("annotation_source__", annotation_name)]] <- "Seurat RNA anchor label transfer"
  obj$annotation_source <- paste0("Seurat RNA anchor label transfer (", annotation_name, ")")
  per_cell <- data.frame(cell = colnames(obj), transferred_label = unname(labels), prediction_score_max = unname(scores), stringsAsFactors = FALSE)
  names(per_cell)[[2]] <- annotation_name
  utils::write.table(per_cell, file.path(tables_dir, paste0("reference_transfer_per_cell__", annotation_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  score_groups <- split(unname(scores), unname(labels))
  label_summary <- do.call(rbind, lapply(names(score_groups), function(label) {
    values <- score_groups[[label]]
    data.frame(label = label, cells = length(values), median_score = stats::median(values), mean_score = mean(values), stringsAsFactors = FALSE)
  }))
  utils::write.table(label_summary, file.path(tables_dir, paste0("reference_transfer_label_summary__", annotation_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  audit <- data.frame(
    reference_file = normalizePath(path, winslash = "/", mustWork = TRUE),
    reference_object = loaded$source_name,
    label_source = if (nzchar(label_column)) label_column else "active identities",
    reference_species_detected = reference_species,
    query_species_detected = query_species,
    reference_mouse_symbol_hits = reference_species_info$mouse_hits,
    reference_human_symbol_hits = reference_species_info$human_hits,
    query_mouse_symbol_hits = query_species_info$mouse_hits,
    query_human_symbol_hits = query_species_info$human_hits,
    reference_converted_to_query_species = converted_reference,
    ortholog_method = if (converted_reference) "MGI one-to-one reference RNA count conversion" else "none",
    reference_genes_before_conversion = reference_genes_before_conversion,
    reference_genes_after_conversion = nrow(reference),
    one_to_one_reference_genes_mapped = converted_genes,
    exact_shared_features_before_conversion = shared_before_conversion,
    reference_cells = ncol(reference),
    reference_labels = length(unique(reference_labels)),
    query_cells = ncol(obj),
    shared_features_used = length(transfer_features),
    pcs_used = length(dims),
    stringsAsFactors = FALSE
  )
  utils::write.table(audit, file.path(tables_dir, paste0("reference_transfer_audit__", annotation_name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  message("Phase 3/3: transferred-label metadata and audit tables are ready.")
  flush.console()
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
if (!inherits(obj, "Seurat")) stop("Annotation input is not a Seurat object; detected class: ", paste(class(obj), collapse = ", "))
annotation_name <- safe_metadata_name(params$annotation_name)
reference_annotation_only <- FALSE
# Rejoin Seurat v5 RNA layers before any annotation scoring or differential
# expression. Values are unchanged; the final object becomes portable and its
# raw counts remain directly accessible.
rna_layers <- assay_layers_safe(obj, "RNA")
if (any(grepl("^(counts|data)\\.", rna_layers))) obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
message("Annotation input loaded: ", ncol(obj), " cells; reductions: ", paste(names(obj@reductions), collapse = ", "))
if (nzchar(params$reference_file) && !identical(tolower(params$reference_file), "none")) {
  obj <- apply_reference_annotation(obj, params$reference_file, annotation_name, params$reference_label_column)
  reference_annotation_only <- identical(stage, "annotate")
} else if (nzchar(params$celltype_file) && !identical(tolower(params$celltype_file), "none")) {
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

if (isTRUE(reference_annotation_only)) {
  # Reference transfer is intentionally metadata-only. The transfer helper
  # already wrote the per-cell labels/scores and audit tables, so avoid
  # regenerating UMAPs, marker tests, dashboards, and composition summaries.
  attr(obj, "codespring_active_annotation") <- annotation_name
  attr(obj, "codespring_integration") <- integration
  attr(obj, "codespring_doublets_removed") <- sum(doublet_summary$removed_doublets)
  # Reference transfer does not alter coordinates, but its new labels must be
  # republished beside those coordinates for the interactive UMAP explorer.
  write_interactive_metadata_tables(obj)
  message("Saving the updated Seurat object; this final write can take several minutes for a large object.")
  flush.console()
  atomic_save_rds(obj, file.path(objects_dir, "processed_seurat.rds"))
  summary_lines <- c(
    paste("engine: seurat"), paste("normalization:", params$normalization), paste("integration:", integration),
    paste("doublet_method:", params$doublet_method), paste("doublets_removed:", sum(doublet_summary$removed_doublets)),
    paste("input_processing_inventory:", file.path("tables", "input_processing_detected.tsv")),
    paste("input_samples:", NROW(samples)), paste("cells_after_qc:", ncol(obj)), paste("clusters:", length(unique(obj$cluster))), paste("active_annotation:", annotation_name),
    paste("annotation_source:", unique(obj$annotation_source)[1]), paste("annotation_output: metadata only"), paste("generated:", as.character(Sys.time()))
  )
  writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
  stage_marker("annotate")
  writeLines(as.character(Sys.time()), file.path(out_dir, "_COMPLETE"))
  message("Reference annotation metadata saved; no UMAP, marker, dashboard, or composition outputs were regenerated.")
  quit(save = "no", status = 0L)
}

message("Annotation labels prepared; rendering UMAP summaries.")
save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "sample_id", shuffle = TRUE), "03_umap_sample.png", 8, 6)
message("Rendered sample UMAP.")
save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "cluster", label = TRUE, repel = TRUE), "04_umap_clusters.png", 8, 6)
message("Rendered cluster UMAP.")
save_plot(categorical_dimplot(obj, reduction = "umap", group.by = annotation_name, label = TRUE, repel = TRUE), paste0("05_umap_", annotation_name, ".png"), 9, 6)
if ("condition" %in% colnames(obj@meta.data) && length(unique(obj$condition)) > 1L) save_plot(categorical_dimplot(obj, reduction = "umap", group.by = "condition", shuffle = TRUE), "06_umap_condition.png", 8, 6)
marker_sets <- attr(obj, "codespring_marker_gene_sets")
if (length(marker_sets)) {
  message("Rendering readable marker-list dot plots and heatmaps.")
  save_marker_annotation_panels(obj, marker_sets, annotation_name)
}

# Cluster markers should use the normalized SCT representation when that is
# the selected workflow.
DefaultAssay(obj) <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
if (isTRUE(params$find_cluster_markers) && ncol(obj) >= 20 && length(unique(obj$cluster)) > 1L) {
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
    save_cluster_marker_heatmaps(obj, markers)
  }
} else if (!isTRUE(params$find_cluster_markers)) {
  message("Skipping optional cluster-marker discovery; reference annotation and its UMAP/table outputs are complete without it.")
}

# Match the Scanpy output schema: retain an input metadata column named
# `cell` as `input_cell`, while always reserving `cell` for the exact final
# barcode/cell identifier.  This prevents duplicate column names in tables
# consumed by the interactive results explorer.
write_interactive_metadata_tables(obj)
# The interactive dashboard requests one normalized gene at a time from the
# processed RDS. Keep the complete symbol list small and engine-consistent.
utils::write.table(data.frame(gene = rownames(obj)), file.path(tables_dir, "dashboard_all_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
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
atomic_save_rds(obj, file.path(objects_dir, "processed_seurat.rds"))
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
  expanded_signatures <- do.call(rbind, lapply(seq_len(NROW(signatures)), function(i) {
    genes <- unique(Filter(nzchar, strsplit(as.character(signatures$gene[[i]]), "[,;[:space:]]+")[[1]]))
    data.frame(signature = as.character(signatures$signature[[i]]), gene = genes, stringsAsFactors = FALSE)
  }))
  if (!identical(params$signature_species, "same")) {
    map_cross_species_genes(expanded_signatures$gene, params$signature_species, params$signature_ortholog_file, rownames(obj), file.path(tables_dir, "signature_ortholog_mapping.tsv"), expanded_signatures$signature)
    audit <- read_delim_safe(file.path(tables_dir, "signature_ortholog_mapping.tsv"))
    expanded_signatures <- data.frame(signature = audit$set[audit$status %in% c("mapped", "same_species")], gene = audit$mapped_gene[audit$status %in% c("mapped", "same_species")], stringsAsFactors = FALSE)
  }
  signatures <- expanded_signatures
  annotation_name <- safe_metadata_name(attr(obj, "codespring_active_annotation") %||% params$annotation_name)
  # Signature scoring is valid before formal annotation. In that case retain
  # the stable cluster labels for its per-cell and summary outputs.
  if (!annotation_name %in% colnames(obj@meta.data)) obj[[annotation_name]] <- obj$cluster
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
  atomic_save_rds(obj, file.path(objects_dir, "processed_seurat.rds"))
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
