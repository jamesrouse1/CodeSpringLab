# CodeSpringLab scRNA-seq workflow

The CodeSpringApp single-cell section runs this workflow as one SLURM job. It
keeps source inputs read-only and writes every result below the project's
`data/scrna/` directory.

## Inputs

Provide a tab-delimited manifest with `sample_id` and `input_path` columns.
Each row can point to one of the following raw-count inputs:

- a Seurat `.rds` object (processed with the Seurat engine);
- an AnnData `.h5ad` object (processed with the Scanpy engine); or
- a filtered 10x matrix directory (processed with the engine selected in the app).

Additional columns such as `condition`, `batch`, and `donor` become cell
metadata. `sample_id` must be unique. Do not mix Seurat `.rds` and AnnData
`.h5ad` inputs in a single run; start from filtered 10x matrices instead.

The Scanpy engine requires raw counts in `X`, `layers['counts']`, or `raw`.
It stops with an explanatory error when only normalized expression is present,
because repeating QC and normalization on processed values is not valid.

For a Seurat RDS, CodeSpring inspects the RNA count layers, SCT/integrated
assays, reductions, cluster metadata, and existing annotation columns. The
detected state is saved in `tables/input_processing_detected.tsv`. Raw RNA
count layers are used for a reproducible rerun of the selected QC and
downstream workflow; if the input already contains a UMAP, it is also rendered
as an input-reference figure before the new embedding is made. Existing
cell-type metadata is retained unless a new marker list or cell-level mapping
is supplied.

## Processing

The workflow performs initial per-sample QC, optional doublet detection,
gene filtering, normalization, highly-variable-gene selection, scaling,
integration when a genuine technical batch column is supplied, PCA,
neighbors/UMAP, clustering, annotation, and cluster markers. It retains raw
counts in the processed output. Seurat defaults to SCTransform v2 and Scanpy
to log-normalization. In automatic mode, integration is used only when the
selected batch column contains more than one value; sample count alone does
not trigger batch correction.

Doublets are detected after initial cell QC and before normalization or
integration. Scanpy uses per-sample Scrublet and requires the Scanpy Scrublet
dependencies (including `scikit-image`). Seurat uses `scDblFinder`, which must
be installed in the R environment (for example, with
`BiocManager::install("scDblFinder")`). Every run writes a per-cell call table,
a per-sample summary, and a score distribution when scores are available.

For annotation, use either:

- a marker file with `cell_type` and `gene` columns; or
- a cell-to-cell-type mapping with `cell` (or `barcode`) and `cell_type`
  columns.

A cell mapping takes precedence. Marker lists are scored per cluster and the
best-scoring label is assigned to each cluster.

## Output layout

- `figures/`: QC and UMAP figures;
- `tables/`: per-cell QC, per-sample QC, doublet calls, highly-variable genes,
  PCA variance, metadata, marker, and annotation tables;
- `objects/`: processed Seurat RDS or AnnData H5AD;
- `run_summary.txt`: processing choices and cell/cluster counts;
- `_COMPLETE`: written only after a successful workflow.
