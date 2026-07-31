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

## Processing

The workflow writes raw-count-preserving processed objects plus QC metrics,
UMAPs, clusters, annotations, and cluster-marker tables. Seurat defaults to
SCTransform v2 with RPCA integration for multiple samples. Scanpy defaults to
log-normalization with scVI integration for multiple samples. Both defaults can
be changed in CodeSpringApp when the study design justifies it.

For annotation, use either:

- a marker file with `cell_type` and `gene` columns; or
- a cell-to-cell-type mapping with `cell` (or `barcode`) and `cell_type`
  columns.

A cell mapping takes precedence. Marker lists are scored per cluster and the
best-scoring label is assigned to each cluster.

## Output layout

- `figures/`: QC and UMAP figures;
- `tables/`: per-cell QC, per-sample QC, metadata, marker, and annotation
  tables;
- `objects/`: processed Seurat RDS or AnnData H5AD;
- `run_summary.txt`: processing choices and cell/cluster counts;
- `_COMPLETE`: written only after a successful workflow.
