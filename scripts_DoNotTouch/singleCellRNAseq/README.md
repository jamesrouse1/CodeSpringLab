# CodeSpringLab scRNA-seq workflow

The CodeSpringApp single-cell section runs this workflow as one SLURM job. It
keeps source inputs read-only and writes every result below the project's
`data/scrna/` directory.

## Inputs

In CodeSpringApp, a manifest is optional: start by selecting one raw-count
input and sample ID, then add or edit further rows in the project-local
manifest. You can instead provide a tab-delimited manifest with `sample_id`
and `input_path` columns when it already exists. Each row can point to one of
the following raw-count inputs:

- a Seurat `.rds` object (processed with the Seurat engine);
- an AnnData `.h5ad` object (processed with the Scanpy engine); or
- a filtered 10x matrix directory (processed with the engine selected in the app).

Additional columns such as `condition`, `batch`, and `donor` become cell
metadata. `sample_id` must be unique. Do not mix Seurat `.rds` and AnnData
`.h5ad` inputs in a single run; start from filtered 10x matrices instead.

## Runtime on the HPC

Starting `./run_codespringweb.sh` only starts the lightweight web interface;
it does **not** load Scanpy or Seurat into the app process. Each scRNA run is
submitted as one SLURM job and the selected engine is passed explicitly to its
job wrapper. The wrapper loads the standard cluster R module for Seurat or the
standard Anaconda module for Scanpy only inside that compute job.

If your HPC Scanpy or Seurat installation is in a different environment,
enter its executable in **Optional custom runtime executable** in the app:

- Seurat: the full path to the appropriate `Rscript` executable.
- Scanpy: the full path to the environment's `python` executable.

The app verifies that path before submission and records it with the run
parameters. The job checks essential packages before loading a large object,
so a missing runtime is reported promptly in the job log.

The Scanpy engine requires raw counts in `X`, `layers['counts']`/`layers['raw']`, or `raw`.
It stops with an explanatory error when only normalized expression is present,
because repeating QC and normalization on processed values is not valid.

For an AnnData `.h5ad`, CodeSpring records whether raw counts came from `X`, a
`counts` layer, or `.raw`, along with existing embeddings, cluster fields, and
annotation metadata. This is written to `tables/input_processing_detected.tsv`.
If a UMAP is present it is rendered as an input-reference figure before a new
embedding is computed. Existing cell-type metadata is retained unless a new
marker list or cell-level mapping is supplied.

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
  PCA variance, metadata, marker, annotation, exact cell-type composition by
  sample, and `umap_coordinates.tsv` tables;
- `objects/`: processed Seurat RDS or AnnData H5AD;
- `run_summary.txt`: processing choices and cell/cluster counts;
- `_COMPLETE`: written only after a successful workflow.

The project output location is selected as **Results root** during project
creation. CodeSpring writes this project beneath that folder as
`<results_root>/<project_name>/data/scrna/`, keeping the processed object,
tables, figures, logs, and temporary job storage together without modifying
the source matrix or source object.

`tables/umap_coordinates.tsv` is intentionally a small, engine-neutral view
of the final embedding plus cell metadata. CodeSpringApp uses it for the
interactive **Explore Cells** view, including metadata colouring and
lasso/box selection, without reopening the full RDS or H5AD in the web
process. The processed object and complete metadata table remain available
for download.

`tables/cell_type_by_sample.tsv` provides exact per-sample cell-type counts
and within-sample proportions. It powers the Results Explorer composition
view and is intentionally calculated before the browser-level UMAP sampling.
