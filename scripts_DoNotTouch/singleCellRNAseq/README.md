# CodeSpringLab scRNA-seq workflow

The CodeSpringApp single-cell section runs this workflow as checkpointed SLURM
stages. It keeps source inputs read-only and writes every result below the project's
`data/scrna/` directory.

## Inputs

In CodeSpringApp, one input needs no manifest: select one raw-count input and
CodeSpring derives its sample ID from the file or folder name. Choose the
multiple-input/integration setup only when combining inputs from different
locations; it reveals a project-local manifest. A multi-input manifest must
have `sample_id` and `input_path` columns. Each row can point to one of
the following raw-count inputs:

- a Seurat `.rds` object (processed with the Seurat engine);
- an AnnData `.h5ad` object (processed with the Scanpy engine); or
- a filtered 10x matrix directory; or
- a folder of gzipped 10x gene-expression R1/R2 FASTQs. FASTQ samples first run
  independently through the cluster's `CellRanger/9.0.1` module and the resulting
  `outs/filtered_feature_bc_matrix` becomes the downstream input.

Additional columns such as `condition`, `batch`, and `donor` become cell
metadata. `sample_id` must be unique. `capture_id` identifies the independent
droplet capture/library used for doublet detection and defaults to `sample_id`.
When multiple hashed or multiplexed biological samples came from the same 10x
channel, give those rows the same `capture_id`. An optional
`expected_doublet_rate` can override automatic rate estimation and must be
constant within a capture. Do not mix Seurat `.rds` and AnnData
`.h5ad` inputs in a single run; start from filtered 10x matrices instead.
For FASTQ inputs, the optional `fastq_sample` column records the filename prefix
passed to `cellranger count --sample`; the app infers it when one prefix is found.

## Runtime on the HPC

Starting `./run_codespringweb.sh` only starts the lightweight web interface;
it does **not** load Scanpy or Seurat into the app process. Each scRNA run is
submitted as a SLURM job and the selected engine is passed explicitly to its
job wrapper. The wrapper loads the versioned cluster Seurat module for Seurat
jobs or the shared versioned container for Scanpy only inside that compute
job. Seurat jobs ignore personal R startup files and user libraries so older
packages cannot override the compatible packages in the tested module.

The job checks essential packages before loading a large object, so a missing
cluster dependency is reported promptly in the job log.

Cell Ranger alignment/counting jobs require a matching `refdata-gex-*`
transcriptome folder. CodeSpringApp automatically looks for the current human
and mouse references below
`/grid/bsr/data/data/bsr_readable_data/references/cellranger` and passes the
selected path to the maintained runner; the runner does not hardcode a species.
Species and transcriptome version are therefore requested only for FASTQ-backed
projects. Processed objects and filtered matrices do not receive a genome-build
setting because alignment has already occurred.
They
run with intronic counting enabled, chemistry auto-detection, BAM creation, and
secondary clustering disabled because CodeSpring performs QC, normalization,
integration, clustering, and annotation in the later checkpointed stages.

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

For FASTQ-backed projects, **Alignment & counting** appears before the standard
stages. The app then exposes **Input inspection**, **QC & doublets**,
**Normalize & PCA**, **Integrate & cluster**, and **Annotate & markers**. Each
stage writes an explicit completion marker plus a checkpoint under
`checkpoints/`; a later stage will not run until its prerequisite checkpoint
exists. Input inspection records pre-existing normalization, reductions,
clusters, and annotations in `tables/input_processing_detected.tsv`, while
the downstream reproducible workflow starts from raw counts. A single RDS or
H5AD that already contains both a UMAP and clusters is also copied into the
project as a reusable continuation object; the source file remains read-only,
and optional reconstruction controls remain available in the app.

The workflow performs initial per-sample QC, optional doublet detection,
gene filtering, normalization, highly-variable-gene selection,
integration when a genuine technical batch column is supplied, PCA,
neighbors/UMAP, clustering, annotation, and cluster markers. It retains raw
counts in the processed output. Seurat defaults to SCTransform v2 and Scanpy
to log-normalization. In automatic mode, Seurat uses RPCA and Scanpy uses
Harmony only when the
selected batch column contains more than one value; sample count alone does
not trigger batch correction.

Scanpy keeps raw counts in `layers["counts"]`, uses count-depth normalization
followed by `log1p`, selects batch-aware highly variable genes from the
log-normalized expression using Scanpy's Seurat-style method, and performs implicitly centered
sparse PCA on those genes. Seurat applies SCTransform v2 separately to each
selected integration unit (or LogNormalize when requested). For anchor-based
RPCA/CCA, inputs sharing the same selected technical-batch value are grouped
before normalization and integration; Harmony instead corrects the selected
metadata column in PCA space.

The normalization/PCA stage always writes an uncorrected UMAP by sample and,
when available, by the selected technical batch. Review that view beside the
final corrected UMAP: a correction should reduce technical separation while
preserving condition- and cell-type-associated structure. Seurat offers RPCA,
CCA, and Harmony; the standard Scanpy container offers Harmony. Harmony theta, lambda, and
iteration settings are recorded in the project parameter table. RPCA is the
conservative automatic Seurat default; CCA is reserved for datasets with
strong shared-state shifts where a stronger correction is justified.

Doublets are detected after only a minimal 200-count prefilter and before the
complete mitochondrial/count/feature QC, normalization, or integration. They
are modeled independently per `capture_id`, not blindly per biological sample.
Scanpy uses capture-aware Scrublet and requires the Scanpy Scrublet
dependencies (including `scikit-image`). Seurat uses `scDblFinder`, which must
be installed in the R environment (for example, with
`BiocManager::install("scDblFinder")`). Doublet rates are estimated per capture
by default; a global manual override remains available. Every run writes a
per-cell call table, a per-capture summary, and a score distribution when scores
are available. Full QC is then applied independently to each sample using the
user-reviewed thresholds.

Mitochondrial percentages remain uncapped in per-cell tables and during
filtering. The mitochondrial axes in QC figures alone use the 99th percentile,
with a 50% maximum display limit, so extreme droplets do not compress the
visible distribution.

Input inspection writes editable QC starting values for each sample. Lower and
upper outliers are estimated with three median absolute deviations (3 MAD),
using log-transformed total counts and detected genes and the original scale for
mitochondrial percentage. The shared app defaults use the most permissive
sample-specific values to reduce the risk of discarding an entire lower-RNA or
higher-mitochondrial population. These remain suggestions: the unfiltered plots,
sample biology, and retained-cell counts should be reviewed before QC is run.

For annotation, use either:

- a marker file with `cell_type` and `gene` columns; or
- a cell-to-cell-type mapping with `cell` (or `barcode`) and `cell_type`
  columns.

A cell mapping takes precedence. Marker lists are scored per cluster and the
best-scoring label is assigned to each cluster.

Marker lists and signature lists can be declared as mouse, human, or already
matching the expression dataset. When conversion is needed, both engines use
the same MGI-style ortholog table accepted by the RNA-seq workflow (for example
`mouse_gene_symbol` and `human_gene_symbol`; common MGI/HGNC aliases are also
recognized). A user-supplied TSV/CSV can replace the bundled table. Ambiguous
source-to-many mappings are excluded, and `marker_ortholog_mapping.tsv` or
`signature_ortholog_mapping.tsv` records every mapped, unchanged, ambiguous,
and unmapped input gene.
Ranked pathway analysis similarly offers automatic species detection or an
explicit mouse/human choice and accepts the same optional ortholog table before
testing human-symbol pathway collections.

## Output layout

- `figures/`: QC, pre-integration UMAP, and final UMAP figures;
- `tables/`: per-cell QC, per-sample QC, doublet calls, highly-variable genes,
  PCA variance, metadata, marker, annotation, exact cell-type composition by
  sample, `preintegration_umap_coordinates.tsv`, and `umap_coordinates.tsv`;
- `objects/`: processed Seurat RDS or AnnData H5AD;
- `checkpoints/`: internal stage checkpoints used to resume the workflow;
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
