#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

while IFS= read -r script; do
  bash -n "$script"
done < <(find "$repo_root/scripts_DoNotTouch" "$repo_root/tests" -type f -name '*.sh' -print)

python3 - "$repo_root" <<'PY'
import ast
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
for path in sorted((root / "scripts_DoNotTouch").rglob("*.py")):
    ast.parse(path.read_text(), filename=str(path))
PY

Rscript -e 'root <- commandArgs(TRUE)[1]; files <- list.files(file.path(root, "scripts_DoNotTouch"), pattern="[.][Rr]$", recursive=TRUE, full.names=TRUE); for (file in files) parse(file=file)' "$repo_root"
Rscript "$repo_root/tests/smoke_test_completed_rnaseq_comparisons.R" "$repo_root"

grep -Fq 'pca_differential_peaks.png' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
grep -Fq 'dba.plotPCA(db, contrast = 1L' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
grep -Fq 'all_differential_peaks.bed' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
grep -Fq 'completed_deseq_comparisons <- function()' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'completed_gsea_comparisons <- function()' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'choices = completed_deseq_columns' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'choices = completed_gsea_columns' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'completed_controls(completed_deseq_catalog' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'completed_controls(completed_gsea_catalog' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'configured_species, "maize"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'maize_nc350_nam1' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'maize_w22_nrgene2' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'optional_quantifiers_enabled' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'uiOutput("rna_overview_sample_progress_ui")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'div(class = "metric-card tone-blue", span("Disk space")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
if grep -Fq 'table_widget("rna_overview_design")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"; then
  echo "RNA-seq Overview must not display the design matrix" >&2
  exit 1
fi
if grep -Fq 'table_widget("rna_overview_status")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"; then
  echo "RNA-seq Overview must use the live sample progress matrix instead of the old pipeline-status table" >&2
  exit 1
fi
grep -Fq 'actionButton("refresh_rna_results", "Refresh results"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq '"log2-transformed DESeq2 normalized counts"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'No DESeq2 normalized-counts file is available for PCA. Finish DESeq2, then click Refresh results.' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
if grep -Fq '"featureCounts count matrix" = "counts"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"; then
  echo "RNA-seq PCA must not offer raw featureCounts values" >&2
  exit 1
fi
grep -Fq 'selectInput("rna_file_tool", "Tool"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'uiOutput("rna_file_sample_ui")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'File = basename(absolute)' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
grep -Fq 'c("Tool", "Sample", "File", "Size", "Modified", "Copy path")' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"
if grep -Fq 'selectInput("rna_file_category"' "$repo_root/scripts_DoNotTouch/Shiny/app_server.R"; then
  echo "RNA file browsing must use the tool filter rather than the old category filter" >&2
  exit 1
fi

# Automatic Seurat integration must remain on the scalable PCA/Harmony path.
# Anchor integration can exceed Matrix's 32-bit sparse-index ceiling on large
# multi-sample datasets even when sufficient RAM is available.
test "$(grep -Fc 'if (identical(integration, "auto")) integration <- if (length(unique(batch_values[nzchar(batch_values)])) > 1L) "harmony" else "none"' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R")" -eq 2
grep -Fq 'batch_column = get("batch_column", "sample_id")' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R"
grep -Fq 'reference_file = get("reference_file", "")' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R"
grep -Fq 'FindTransferAnchors(' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R"
grep -Fq 'reference_transfer_per_cell__' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R"
grep -Fq 'run_umap_with_allocated_threads <- function' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R"
test "$(grep -Fc 'run_umap_with_allocated_threads(' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R")" -eq 2
grep -Fq '#SBATCH --cpus-per-task=8' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/qsub_scrna_pipeline.sh"
grep -Fq '#SBATCH --mem=64G' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/qsub_scrna_pipeline.sh"
grep -Fq 'export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/qsub_scrna_pipeline.sh"
grep -Fq '#SBATCH --cpus-per-task=20' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/qsub_cellranger_count.sh"
grep -Fq 'export CELLRANGER_LOCALMEM_GB=$((allocated_mem_gb - 6))' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/qsub_cellranger_count.sh"
test "$(grep -Fc 'write_interactive_metadata_tables(obj)' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_seurat.R")" -ge 2
grep -Fq '"batch_column": get("batch_column", "sample_id")' "$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline_scanpy.py"

python3 - "$repo_root" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
manifest = root / "scripts_DoNotTouch/test/manifest_atac/design_matrix.txt"
fastq_dir = root / "scripts_DoNotTouch/test/fastq_atac"
with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

required = {"sample", "condition", "replicate", "filename"}
if not rows or not required.issubset(rows[0]):
    raise SystemExit("ATAC example manifest is missing required sample metadata columns")

samples = [row["sample"].strip() for row in rows]
if any(not sample for sample in samples) or len(samples) != len(set(samples)):
    raise SystemExit("ATAC example sample identifiers must be non-empty and unique")

condition_counts = {}
for row in rows:
    condition = row["condition"].strip()
    replicate = row["replicate"].strip()
    if not condition or not replicate:
        raise SystemExit("ATAC example condition and replicate values must be non-empty")
    condition_counts[condition] = condition_counts.get(condition, 0) + 1
    files = [item.strip() for item in row["filename"].split(",") if item.strip()]
    if len(files) != 2 or not all((fastq_dir / item).is_file() for item in files):
        raise SystemExit(f"ATAC example FASTQ pair does not resolve for {row['sample']}")

if len(condition_counts) < 2 or any(count < 2 for count in condition_counts.values()):
    raise SystemExit("ATAC example needs at least two replicates in each comparison condition")
PY

python3 - "$repo_root" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
manifest = root / "scripts_DoNotTouch/test/manifest_chip/design_matrix.txt"
fastq_dir = root / "scripts_DoNotTouch/test/fastq_chip"
with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

targets = [row for row in rows if row["reference"].strip().lower() == "chip"]
inputs = [row for row in rows if row["reference"].strip().lower() == "input"]
if len(targets) != 4 or len(inputs) != 1:
    raise SystemExit("ChIP example must contain four ChIP libraries and one shared input library")
shared_input = inputs[0]["sample"].strip()
if not shared_input or any(row["control_sample"].strip() != shared_input for row in targets):
    raise SystemExit("Every ChIP example target must reference the one shared input library")
if any(not (fastq_dir / row["filename"].strip()).is_file() for row in rows):
    raise SystemExit("A ChIP example FASTQ named in the design matrix is missing")
PY

echo "Shell, Python, and R source syntax checks passed."
