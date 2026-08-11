#!/usr/bin/env bash
set -Eeuo pipefail

sample="${1:-}"
fastq_dir="${2:-}"
fastq_sample="${3:-}"
transcriptome="${4:-}"
output_root="${5:-}"
expect_cells="${6:-0}"
stage_parent="${7:-$output_root}"

if [[ -z "$sample" || -z "$fastq_dir" || -z "$transcriptome" || -z "$output_root" ]]; then
  echo "ERROR: usage: cellranger_count.sh <sample> <FASTQ folder> <FASTQ sample prefix or blank> <transcriptome> <output root> [expected cells] [staging parent]" >&2
  exit 2
fi

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    if [[ -r "$module_init" ]]; then
      set +u
      # shellcheck disable=SC1090
      source "$module_init"
      set -u
      break
    fi
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load CellRanger/9.0.1
command -v cellranger >/dev/null 2>&1 || { echo "ERROR: cellranger was not found after loading CellRanger/9.0.1." >&2; exit 127; }

[[ -d "$fastq_dir" && -r "$fastq_dir" ]] || { echo "ERROR: FASTQ folder is missing or unreadable: $fastq_dir" >&2; exit 2; }
[[ -d "$transcriptome" && -s "$transcriptome/reference.json" ]] || { echo "ERROR: Cell Ranger transcriptome is missing reference.json: $transcriptome" >&2; exit 2; }
if ! find "$fastq_dir" -maxdepth 1 -type f \( -iname '*.fastq.gz' -o -iname '*.fq.gz' \) -print -quit | grep -q .; then
  echo "ERROR: no gzipped FASTQs were found in $fastq_dir" >&2
  exit 2
fi

mkdir -p "$output_root" "$stage_parent"
# A stable per-sample staging directory lets Cell Ranger resume an incomplete
# pipestance on retry instead of duplicating tens of gigabytes under a new
# SLURM-job-specific path.
stage_root="$stage_parent/.staging_${sample}"
run_id="cellranger_run"
mkdir -p "$stage_root"
export TMPDIR="$stage_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"

min_free_gb="${CELLRANGER_MIN_FREE_GB:-150}"
if [[ "$min_free_gb" =~ ^[0-9]+$ ]] && ((min_free_gb > 0)); then
  available_kb="$(df -Pk "$stage_parent" | awk 'NR == 2 {print $4}')"
  required_kb=$((min_free_gb * 1024 * 1024))
  if [[ ! "$available_kb" =~ ^[0-9]+$ ]] || ((available_kb < required_kb)); then
    available_gb=$(( ${available_kb:-0} / 1024 / 1024 ))
    echo "ERROR: Cell Ranger needs at least ${min_free_gb} GB free in the staging filesystem before starting; only approximately ${available_gb} GB is available at $stage_parent." >&2
    echo "Existing resumable staging data were preserved at $stage_root." >&2
    exit 28
  fi
fi

args=(count
  "--id=$run_id"
  "--transcriptome=$transcriptome"
  "--fastqs=$fastq_dir"
  # CodeSpringApp consumes the filtered feature-barcode matrix, not the BAM.
  # Disabling BAM creation substantially reduces runtime, I/O, and disk use.
  "--create-bam=false"
  "--include-introns=true"
  "--nosecondary"
  "--localcores=${SLURM_CPUS_PER_TASK:-8}"
  "--localmem=${CELLRANGER_LOCALMEM_GB:-60}"
)
[[ -n "$fastq_sample" ]] && args+=("--sample=$fastq_sample")
if [[ "$expect_cells" =~ ^[0-9]+$ ]] && ((expect_cells > 0)); then args+=("--expect-cells=$expect_cells"); fi

(cd "$stage_root" && cellranger "${args[@]}")
staged="$stage_root/$run_id"
filtered="$staged/outs/filtered_feature_bc_matrix"
for required in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz; do
  [[ -s "$filtered/$required" ]] || { echo "ERROR: Cell Ranger did not create $filtered/$required" >&2; exit 1; }
done

final="$output_root/$sample"
finalizing="$output_root/.finalizing_${sample}_${SLURM_JOB_ID:-$$}"
rm -rf -- "$finalizing"
mkdir -p "$finalizing/outs"

# Retain only the downstream matrix and small human-readable summaries. The
# large Martian pipestance/intermediate tree is deleted after finalization.
mv -- "$filtered" "$finalizing/outs/filtered_feature_bc_matrix"
for artifact in metrics_summary.csv web_summary.html; do
  if [[ -s "$staged/outs/$artifact" ]]; then
    mv -- "$staged/outs/$artifact" "$finalizing/outs/$artifact"
  fi
done
for metadata in _invocation _cmdline arguments.txt; do
  if [[ -s "$staged/$metadata" ]]; then
    mv -- "$staged/$metadata" "$finalizing/$metadata"
  fi
done
cat > "$finalizing/${sample}_cellranger_summary.txt" <<EOF
sample_id	$sample
fastq_folder	$fastq_dir
fastq_sample	$fastq_sample
transcriptome	$transcriptome
expected_cells	$expect_cells
filtered_matrix	$final/outs/filtered_feature_bc_matrix
create_bam	false
staging_directory	$stage_root
EOF
touch "$finalizing/_CELLRANGER_COMPLETE"
rm -rf -- "$final"
mv -- "$finalizing" "$final"
rm -rf -- "$stage_root"
