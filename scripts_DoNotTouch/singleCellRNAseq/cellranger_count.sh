#!/usr/bin/env bash
set -Eeuo pipefail

sample="${1:-}"
fastq_dir="${2:-}"
fastq_sample="${3:-}"
transcriptome="${4:-}"
output_root="${5:-}"
expect_cells="${6:-0}"

if [[ -z "$sample" || -z "$fastq_dir" || -z "$transcriptome" || -z "$output_root" ]]; then
  echo "ERROR: usage: cellranger_count.sh <sample> <FASTQ folder> <FASTQ sample prefix or blank> <transcriptome> <output root> [expected cells]" >&2
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

mkdir -p "$output_root"
stage_root="$output_root/.staging_${sample}_${SLURM_JOB_ID:-$$}"
run_id="cellranger_run"
rm -rf -- "$stage_root"
mkdir -p "$stage_root"
# Cell Ranger can create large temporary files. Keep them with the project
# instead of relying on a compute node's small shared /tmp filesystem.
export TMPDIR="$stage_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"

args=(count
  "--id=$run_id"
  "--transcriptome=$transcriptome"
  "--fastqs=$fastq_dir"
  "--create-bam=true"
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
rm -rf -- "$final"
mv -- "$staged" "$final"
rm -rf -- "$stage_root"
cat > "$final/${sample}_cellranger_summary.txt" <<EOF
sample_id	$sample
fastq_folder	$fastq_dir
fastq_sample	$fastq_sample
transcriptome	$transcriptome
expected_cells	$expect_cells
filtered_matrix	$final/outs/filtered_feature_bc_matrix
EOF
touch "$final/_CELLRANGER_COMPLETE"
