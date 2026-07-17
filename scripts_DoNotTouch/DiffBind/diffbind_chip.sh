#!/usr/bin/env bash
set -Eeuo pipefail

r_script="${1:?ERROR: ChIP DiffBind R script is required}"
outdir="${2:?ERROR: output directory is required}"
sample_sheet="${3:?ERROR: sample sheet is required}"
reference="${4:?ERROR: reference condition is required}"
comparison="${5:?ERROR: comparison condition is required}"
genome="${6:?ERROR: genome is required}"
blacklist="${7:-none}"

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
[[ -s "$r_script" ]] || { echo "ERROR: ChIP DiffBind R script is missing: $r_script" >&2; exit 2; }
[[ -s "$sample_sheet" ]] || { echo "ERROR: ChIP DiffBind sample sheet is missing: $sample_sheet" >&2; exit 2; }

mkdir -p "$outdir"
run_started="${outdir}/_RUN_STARTED"
complete_marker="${outdir}/_COMPLETE"
rm -f "$complete_marker"
printf 'status\trunning\nstarted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$run_started"
tmp_dir="${outdir}/.diffbind_tmp_${SLURM_JOB_ID:-$$}"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir" TMP="$tmp_dir" TEMP="$tmp_dir"
cleanup() {
  rm -rf "$tmp_dir"
  if [[ -e "$run_started" ]]; then
    rm -f "${outdir}/_DIFFBIND_COMPLETE"
  fi
}
trap cleanup EXIT

module load EBModules
module load R/4.3.2-gfbf-2023a
rm -f "${outdir}/_DIFFBIND_COMPLETE"
Rscript "$r_script" "$outdir" "$sample_sheet" "$reference" "$comparison" "$genome" "$blacklist"
[[ -s "${outdir}/_DIFFBIND_COMPLETE" ]] || { echo "ERROR: ChIP DiffBind did not finish its statistical analysis." >&2; exit 1; }

prefix="DifferentialPeaks_${comparison}_vs_${reference}_ref"
bed_file="${outdir}/${prefix}.with_stats.bed"
annotation="${outdir}/${prefix}_annotated_with_stats.txt"
if [[ ! -s "$bed_file" ]]; then
  printf 'Status\tNo differential peaks passed the DiffBind reporting threshold\nComparison\t%s vs %s\n' "$comparison" "$reference" > "$annotation"
  rm -f "${outdir}/_DIFFBIND_COMPLETE"
  rm -f "$run_started"
  printf 'status\tcomplete\ncompleted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$complete_marker"
  exit 0
fi

case "$genome" in
  mouse) homer_genome="mm39" ;;
  human) homer_genome="hg38" ;;
  *) echo "ERROR: unsupported Homer genome: $genome" >&2; exit 2 ;;
esac
module unload R/4.3.2-gfbf-2023a >/dev/null 2>&1 || true
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a
export PATH="$PATH:/grid/bsr/data/data/utama/tools/homer/bin"
annotatePeaks.pl "$bed_file" "$homer_genome" > "$annotation"
Rscript "$(cd -- "$(dirname -- "$0")/../Homer" && pwd)/expand_peakid_stats.R" "$annotation"
[[ -s "$annotation" ]] || { echo "ERROR: ChIP differential-peak annotation is empty." >&2; exit 1; }
rm -f "${outdir}/_DIFFBIND_COMPLETE"
rm -f "$run_started"
printf 'status\tcomplete\ncompleted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$complete_marker"
