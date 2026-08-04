#!/usr/bin/env bash
set -Eeuo pipefail

sample="${1:?ERROR: sample is required}"
target_bam="${2:?ERROR: target BAM is required}"
control_bam="${3:?ERROR: matched input BAM is required}"
input_format="${4:?ERROR: MACS2 input format is required}"
genome_size="${5:?ERROR: MACS2 genome size is required}"
qvalue="${6:-0.01}"
peak_type="${7:-narrow}"
outdir="${8:?ERROR: output directory is required}"

if [[ "$input_format" != "BAM" && "$input_format" != "BAMPE" ]]; then
  echo "ERROR: ChIP-seq MACS2 format must be BAM or BAMPE, not: $input_format" >&2
  exit 2
fi
if [[ "$peak_type" != "narrow" && "$peak_type" != "broad" ]]; then
  echo "ERROR: ChIP-seq peak type must be narrow or broad, not: $peak_type" >&2
  exit 2
fi
[[ -s "$target_bam" ]] || { echo "ERROR: ChIP-seq target BAM is missing or empty: $target_bam" >&2; exit 1; }
[[ -s "$control_bam" ]] || { echo "ERROR: matched input BAM is missing or empty: $control_bam" >&2; exit 1; }

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load MACS2/2.2.9.1-foss-2022b
command -v macs2 >/dev/null 2>&1 || { echo "ERROR: macs2 was not found after loading its module." >&2; exit 127; }

mkdir -p "$outdir"
tmp_dir="${outdir}/.macs2_tmp_${SLURM_JOB_ID:-$$}"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir" TMP="$tmp_dir" TEMP="$tmp_dir"
cleanup() { rm -rf "$tmp_dir"; }
trap cleanup EXIT

args=(callpeak -t "$target_bam" -c "$control_bam" -f "$input_format" -g "$genome_size" -q "$qvalue" --keep-dup all -n "$sample" --outdir "$outdir" --bdg)
if [[ "$peak_type" == "broad" ]]; then
  args+=(--broad)
else
  args+=(--call-summits)
fi

run_log="${outdir}/${sample}_macs2.log"
summary="${outdir}/${sample}_macs2_summary.txt"
complete_marker="${outdir}/${sample}_macs2_complete.txt"
rm -f "$complete_marker"
macs_status=0
model_mode="estimated"
macs2 "${args[@]}" 2> "$run_log" || macs_status=$?
if ((macs_status != 0)) && [[ "$input_format" == "BAM" ]] && grep -Eq 'Too few paired peaks|can not build the model|cannot build the model' "$run_log"; then
  model_mode="fixed_extension_147_fallback"
  fallback_log="${run_log}.fallback"
  printf '\nINFO: MACS2 could not estimate a strand-shift model; retrying with --nomodel --extsize 147.\n' >> "$run_log"
  macs_status=0
  macs2 "${args[@]}" --nomodel --extsize 147 2> "$fallback_log" || macs_status=$?
  cat "$fallback_log" >> "$run_log"
  rm -f "$fallback_log"
fi
cat "$run_log" >&2
if ((macs_status != 0)); then
  echo "ERROR: MACS2 exited with status ${macs_status} for ${sample}. See ${run_log}." >&2
  exit "$macs_status"
fi

fatal_pattern='Traceback \(most recent call last\):|Exception ignored in:|KeyError:|ValueError:|TypeError:|OSError:|No space left on device|Killed|Segmentation fault'
if grep -Eq "$fatal_pattern" "$run_log"; then
  echo "ERROR: MACS2 reported an internal exception for ${sample}. See ${run_log}." >&2
  exit 1
fi

peak_file="${outdir}/${sample}_peaks.$([[ "$peak_type" == "broad" ]] && printf broadPeak || printf narrowPeak)"
[[ -f "$peak_file" ]] || { echo "ERROR: expected MACS2 peak file was not created: $peak_file" >&2; exit 1; }
[[ -s "${outdir}/${sample}_peaks.xls" ]] || { echo "ERROR: MACS2 peaks table is missing or empty for ${sample}." >&2; exit 1; }
peak_count="$(wc -l < "$peak_file" | tr -d '[:space:]')"
caller_version="$(macs2 --version 2>&1 | head -n 1)"
{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tcomplete\n'
  printf 'caller\tmacs2\n'
  printf 'caller_version\t%s\n' "$caller_version"
  printf 'target_bam\t%s\n' "$target_bam"
  printf 'control_bam\t%s\n' "$control_bam"
  printf 'input_format\t%s\n' "$input_format"
  printf 'genome_size\t%s\n' "$genome_size"
  printf 'qvalue\t%s\n' "$qvalue"
  printf 'peak_type\t%s\n' "$peak_type"
  printf 'model_mode\t%s\n' "$model_mode"
  printf 'fixed_extension_bp\t%s\n' "$([[ "$model_mode" == "fixed_extension_147_fallback" ]] && printf 147 || printf NA)"
  printf 'peak_file\t%s\n' "$peak_file"
  printf 'peak_count\t%s\n' "$peak_count"
  printf 'run_log\t%s\n' "$run_log"
} > "$summary"
marker_tmp="${complete_marker}.tmp.$$"
{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tcomplete\n'
  printf 'peak_file\t%s\n' "$peak_file"
  printf 'peak_count\t%s\n' "$peak_count"
  printf 'summary\t%s\n' "$summary"
} > "$marker_tmp"
mv "$marker_tmp" "$complete_marker"
