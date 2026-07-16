#!/usr/bin/env bash
set -euo pipefail

sample="$1"
target_bam="$2"
control_bam="$3"
genome_size="$4"
qval="${5:-0.01}"
peak_type="${6:-narrow}"
outdir="$7"
project_name="$8"

module load EBModules
module load MACS2/2.2.9.1-foss-2022b

if [[ ! -s "$target_bam" ]]; then
  echo "ERROR: CUT&RUN target BAM is missing or empty: ${target_bam}" >&2
  exit 1
fi
if [[ "$control_bam" != "none" && ! -s "$control_bam" ]]; then
  echo "ERROR: CUT&RUN control BAM is missing or empty: ${control_bam}" >&2
  exit 1
fi
if [[ "$peak_type" != "narrow" && "$peak_type" != "broad" ]]; then
  echo "ERROR: CUT&RUN MACS peak type must be narrow or broad, not: ${peak_type}" >&2
  exit 1
fi

mkdir -p "$outdir"

# Python's tempfile module otherwise uses the compute node's shared /tmp,
# which can fill when many MACS jobs run concurrently. Keep temporary files on
# the project filesystem and remove them on both successful and failed exits.
tmp_dir="${outdir}/.macs2_tmp"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"
export TMP="$tmp_dir"
export TEMP="$tmp_dir"
cleanup_tmp() {
  rm -rf "$tmp_dir"
}
trap cleanup_tmp EXIT

args=(callpeak -t "$target_bam" -f BAMPE -g "$genome_size" -q "$qval" --keep-dup all -n "$sample" --outdir "$outdir" --bdg)
if [[ "$control_bam" != "none" && -s "$control_bam" ]]; then
  args+=(-c "$control_bam")
fi
if [[ "$peak_type" == "broad" ]]; then
  args+=(--broad)
fi

run_log="${outdir}/${sample}_macs2.log"
summary="${outdir}/${sample}_macs2_summary.txt"
complete_marker="${outdir}/${sample}_macs2_complete.txt"
rm -f "$complete_marker"
if [[ "$peak_type" == "broad" ]]; then
  peak_file="${outdir}/${sample}_peaks.broadPeak"
else
  peak_file="${outdir}/${sample}_peaks.narrowPeak"
fi

# MACS2 can emit an ignored internal exception, create partial peak files, and
# still exit successfully. Capture stderr and reject such chromosome-level
# failures instead of reporting a partial result as complete.
macs_status=0
macs2 "${args[@]}" 2> "$run_log" || macs_status=$?
cat "$run_log" >&2
if ((macs_status != 0)); then
  echo "ERROR: MACS2 exited with status ${macs_status} for ${sample}. See ${run_log}." >&2
  exit "$macs_status"
fi

fatal_pattern='Traceback \(most recent call last\):|Exception ignored in:|KeyError:|ValueError:|TypeError:|OSError:|No space left on device|Killed|Segmentation fault'
if grep -Eq "$fatal_pattern" "$run_log"; then
  echo "ERROR: MACS2 reported an internal peak-calling exception for ${sample}. See ${run_log}." >&2
  exit 1
fi

if [[ ! -f "$peak_file" ]]; then
  echo "ERROR: MACS2 completed without creating the expected peak file: ${peak_file}" >&2
  exit 1
fi

if [[ ! -s "${outdir}/${sample}_peaks.xls" ]]; then
  echo "ERROR: MACS2 completed without creating a non-empty peaks table for ${sample}." >&2
  exit 1
fi

peak_count="$(wc -l < "$peak_file" | tr -d '[:space:]')"
caller_version="$(macs2 --version 2>&1 | head -n 1)"

{
  echo -e "sample\t${sample}"
  echo -e "status\tcomplete"
  echo -e "caller\tmacs2"
  echo -e "caller_version\t${caller_version}"
  echo -e "target_bam\t${target_bam}"
  echo -e "control_bam\t${control_bam}"
  echo -e "genome_size\t${genome_size}"
  echo -e "qval\t${qval}"
  echo -e "peak_type\t${peak_type}"
  echo -e "peak_file\t${peak_file}"
  echo -e "peak_count\t${peak_count}"
  echo -e "run_log\t${run_log}"
} > "$summary"

marker_tmp="${complete_marker}.tmp.$$"
{
  echo -e "sample\t${sample}"
  echo -e "status\tcomplete"
  echo -e "caller_version\t${caller_version}"
  echo -e "peak_file\t${peak_file}"
  echo -e "peak_count\t${peak_count}"
  echo -e "summary\t${summary}"
} > "$marker_tmp"
mv "$marker_tmp" "$complete_marker"
