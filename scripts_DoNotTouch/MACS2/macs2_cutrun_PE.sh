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

mkdir -p "$outdir"

args=(callpeak -t "$target_bam" -f BAMPE -g "$genome_size" -q "$qval" --keep-dup all -n "$sample" --outdir "$outdir" --bdg)
if [[ "$control_bam" != "none" && -s "$control_bam" ]]; then
  args+=(-c "$control_bam")
fi
if [[ "$peak_type" == "broad" ]]; then
  args+=(--broad)
fi

macs2 "${args[@]}"

{
  echo -e "sample\t${sample}"
  echo -e "target_bam\t${target_bam}"
  echo -e "control_bam\t${control_bam}"
  echo -e "genome_size\t${genome_size}"
  echo -e "qval\t${qval}"
  echo -e "peak_type\t${peak_type}"
} > "${outdir}/${sample}_macs2_summary.txt"
