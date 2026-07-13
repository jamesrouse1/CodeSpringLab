#!/usr/bin/env bash
set -euo pipefail

target_bedgraph="$1"
control_bedgraph="$2"
norm="${3:-norm}"
stringency="${4:-stringent}"
out_prefix="$5"
project_name="$6"

seacr_script="../scripts_DoNotTouch/SEACR/SEACR_1.3.sh"
if [[ ! -s "$seacr_script" ]]; then
  echo "ERROR: SEACR was not found at $seacr_script" >&2
  echo "Install once from the CodeSpringLab folder with:" >&2
  echo "  bash ../scripts_DoNotTouch/SEACR/download_seacr.sh" >&2
  exit 2
fi

mkdir -p "$(dirname "$out_prefix")"
module load EBModules || true
module load BEDTools/2.30.0-GCC-10.3.0 || true
module load R/4.1.2-foss-2021a || true

if [[ "$control_bedgraph" == "none" || ! -s "$control_bedgraph" ]]; then
  bash "$seacr_script" "$target_bedgraph" 0.01 "$norm" "$stringency" "$out_prefix"
else
  bash "$seacr_script" "$target_bedgraph" "$control_bedgraph" "$norm" "$stringency" "$out_prefix"
fi

echo -e "target_bedgraph\t${target_bedgraph}" > "${out_prefix}_seacr_summary.txt"
echo -e "control_bedgraph\t${control_bedgraph}" >> "${out_prefix}_seacr_summary.txt"
echo -e "normalization\t${norm}" >> "${out_prefix}_seacr_summary.txt"
echo -e "stringency\t${stringency}" >> "${out_prefix}_seacr_summary.txt"
