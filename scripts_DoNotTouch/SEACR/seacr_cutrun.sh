#!/usr/bin/env bash
set -euo pipefail

target_bedgraph="$1"
control_bedgraph="$2"
norm="${3:-norm}"
stringency="${4:-stringent}"
out_prefix="$5"
project_name="$6"
target_fragments="${7:-none}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
seacr_script="${script_dir}/SEACR_1.3.sh"
seacr_r_script="${script_dir}/SEACR_1.3.R"
if [[ ! -s "$seacr_script" || ! -s "$seacr_r_script" ]]; then
  echo "ERROR: The bundled SEACR 1.3 shell and R scripts were not both found in ${script_dir}" >&2
  echo "Install once from the CodeSpringLab folder with:" >&2
  echo "  bash ../scripts_DoNotTouch/SEACR/download_seacr.sh" >&2
  exit 2
fi

mkdir -p "$(dirname "$out_prefix")"
module load EBModules
module load BEDTools/2.30.0-GCC-10.3.0
module load R/4.1.2-foss-2021a

if [[ "$control_bedgraph" == "none" || ! -s "$control_bedgraph" ]]; then
  bash "$seacr_script" "$target_bedgraph" 0.01 "$norm" "$stringency" "$out_prefix"
else
  bash "$seacr_script" "$target_bedgraph" "$control_bedgraph" "$norm" "$stringency" "$out_prefix"
fi

echo -e "target_bedgraph\t${target_bedgraph}" > "${out_prefix}_seacr_summary.txt"
echo -e "control_bedgraph\t${control_bedgraph}" >> "${out_prefix}_seacr_summary.txt"
echo -e "normalization\t${norm}" >> "${out_prefix}_seacr_summary.txt"
echo -e "stringency\t${stringency}" >> "${out_prefix}_seacr_summary.txt"

peak_bed="${out_prefix}.${stringency}.bed"
echo -e "peak_bed\t${peak_bed}" >> "${out_prefix}_seacr_summary.txt"
echo -e "peak_count\t$(if [[ -s "$peak_bed" ]]; then wc -l < "$peak_bed"; else echo 0; fi)" >> "${out_prefix}_seacr_summary.txt"
if [[ -s "$peak_bed" && "$target_fragments" != "none" && -s "$target_fragments" ]]; then
  total_fragments="$(wc -l < "$target_fragments")"
  fragments_in_peaks="$(bedtools intersect -nonamecheck -u -a "$target_fragments" -b "$peak_bed" | wc -l)"
  frip="$(awk -v a="$fragments_in_peaks" -v b="$total_fragments" 'BEGIN{if (b>0) printf "%.6f", a/b; else print "NA"}')"
  echo -e "target_fragments\t${target_fragments}" >> "${out_prefix}_seacr_summary.txt"
  echo -e "total_fragments\t${total_fragments}" >> "${out_prefix}_seacr_summary.txt"
  echo -e "fragments_in_peaks\t${fragments_in_peaks}" >> "${out_prefix}_seacr_summary.txt"
  echo -e "frip\t${frip}" >> "${out_prefix}_seacr_summary.txt"
fi
