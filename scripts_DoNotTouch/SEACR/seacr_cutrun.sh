#!/usr/bin/env bash
set -euo pipefail

target_bedgraph="$1"
control_bedgraph="$2"
norm="${3:-norm}"
stringency="${4:-stringent}"
out_prefix="$5"
project_name="$6"
target_fragments="${7:-none}"

[[ -s "$target_bedgraph" ]] || { echo "ERROR: target bedGraph is missing or empty: $target_bedgraph" >&2; exit 2; }
[[ "$norm" == "norm" || "$norm" == "non" ]] || { echo "ERROR: SEACR normalization must be norm or non, not: $norm" >&2; exit 2; }
[[ "$stringency" == "stringent" || "$stringency" == "relaxed" ]] || { echo "ERROR: SEACR stringency must be stringent or relaxed, not: $stringency" >&2; exit 2; }
if [[ "$control_bedgraph" != "none" && ! -s "$control_bedgraph" ]]; then
  echo "ERROR: control bedGraph is missing or empty: $control_bedgraph" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
seacr_script="${script_dir}/SEACR_1.3.sh"
seacr_r_script="${script_dir}/SEACR_1.3.R"
if [[ ! -s "$seacr_script" || ! -s "$seacr_r_script" ]]; then
  echo "ERROR: The bundled SEACR 1.3 shell and R scripts were not both found in ${script_dir}" >&2
  echo "Install once from the CodeSpringLab folder with:" >&2
  echo "  bash ../scripts_DoNotTouch/SEACR/download_seacr.sh" >&2
  exit 2
fi

out_dir="$(dirname "$out_prefix")"
mkdir -p "$out_dir"
out_dir="$(cd -- "$out_dir" && pwd)"
out_prefix="$out_dir/$(basename "$out_prefix")"
target_bedgraph="$(cd -- "$(dirname "$target_bedgraph")" && pwd)/$(basename "$target_bedgraph")"
if [[ "$control_bedgraph" != "none" ]]; then
  control_bedgraph="$(cd -- "$(dirname "$control_bedgraph")" && pwd)/$(basename "$control_bedgraph")"
fi
seacr_tmp="$(mktemp -d "$out_dir/.seacr_tmp.XXXXXX")"
cleanup() { rm -rf -- "$seacr_tmp"; }
trap cleanup EXIT
ln -s "$target_bedgraph" "$seacr_tmp/target.bedgraph"
if [[ "$control_bedgraph" != "none" ]]; then
  ln -s "$control_bedgraph" "$seacr_tmp/control.bedgraph"
fi
## SLURM jobs use a minimal environment, so the module function may not be
## initialized even though it is available in an interactive login shell.
if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    if [[ -s "$module_init" ]]; then
      # shellcheck disable=SC1090
      source "$module_init"
      break
    fi
  done
fi
if ! type module >/dev/null 2>&1; then
  echo "ERROR: the cluster module command is unavailable in this SLURM job." >&2
  exit 127
fi

module load EBModules
module load BEDTools/2.30.0-GCC-10.3.0
module load R/4.1.2-foss-2021a

if [[ "$control_bedgraph" == "none" || ! -s "$control_bedgraph" ]]; then
  (cd "$seacr_tmp" && bash "$seacr_script" target.bedgraph 0.01 "$norm" "$stringency" result)
else
  (cd "$seacr_tmp" && bash "$seacr_script" target.bedgraph control.bedgraph "$norm" "$stringency" result)
fi
generated_peak="$seacr_tmp/result.${stringency}.bed"
[[ -f "$generated_peak" ]] || { echo "ERROR: SEACR did not create the expected peak file: $generated_peak" >&2; exit 1; }
mv "$generated_peak" "${out_prefix}.${stringency}.bed"

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
