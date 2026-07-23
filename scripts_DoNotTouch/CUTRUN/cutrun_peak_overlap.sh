#!/usr/bin/env bash
set -euo pipefail

# Build the *genomic intersection* of two completed per-sample peak sets.
# The output is a three-column BED containing only bases supported by both
# callers/settings, merged when adjacent intersection fragments touch.
source_a="$1"
source_b="$2"
out_bed="$3"
overlap_name="$4"
sample="$5"
minimum_reciprocal_overlap="${6:-0}"

export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin${PATH:+:$PATH}"

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
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools was not found after loading BEDTools." >&2; exit 127; }

[[ -s "$source_a" ]] || { echo "ERROR: first peak file is missing or empty: $source_a" >&2; exit 2; }
[[ -s "$source_b" ]] || { echo "ERROR: second peak file is missing or empty: $source_b" >&2; exit 2; }
[[ "$minimum_reciprocal_overlap" =~ ^(0|0\.[0-9]+|1(\.0+)?)$ ]] || { echo "ERROR: minimum reciprocal overlap must be between 0 and 1." >&2; exit 2; }

out_dir="$(dirname "$out_bed")"
mkdir -p "$out_dir"
out_dir="$(cd "$out_dir" && pwd)"
out_bed="$out_dir/$(basename "$out_bed")"
summary="${out_bed%.bed}_summary.txt"
no_peaks_marker="${out_bed%.bed}_no_called_overlap_peaks.txt"
tmp_dir="$(mktemp -d "$out_dir/.peak_overlap_tmp.XXXXXX")"
cleanup() { rm -rf "$tmp_dir"; }
trap cleanup EXIT

# BEDTools operates on sorted coordinate inputs.  It emits the actual shared
# genomic portions; merge prevents duplicate rows when one peak overlaps more
# than one peak in the other caller's set.
awk 'BEGIN{OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3>$2 {print $1,$2,$3}' "$source_a" \
  | LC_ALL=C sort -k1,1 -k2,2n > "$tmp_dir/source_a.bed"
awk 'BEGIN{OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3>$2 {print $1,$2,$3}' "$source_b" \
  | LC_ALL=C sort -k1,1 -k2,2n > "$tmp_dir/source_b.bed"

intersect_args=(-a "$tmp_dir/source_a.bed" -b "$tmp_dir/source_b.bed")
if awk -v value="$minimum_reciprocal_overlap" 'BEGIN{exit !(value>0)}'; then
  intersect_args+=(-f "$minimum_reciprocal_overlap" -r)
fi
bedtools intersect "${intersect_args[@]}" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$tmp_dir/overlap.bed"
mv "$tmp_dir/overlap.bed" "$out_bed"

source_a_count="$(wc -l < "$tmp_dir/source_a.bed" | tr -d '[:space:]')"
source_b_count="$(wc -l < "$tmp_dir/source_b.bed" | tr -d '[:space:]')"
overlap_count="$(wc -l < "$out_bed" | tr -d '[:space:]')"
{
  echo -e "sample\t${sample}"
  echo -e "overlap_name\t${overlap_name}"
  echo -e "source_a\t${source_a}"
  echo -e "source_b\t${source_b}"
  echo -e "minimum_reciprocal_overlap\t${minimum_reciprocal_overlap}"
  echo -e "source_a_peaks\t${source_a_count}"
  echo -e "source_b_peaks\t${source_b_count}"
  echo -e "overlap_peaks\t${overlap_count}"
  echo -e "overlap_bed\t${out_bed}"
} > "$summary"

if [[ "$overlap_count" -eq 0 ]]; then
  cat > "$no_peaks_marker" <<EOF
Peak overlap completed successfully but no shared intervals were called.
Sample: ${sample}
Source A: ${source_a}
Source B: ${source_b}
Minimum reciprocal overlap: ${minimum_reciprocal_overlap}
EOF
else
  rm -f "$no_peaks_marker"
fi

