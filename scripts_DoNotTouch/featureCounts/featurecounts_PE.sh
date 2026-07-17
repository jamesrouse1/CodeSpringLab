#!/usr/bin/env bash
set -Eeuo pipefail

bam="${1:-}"
gtf="${2:-}"
feature="${3:-}"
out_prefix="${4:-}"
strand_bed="${5:-}"
if [[ -z "$bam" || -z "$gtf" || -z "$feature" || -z "$out_prefix" || -z "$strand_bed" ]]; then
  echo "ERROR: usage: featurecounts_PE.sh <BAM> <GTF> <feature attribute> <output prefix> <strand BED>" >&2
  exit 2
fi

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load RSeQC/4.0.0-foss-2021b
module load Subread
command -v infer_experiment.py >/dev/null 2>&1 || { echo "ERROR: infer_experiment.py was not found after loading RSeQC." >&2; exit 127; }
command -v featureCounts >/dev/null 2>&1 || { echo "ERROR: featureCounts was not found after loading Subread." >&2; exit 127; }
[[ -s "$bam" ]] || { echo "ERROR: input BAM is missing or empty: $bam" >&2; exit 2; }
[[ -s "$gtf" ]] || { echo "ERROR: annotation GTF is missing or empty: $gtf" >&2; exit 2; }
[[ -s "$strand_bed" ]] || { echo "ERROR: strand-detection BED is missing or empty: $strand_bed" >&2; exit 2; }

mkdir -p "$(dirname "$out_prefix")"
strand_report="${out_prefix}_strand.txt"
counts_file="${out_prefix}_counts.txt"
rm -f -- "$counts_file" "${counts_file}.summary"
infer_experiment.py -r "$strand_bed" -i "$bam" > "$strand_report"
fw="$(awk 'NR==5 {print $NF}' "$strand_report")"
rv="$(awk 'NR==6 {print $NF}' "$strand_report")"
if ! awk -v fw="$fw" -v rv="$rv" 'BEGIN {exit !(fw ~ /^[0-9.]+$/ && rv ~ /^[0-9.]+$/ && fw + rv > 0)}'; then
  echo "ERROR: could not determine library strandedness from $strand_report" >&2
  exit 1
fi
strand_idx="$(awk -v fw="$fw" -v rv="$rv" 'BEGIN {d=(fw-rv)*100/(fw+rv); if (d<0) a=-d; else a=d; if (a<20) print 0; else if (d>0) print 1; else print 2}')"

featureCounts -a "$gtf" \
  -T "${SLURM_CPUS_PER_TASK:-2}" \
  -p --countReadPairs \
  -t exon \
  -g "$feature" \
  -s "$strand_idx" \
  -Q 12 \
  --minOverlap 1 \
  -C \
  -o "$counts_file" \
  "$bam"

[[ -s "$counts_file" ]] || { echo "ERROR: featureCounts did not create a non-empty count table: $counts_file" >&2; exit 1; }

