#!/usr/bin/env bash
set -Eeuo pipefail

genome_bam="${1:-}"
rsem_index="${2:-}"
feature="${3:-}"
out_prefix="${4:-}"
strand_bed="${5:-}"
transcript_bam="${6:-}"
if [[ -z "$genome_bam" || -z "$rsem_index" || -z "$out_prefix" || -z "$strand_bed" || -z "$transcript_bam" ]]; then
  echo "ERROR: usage: RSEM_PE.sh <genome BAM> <RSEM index prefix> <feature> <output prefix> <strand BED> <transcriptome BAM>" >&2
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
module load STAR/2.7.10a-GCC-10.3.0
module load RSEM/1.3.3-foss-2019b
module load SAMtools/1.16.1-GCC-11.3.0
command -v infer_experiment.py >/dev/null 2>&1 || { echo "ERROR: infer_experiment.py was not found after loading RSeQC." >&2; exit 127; }
command -v rsem-calculate-expression >/dev/null 2>&1 || { echo "ERROR: rsem-calculate-expression was not found after loading RSEM." >&2; exit 127; }
[[ -s "$genome_bam" ]] || { echo "ERROR: genome-coordinate BAM is missing or empty: $genome_bam" >&2; exit 2; }
[[ -s "$transcript_bam" ]] || { echo "ERROR: transcriptome BAM is missing or empty: $transcript_bam" >&2; exit 2; }
[[ -s "$strand_bed" ]] || { echo "ERROR: strand-detection BED is missing or empty: $strand_bed" >&2; exit 2; }
if [[ ! -s "${rsem_index}.grp" && ! -s "${rsem_index}.ti" ]]; then
  echo "ERROR: RSEM index prefix does not have a readable .grp or .ti file: $rsem_index" >&2
  exit 2
fi

mkdir -p "$(dirname "$out_prefix")"
strand_report="${out_prefix}_strand.txt"
genes_result="${out_prefix}.genes.results"
isoforms_result="${out_prefix}.isoforms.results"
rm -f -- "$genes_result" "$isoforms_result"
infer_experiment.py -r "$strand_bed" -i "$genome_bam" > "$strand_report"
fw="$(awk 'NR==5 {print $NF}' "$strand_report")"
rv="$(awk 'NR==6 {print $NF}' "$strand_report")"
if ! awk -v fw="$fw" -v rv="$rv" 'BEGIN {exit !(fw ~ /^[0-9.]+$/ && rv ~ /^[0-9.]+$/ && fw + rv > 0)}'; then
  echo "ERROR: could not determine RSEM library strandedness from $strand_report" >&2
  exit 1
fi
strand="$(awk -v fw="$fw" -v rv="$rv" 'BEGIN {d=(fw-rv)*100/(fw+rv); if (d<0) a=-d; else a=d; if (a<50) print "none"; else if (d>0) print "forward"; else print "reverse"}')"

ulimit -n 10000
rsem-calculate-expression --paired-end -p "${SLURM_CPUS_PER_TASK:-8}" --strandedness "$strand" \
  --bam "$transcript_bam" "$rsem_index" "$out_prefix"
[[ -s "$genes_result" ]] || { echo "ERROR: RSEM did not create a non-empty genes result: $genes_result" >&2; exit 1; }
[[ -s "$isoforms_result" ]] || { echo "ERROR: RSEM did not create a non-empty isoforms result: $isoforms_result" >&2; exit 1; }
rm -f -- "${out_prefix}.transcript.bam"
rm -rf -- "${out_prefix}.stat"

