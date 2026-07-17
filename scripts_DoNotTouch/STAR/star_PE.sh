#!/usr/bin/env bash
set -Eeuo pipefail

out_prefix="${1:-}"
genome_dir="${2:-}"
read1="${3:-}"
read2="${4:-}"
if [[ -z "$out_prefix" || -z "$genome_dir" || -z "$read1" || -z "$read2" ]]; then
  echo "ERROR: usage: star_PE.sh <output prefix> <STAR index> <R1 FASTQ(s)> <R2 FASTQ(s)>" >&2
  exit 2
fi

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load STAR/2.7.10a-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
command -v STAR >/dev/null 2>&1 || { echo "ERROR: STAR was not found after loading its module." >&2; exit 127; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools was not found after loading its module." >&2; exit 127; }

[[ -d "$genome_dir" ]] || { echo "ERROR: STAR genome index directory is missing: $genome_dir" >&2; exit 2; }
validate_fastq_list() {
  local csv="$1" label="$2" file
  local -a files
  IFS=',' read -r -a files <<< "$csv"
  ((${#files[@]})) || { echo "ERROR: no $label FASTQs were supplied." >&2; return 2; }
  for file in "${files[@]}"; do
    [[ -s "$file" ]] || { echo "ERROR: $label FASTQ is missing or empty: $file" >&2; return 2; }
  done
}
validate_fastq_list "$read1" "R1"
validate_fastq_list "$read2" "R2"

all_reads="${read1},${read2}"
IFS=',' read -r -a read_files <<< "$all_reads"
compressed=0
plain=0
for file in "${read_files[@]}"; do
  if [[ "$file" == *.gz ]]; then compressed=1; else plain=1; fi
done
if ((compressed && plain)); then
  echo "ERROR: compressed and uncompressed FASTQs cannot be mixed in one STAR job." >&2
  exit 2
fi
read_command="cat"
((compressed)) && read_command="zcat"

mkdir -p "$(dirname "$out_prefix")"
aligned_bam="${out_prefix}Aligned.sortedByCoord.out.bam"
aligned_bai="${aligned_bam}.bai"
rm -f -- "$aligned_bam" "$aligned_bai"
ulimit -n 10000

STAR --runThreadN "${SLURM_CPUS_PER_TASK:-4}" \
  --quantMode TranscriptomeSAM \
  --outFileNamePrefix "$out_prefix" \
  --genomeLoad NoSharedMemory \
  --genomeDir "$genome_dir" \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 2 \
  --outSAMunmapped None \
  --outSAMstrandField None \
  --readFilesCommand "$read_command" \
  --readFilesIn "$read1" "$read2"

[[ -s "$aligned_bam" ]] || { echo "ERROR: STAR did not create a non-empty aligned BAM: $aligned_bam" >&2; exit 1; }
samtools quickcheck -v "$aligned_bam"
samtools index -b "$aligned_bam" "$aligned_bai"
[[ -s "$aligned_bai" ]] || { echo "ERROR: samtools did not create a BAM index: $aligned_bai" >&2; exit 1; }
