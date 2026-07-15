#!/usr/bin/env bash
set -euo pipefail

reads_csv="${1:-}"
outdir="${2:-}"
if [[ -z "$reads_csv" || -z "$outdir" ]]; then
  echo "ERROR: usage: fastqc.sh <comma-separated FASTQs> <output directory>" >&2
  return 2 2>/dev/null || exit 2
fi

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
  return 127 2>/dev/null || exit 127
fi

module load EBModules
module load FastQC/0.12.1-Java-11
module load Bowtie2/2.4.4-GCC-10.3.0

if ! command -v fastqc >/dev/null 2>&1; then
  echo "ERROR: fastqc was not found after loading FastQC/0.12.1-Java-11." >&2
  return 127 2>/dev/null || exit 127
fi
if ! command -v bowtie2 >/dev/null 2>&1; then
  echo "ERROR: bowtie2 was not found after loading Bowtie2/2.4.4-GCC-10.3.0." >&2
  return 127 2>/dev/null || exit 127
fi

fastq_screen="${FASTQ_SCREEN_BIN:-/grid/bsr/data/data/utama/tools/FastQ-Screen-0.15.2/fastq_screen}"
if [[ ! -x "$fastq_screen" ]]; then
  echo "ERROR: FastQ Screen executable is missing or not executable: $fastq_screen" >&2
  return 127 2>/dev/null || exit 127
fi

mkdir -p "$outdir"
IFS=',' read -r -a reads <<< "$reads_csv"
if [[ ${#reads[@]} -eq 0 ]]; then
  echo "ERROR: no FASTQ files were supplied." >&2
  return 2 2>/dev/null || exit 2
fi

for read in "${reads[@]}"; do
  if [[ ! -s "$read" ]]; then
    echo "ERROR: FASTQ file is missing or empty: $read" >&2
    return 2 2>/dev/null || exit 2
  fi
  echo "Running FastQC: $read"
  fastqc -t "${SLURM_CPUS_PER_TASK:-4}" -o "$outdir" "$read"

  echo "Running FastQ Screen: $read"
  "$fastq_screen" -threads "${SLURM_CPUS_PER_TASK:-4}" -outdir "$outdir" "$read"
done
