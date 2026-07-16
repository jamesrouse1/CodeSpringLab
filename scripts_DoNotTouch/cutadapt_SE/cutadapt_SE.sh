#!/usr/bin/env bash
set -Eeuo pipefail

min_length="${1:?ERROR: minimum read length is required}"
adapter="${2:?ERROR: read adapter is required}"
output="${4:?ERROR: output FASTQ is required}"
reads_csv="${6:?ERROR: input FASTQ is required}"

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load cutadapt/4.4-GCCcore-12.2.0
command -v cutadapt >/dev/null 2>&1 || { echo "ERROR: cutadapt was not found after loading its module." >&2; exit 127; }

mkdir -p "$(dirname "$output")"
stream_fastqs() {
  local csv="$1" file
  local -a inputs
  IFS=',' read -r -a inputs <<< "$csv"
  for file in "${inputs[@]}"; do
    [[ -s "$file" ]] || { echo "ERROR: pooled FASTQ input is missing or empty: $file" >&2; return 2; }
    if [[ "$file" == *.gz ]]; then gzip -cd -- "$file"; else command cat -- "$file"; fi
  done
}

if [[ "$reads_csv" == *,* ]]; then
  tmpdir="$(mktemp -d "$(dirname "$output")/.cutadapt_tmp.XXXXXX")"
  fifo="$tmpdir/R1.fastq"
  mkfifo "$fifo"
  cleanup() {
    [[ -n "${stream_pid:-}" ]] && kill "$stream_pid" 2>/dev/null || true
    rm -rf "$tmpdir"
  }
  trap cleanup EXIT
  stream_fastqs "$reads_csv" > "$fifo" & stream_pid=$!
  cutadapt -j "${SLURM_CPUS_PER_TASK:-4}" -m "$min_length" -a "$adapter" -o "$output" "$fifo"
  wait "$stream_pid"
else
  [[ -s "$reads_csv" ]] || { echo "ERROR: FASTQ input is missing or empty: $reads_csv" >&2; exit 2; }
  cutadapt -j "${SLURM_CPUS_PER_TASK:-4}" -m "$min_length" -a "$adapter" -o "$output" "$reads_csv"
fi
