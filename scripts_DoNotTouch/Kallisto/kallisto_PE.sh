#!/usr/bin/env bash
set -Eeuo pipefail

outdir="${1:-}"
index="${2:-}"
read1="${3:-}"
read2="${4:-}"
if [[ -z "$outdir" || -z "$index" || -z "$read1" || -z "$read2" ]]; then
  echo "ERROR: usage: kallisto_PE.sh <output directory> <index> <R1 FASTQ(s)> <R2 FASTQ(s)>" >&2
  exit 2
fi

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load kallisto/0.46.1-foss-2019b
command -v kallisto >/dev/null 2>&1 || { echo "ERROR: kallisto was not found after loading its module." >&2; exit 127; }
[[ -s "$index" ]] || { echo "ERROR: Kallisto index is missing or empty: $index" >&2; exit 2; }

stream_fastqs() {
  local csv="$1" file
  local -a inputs
  IFS=',' read -r -a inputs <<< "$csv"
  for file in "${inputs[@]}"; do
    [[ -s "$file" ]] || { echo "ERROR: pooled FASTQ input is missing or empty: $file" >&2; return 2; }
    if [[ "$file" == *.gz ]]; then gzip -cd -- "$file"; else cat -- "$file"; fi
  done
}
validate_fastqs() {
  local csv="$1" file
  local -a inputs
  IFS=',' read -r -a inputs <<< "$csv"
  for file in "${inputs[@]}"; do
    [[ -s "$file" ]] || { echo "ERROR: FASTQ input is missing or empty: $file" >&2; return 2; }
  done
}
validate_fastqs "$read1"
validate_fastqs "$read2"

mkdir -p "$outdir"
rm -f -- "$outdir/abundance.tsv" "$outdir/abundance.h5" "$outdir/run_info.json"
ulimit -n 10000
if [[ "$read1" == *,* || "$read2" == *,* ]]; then
  tmpdir="$(mktemp -d "$outdir/.kallisto_tmp.XXXXXX")"
  fifo1="$tmpdir/R1.fastq"
  fifo2="$tmpdir/R2.fastq"
  mkfifo "$fifo1" "$fifo2"
  cleanup() {
    jobs -pr | xargs -r kill 2>/dev/null || true
    rm -rf -- "$tmpdir"
  }
  trap cleanup EXIT
  stream_fastqs "$read1" > "$fifo1" & pid1=$!
  stream_fastqs "$read2" > "$fifo2" & pid2=$!
  kallisto quant -i "$index" -t "${SLURM_CPUS_PER_TASK:-4}" -o "$outdir" -b 100 "$fifo1" "$fifo2"
  wait "$pid1" "$pid2"
else
  kallisto quant -i "$index" -t "${SLURM_CPUS_PER_TASK:-4}" -o "$outdir" -b 100 "$read1" "$read2"
fi
[[ -s "$outdir/abundance.tsv" ]] || { echo "ERROR: Kallisto did not create a non-empty abundance.tsv in $outdir" >&2; exit 1; }
[[ -s "$outdir/run_info.json" ]] || { echo "ERROR: Kallisto did not create run_info.json in $outdir" >&2; exit 1; }
