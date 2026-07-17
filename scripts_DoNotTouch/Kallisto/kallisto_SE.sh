#!/usr/bin/env bash
set -Eeuo pipefail

outdir="${1:-}"
index="${2:-}"
reads_csv="${3:-}"
if [[ -z "$outdir" || -z "$index" || -z "$reads_csv" ]]; then
  echo "ERROR: usage: kallisto_SE.sh <output directory> <index> <FASTQ(s)>" >&2
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
IFS=',' read -r -a input_reads <<< "$reads_csv"
for file in "${input_reads[@]}"; do
  [[ -s "$file" ]] || { echo "ERROR: FASTQ input is missing or empty: $file" >&2; exit 2; }
done

mkdir -p "$outdir"
rm -f -- "$outdir/abundance.tsv" "$outdir/abundance.h5" "$outdir/run_info.json"
ulimit -n 10000
if [[ "$reads_csv" == *,* ]]; then
  tmpdir="$(mktemp -d "$outdir/.kallisto_tmp.XXXXXX")"
  fifo="$tmpdir/R1.fastq"
  mkfifo "$fifo"
  cleanup() {
    [[ -n "${stream_pid:-}" ]] && kill "$stream_pid" 2>/dev/null || true
    rm -rf -- "$tmpdir"
  }
  trap cleanup EXIT
  stream_fastqs "$reads_csv" > "$fifo" & stream_pid=$!
  kallisto quant -i "$index" -t "${SLURM_CPUS_PER_TASK:-4}" -o "$outdir" -b 100 --single -l 180 -s 20 "$fifo"
  wait "$stream_pid"
else
  kallisto quant -i "$index" -t "${SLURM_CPUS_PER_TASK:-4}" -o "$outdir" -b 100 --single -l 180 -s 20 "$reads_csv"
fi
[[ -s "$outdir/abundance.tsv" ]] || { echo "ERROR: Kallisto did not create a non-empty abundance.tsv in $outdir" >&2; exit 1; }
[[ -s "$outdir/run_info.json" ]] || { echo "ERROR: Kallisto did not create run_info.json in $outdir" >&2; exit 1; }
