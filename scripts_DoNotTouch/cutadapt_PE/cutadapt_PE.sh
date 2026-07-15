#!/usr/bin/env bash
set -euo pipefail

module load EBModules
module load cutadapt/4.4-GCCcore-12.2.0

#module load Anaconda3/2022.05
#conda init bash
#conda activate cutadapt


min_length="$1"
adapter1="$2"
adapter2="$3"
output1="$4"
output2="$5"
read1="$6"
read2="$7"

stream_fastqs() {
  local csv="$1"
  local file
  local -a inputs
  IFS=',' read -r -a inputs <<< "$csv"
  for file in "${inputs[@]}"; do
    if [[ ! -s "$file" ]]; then
      echo "ERROR: pooled FASTQ input is missing or empty: $file" >&2
      return 2
    fi
    if [[ "$file" == *.gz ]]; then gzip -cd -- "$file"; else cat -- "$file"; fi
  done
}

mkdir -p "$(dirname "$output1")" "$(dirname "$output2")"
if [[ "$read1" == *,* || "$read2" == *,* ]]; then
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/codespring_cutadapt.XXXXXX")"
  fifo1="$tmpdir/R1.fastq"
  fifo2="$tmpdir/R2.fastq"
  mkfifo "$fifo1" "$fifo2"
  cleanup() {
    jobs -pr | xargs -r kill 2>/dev/null || true
    rm -rf "$tmpdir"
  }
  trap cleanup EXIT
  stream_fastqs "$read1" > "$fifo1" & pid1=$!
  stream_fastqs "$read2" > "$fifo2" & pid2=$!
  cutadapt -j 4 -m "$min_length" -a "$adapter1" -A "$adapter2" -o "$output1" -p "$output2" "$fifo1" "$fifo2"
  wait "$pid1" "$pid2"
else
  cutadapt -j 4 -m "$min_length" -a "$adapter1" -A "$adapter2" -o "$output1" -p "$output2" "$read1" "$read2"
fi

#conda deactivate
