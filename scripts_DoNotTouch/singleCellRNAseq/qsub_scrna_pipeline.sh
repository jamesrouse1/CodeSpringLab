#!/usr/bin/env bash
#SBATCH --job-name=codespring_scrna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --export=NONE

set -euo pipefail

runner="$1"
engine="$2"
samples="$3"
out_dir="$4"
params="$5"

if [[ ! -x "$runner" ]]; then
  echo "ERROR: scRNA-seq runner is missing or not executable: $runner" >&2
  exit 2
fi
"$runner" "$engine" "$samples" "$out_dir" "$params"
