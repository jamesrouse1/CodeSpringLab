#!/usr/bin/env bash
#SBATCH --job-name=cutrun_summary
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=02:00:00

set -euo pipefail

runner="$5"
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN summary runner was not found: $runner" >&2
  exit 2
fi

source "$runner" "$1" "$2" "$3" "$4"
