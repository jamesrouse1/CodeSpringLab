#!/usr/bin/env bash
#SBATCH --job-name=peak_annotation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --export=NONE

set -Eeuo pipefail

runner="${4:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  echo "ERROR: peak-annotation runner was not found at: ${runner:-<missing>}" >&2
  exit 2
fi

exec bash "$runner" "$1" "$2" "$3"
