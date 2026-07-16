#!/usr/bin/env bash
#SBATCH --job-name=cutadapt
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=12:00:00

set -euo pipefail
runner="${9:-}"
[[ -s "$runner" ]] || { echo "ERROR: single-end Cutadapt runner not found: ${runner:-<missing>}" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7"
