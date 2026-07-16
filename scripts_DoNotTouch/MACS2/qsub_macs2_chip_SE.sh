#!/usr/bin/env bash
#SBATCH --job-name=chip_macs2
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --time=12:00:00

set -euo pipefail
runner="${9:-}"
[[ -s "$runner" ]] || { echo "ERROR: ChIP-seq MACS2 runner not found: ${runner:-<missing>}" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
