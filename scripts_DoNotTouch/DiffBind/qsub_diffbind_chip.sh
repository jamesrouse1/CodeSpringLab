#!/usr/bin/env bash
#SBATCH --job-name=chip_diffbind
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail
runner="${8:-}"
[[ -s "$runner" ]] || { echo "ERROR: ChIP DiffBind runner not found: ${runner:-<missing>}" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7"
