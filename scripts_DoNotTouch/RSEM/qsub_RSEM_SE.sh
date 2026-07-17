#!/usr/bin/env bash
#SBATCH --job-name=rsem
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -Eeuo pipefail
runner="${8:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/RSEM/RSEM_SE.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: single-end RSEM runner was not found: $runner" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6"
