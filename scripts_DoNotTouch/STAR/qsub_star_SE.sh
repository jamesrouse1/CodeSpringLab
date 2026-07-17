#!/usr/bin/env bash
#SBATCH --job-name=star
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -Eeuo pipefail
runner="${5:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/STAR/star_SE.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: single-end STAR runner was not found: $runner" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3"
