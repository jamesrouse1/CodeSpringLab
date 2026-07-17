#!/usr/bin/env bash
#SBATCH --job-name=featurecounts
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -Eeuo pipefail
runner="${7:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/featureCounts/featurecounts_PE.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: paired-end featureCounts runner was not found: $runner" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5"
