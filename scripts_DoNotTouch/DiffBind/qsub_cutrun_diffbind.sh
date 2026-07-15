#!/usr/bin/env bash
#SBATCH --job-name=cutrun_diffbind
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail

runner="${7:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/DiffBind/cutrun_diffbind.sh"
fi
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN DiffBind runner was not found at: $runner" >&2
  exit 2
fi

source "$runner" "$1" "$2" "$3" "$4" "$5" "$6"
