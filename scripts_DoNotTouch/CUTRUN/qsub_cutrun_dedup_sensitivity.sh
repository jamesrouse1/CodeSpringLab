#!/usr/bin/env bash
#SBATCH --job-name=cutrun_dedup_sensitivity
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=12:00:00

runner="${12:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/CUTRUN/cutrun_dedup_sensitivity.sh"
fi
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN deduplicated-target sensitivity runner was not found at: $runner" >&2
  exit 2
fi

source "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}"
