#!/usr/bin/env bash
#SBATCH --job-name=cutrun_seacr
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=12:00:00

seacr_runner="${8:-}"
if [[ -z "$seacr_runner" || ! -s "$seacr_runner" ]]; then
  seacr_runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/SEACR/seacr_cutrun.sh"
fi
if [[ ! -s "$seacr_runner" ]]; then
  echo "ERROR: SEACR runner was not found at: $seacr_runner" >&2
  exit 2
fi
source "$seacr_runner" "$1" "$2" "$3" "$4" "$5" "$6" "${7:-none}"
