#!/usr/bin/env bash
#SBATCH --job-name=cutrun_diffbind
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail

minimum_peaks_per_sample="${11:-1}"
peak_source="${12:-legacy}"
runner="${13:-}"
if [[ -z "$runner" && -n "${11:-}" && -s "${11:-}" ]]; then
  # Backward compatibility for jobs submitted before peak-source selection was
  # added, where argument 11 was the runner path.
  runner="${11}"
  minimum_peaks_per_sample=1
  peak_source=legacy
fi
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/DiffBind/cutrun_diffbind.sh"
fi
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN DiffBind runner was not found at: $runner" >&2
  exit 2
fi

source "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "${7:-none}" "${8:-}" "${9:-}" "${10:-}" "$minimum_peaks_per_sample" "$peak_source"
