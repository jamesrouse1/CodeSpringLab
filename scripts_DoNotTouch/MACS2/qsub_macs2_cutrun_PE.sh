#!/usr/bin/env bash
#SBATCH --job-name=cutrun_macs2
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --time=12:00:00

runner="${9:-}"
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN MACS runner was not found at: ${runner:-<missing>}" >&2
  exit 2
fi

exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
