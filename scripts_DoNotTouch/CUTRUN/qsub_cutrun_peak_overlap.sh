#!/usr/bin/env bash
#SBATCH --job-name=cutrun_peak_overlap
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=02:00:00

set -euo pipefail
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin${PATH:+:$PATH}"

runner="${7:-}"
if [[ ! -s "$runner" ]]; then
  echo "ERROR: CUT&RUN peak-overlap runner was not found at: ${runner:-<missing>}" >&2
  exit 2
fi

exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6"
