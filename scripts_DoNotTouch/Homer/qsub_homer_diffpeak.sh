#!/usr/bin/env bash
#SBATCH --job-name=homer_diffpeak
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail
runner="${8:-${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/Homer/homer_diffpeak.sh}"
[[ -s "$runner" ]] || { echo "ERROR: HOMER differential runner is missing: $runner" >&2; exit 2; }
exec bash "$runner" "${1:?}" "${2:?}" "${3:?}" "${4:?}" "${5:?}" "${6:?}"

# Ori script with grid qsub
#qsub  -l mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${7}/log/output_homer_diffpeak.txt -e ../../csl_results/${7}/log/error_homer_diffpeak.txt ../scripts_DoNotTouch/Homer/homer_diffpeak.sh $1 $2 $3 $4 $5 $6
