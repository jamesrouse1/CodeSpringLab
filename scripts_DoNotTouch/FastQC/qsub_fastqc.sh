#!/usr/bin/env bash
#SBATCH --job-name=fastQC
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00

set -euo pipefail

runner="${4:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/FastQC/fastqc.sh"
fi
if [[ ! -s "$runner" ]]; then
  echo "ERROR: FastQC runner was not found at: $runner" >&2
  exit 2
fi

source "$runner" "$1" "$2"

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${3}/log/output_fastQC.txt -e ../../csl_results/${3}/log/error_fastQC.txt -V ../scripts_DoNotTouch/FastQC/fastqc.sh $1 $2
