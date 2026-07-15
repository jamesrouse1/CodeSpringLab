#!/bin/bash
#SBATCH --job-name=atac_postprocess
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=18:00:00
#SBATCH --export=NONE

set -euo pipefail

runner="${5:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/bowtie2/postprocess_atac_PE.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: ATAC post-alignment runner not found: $runner" >&2; exit 2; }

source "$runner" "$1" "$2" "$3" "$4"
