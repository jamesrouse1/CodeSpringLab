#!/bin/bash
#SBATCH --job-name=atac_bigwig
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

set -euo pipefail

runner="${3:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/bowtie2/atac_bigwig.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: ATAC bigWig runner not found: $runner" >&2; exit 2; }

source "$runner" "$1" "$2"
