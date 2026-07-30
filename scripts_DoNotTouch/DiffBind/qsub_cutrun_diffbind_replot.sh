#!/usr/bin/env bash
#SBATCH --job-name=cutrun_diffbind_pca
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --time=04:00:00

set -euo pipefail

r_script="$1"
object_path="$2"
out_dir="$3"

module load EBModules
module load R/4.3.2-gfbf-2023a

Rscript "$r_script" "$object_path" "$out_dir"
