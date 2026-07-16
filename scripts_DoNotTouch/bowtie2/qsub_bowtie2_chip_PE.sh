#!/usr/bin/env bash
#SBATCH --job-name=chip_bowtie2
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail
runner="${8:-}"
[[ -s "$runner" ]] || { echo "ERROR: ChIP-seq paired-end Bowtie2 runner not found: ${runner:-<missing>}" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7"
