#!/usr/bin/env bash
#SBATCH --job-name=cutrun_bowtie2
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

# Keep the same argument layout as paired-end. $4 is an unused R2 placeholder.
source ../scripts_DoNotTouch/bowtie2/bowtie2_cutrun_SE.sh "$1" "$2" "$3" "$5" "$6" "$7" "$8" "$9" "${10}"
