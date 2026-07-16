#!/usr/bin/env bash
#SBATCH --job-name=cutrun_bowtie2
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

source ../scripts_DoNotTouch/bowtie2/bowtie2_cutrun_PE.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11:-CPM}" "${12:-none}" "${13:-spikein}" "${14:-1000}" "${15:-full}"
