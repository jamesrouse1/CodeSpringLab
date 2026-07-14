#!/usr/bin/env bash
#SBATCH --job-name=cutrun_seacr
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=12:00:00

source ../scripts_DoNotTouch/SEACR/seacr_cutrun.sh "$1" "$2" "$3" "$4" "$5" "$6" "${7:-none}"
