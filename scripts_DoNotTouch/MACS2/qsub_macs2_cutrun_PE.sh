#!/usr/bin/env bash
#SBATCH --job-name=cutrun_macs2
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --time=12:00:00

source ../scripts_DoNotTouch/MACS2/macs2_cutrun_PE.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
