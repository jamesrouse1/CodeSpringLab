#!/usr/bin/env bash
#SBATCH --job-name=cutrun_peakqc
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=12:00:00

source ../scripts_DoNotTouch/CUTRUN/cutrun_peak_qc.sh "$1" "$2" "$3" "$4"
