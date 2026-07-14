#!/usr/bin/env bash
#SBATCH --job-name=cutrun_seacr
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=12:00:00

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/seacr_cutrun.sh" "$1" "$2" "$3" "$4" "$5" "$6" "${7:-none}"
