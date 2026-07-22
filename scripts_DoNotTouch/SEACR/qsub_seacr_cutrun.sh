#!/usr/bin/env bash
#SBATCH --job-name=cutrun_seacr
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=12:00:00

seacr_runner="${8:-}"
if [[ -z "$seacr_runner" || ! -s "$seacr_runner" ]]; then
  seacr_runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/SEACR/seacr_cutrun.sh"
fi
if [[ ! -s "$seacr_runner" ]]; then
  echo "ERROR: SEACR runner was not found at: $seacr_runner" >&2
  exit 2
fi

## --export=NONE can omit the cluster module initialization on compute nodes.
if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    if [[ -s "$module_init" ]]; then
      # shellcheck disable=SC1090
      source "$module_init"
      break
    fi
  done
fi

source "$seacr_runner" "$1" "$2" "$3" "$4" "$5" "$6" "${7:-none}"
