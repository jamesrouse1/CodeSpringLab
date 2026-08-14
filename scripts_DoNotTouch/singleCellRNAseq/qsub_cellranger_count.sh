#!/usr/bin/env bash
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE

set -Eeuo pipefail

# A clean SLURM environment may omit the basic system PATH required by the
# cluster module initializer.
export PATH="${PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    if [[ -r "$module_init" ]]; then
      set +u
      # shellcheck disable=SC1090
      source "$module_init"
      set -u
      break
    fi
  done
fi
stage_parent="${7:-}"
runner="${8:-}"
# Backward compatibility for direct calls that used argument 7 for the runner
# before a separate staging filesystem became configurable.
if [[ -z "$runner" && -s "$stage_parent" ]]; then
  runner="$stage_parent"
  stage_parent="$5"
fi
[[ -s "$runner" ]] || { echo "ERROR: Cell Ranger runner was not found: $runner" >&2; exit 2; }
allocated_mem_mb="${SLURM_MEM_PER_NODE:-98304}"
[[ "$allocated_mem_mb" =~ ^[0-9]+$ ]] || allocated_mem_mb=98304
allocated_mem_gb=$((allocated_mem_mb / 1024))
(( allocated_mem_gb > 8 )) || allocated_mem_gb=64
export CELLRANGER_LOCALMEM_GB=$((allocated_mem_gb - 6))
export CELLRANGER_MIN_FREE_GB="${CELLRANGER_MIN_FREE_GB:-150}"
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$stage_parent"
