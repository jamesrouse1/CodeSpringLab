#!/usr/bin/env bash
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
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
runner="${7:-}"
[[ -s "$runner" ]] || { echo "ERROR: Cell Ranger runner was not found: $runner" >&2; exit 2; }
export CELLRANGER_LOCALMEM_GB=60
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6"
