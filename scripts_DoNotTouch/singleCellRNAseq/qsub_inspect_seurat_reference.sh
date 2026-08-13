#!/usr/bin/env bash
#SBATCH --job-name=codespring_reference_labels
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --export=NONE

set -euo pipefail

reference_path="$1"
output_path="$2"
runner="$3"

export PATH="${PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"
if ! command -v module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$module_init" ]]; then
      set +u
      # shellcheck disable=SC1090
      source "$module_init"
      set -u
      break
    fi
  done
fi

module purge >/dev/null 2>&1 || true
module load EB5Modules
module load "${CSL_SEURAT_MODULE:-Seurat/5.4.0-foss-2024a-R-4.4.2}"
export R_LIBS_USER="${CSL_R_LIBS_USER:-$(dirname "$output_path")/.codespring_unused_user_library}"
export R_ENVIRON_USER="${CSL_R_ENVIRON_USER:-/dev/null}"
export R_PROFILE_USER="${CSL_R_PROFILE_USER:-/dev/null}"
mkdir -p "$(dirname "$output_path")"
printf 'Reference inspection started: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
printf 'Reference: %s\n' "$reference_path"
Rscript "$runner" "$reference_path" "$output_path"
printf 'Reference inspection finished: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
