#!/usr/bin/env bash
#SBATCH --job-name=csl_scrna_query
#SBATCH --partition=cpuq
#SBATCH --qos=cpu_snice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

# Run read-only queries against a completed single-cell object on a compute
# node. Loading a multi-gigabyte Seurat/H5AD object must never happen inside
# the long-lived CodeSpringApp web process on a login node.
#
# Arguments:
#   engine task object_path output_path value scanpy_container_path

engine="${1:?engine is required}"
task="${2:?task is required}"
object_path="${3:?object path is required}"
output_path="${4:?output path is required}"
value="${5:-}"
scanpy_container_path="${6:-${CSL_SCANPY_SIF:-}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -r "$object_path" ]] || { echo "ERROR: Object is not readable: $object_path" >&2; exit 2; }
mkdir -p "$(dirname "$output_path")"

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

case "${engine,,}" in
  seurat)
    module purge >/dev/null 2>&1 || true
    module load EB5Modules >/dev/null 2>&1
    module load "${CSL_SEURAT_MODULE:-Seurat/5.4.0-foss-2024a-R-4.4.2}" >/dev/null 2>&1
    export R_LIBS_USER="${CSL_R_LIBS_USER:-$(dirname "$output_path")/.codespring_unused_user_library}"
    export R_ENVIRON_USER="${CSL_R_ENVIRON_USER:-/dev/null}"
    export R_PROFILE_USER="${CSL_R_PROFILE_USER:-/dev/null}"
    if [[ "$task" == "gene" ]]; then
      [[ -n "$value" ]] || { echo "ERROR: A gene is required." >&2; exit 2; }
      Rscript "$script_dir/extract_seurat_marker_expression.R" "$object_path" "$value" "$output_path"
    elif [[ "$task" == "genes" ]]; then
      Rscript "$script_dir/list_seurat_genes.R" "$object_path" "$output_path"
    else
      echo "ERROR: Unsupported object-query task: $task" >&2
      exit 2
    fi
    ;;
  scanpy)
    module load singularity/3.6.3 >/dev/null 2>&1 || true
    singularity_executable="$(command -v singularity || true)"
    [[ -x "$singularity_executable" ]] || { echo "ERROR: Singularity is unavailable." >&2; exit 2; }
    [[ -r "$scanpy_container_path" ]] || { echo "ERROR: Scanpy container is unavailable: $scanpy_container_path" >&2; exit 2; }
    helper="$script_dir/extract_h5ad_marker_expression.py"
    helper_args=("$object_path" "$value" "$output_path")
    if [[ "$task" == "genes" ]]; then
      helper="$script_dir/list_h5ad_genes.py"
      helper_args=("$object_path" "$output_path")
    elif [[ "$task" != "gene" ]]; then
      echo "ERROR: Unsupported object-query task: $task" >&2
      exit 2
    fi
    bind_paths="$(dirname "$object_path"),$(dirname "$output_path"),$script_dir"
    SINGULARITYENV_TMPDIR="${TMPDIR:-$(dirname "$output_path")}" \
      SINGULARITYENV_NUMBA_CACHE_DIR="${TMPDIR:-$(dirname "$output_path")}/numba" \
      "$singularity_executable" exec --cleanenv --bind "$bind_paths" "$scanpy_container_path" \
      python "$helper" "${helper_args[@]}"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac

[[ -s "$output_path" ]] || { echo "ERROR: Object query did not create a non-empty output." >&2; exit 3; }
echo "Completed $task query for $(basename "$object_path")."
