#!/usr/bin/env bash
set -euo pipefail

# Arguments: engine samples.tsv output_dir params.tsv [stage]
engine="$1"
samples="$2"
out_dir="$3"
params="$4"
stage="${5:-all}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# This script is executed as a new shell by the SLURM wrapper. Shell functions
# such as `module` are not guaranteed to cross that boundary, even when the
# wrapper already initialized modules. Initialize them again here so a normal
# server run cannot fail with `module: command not found`.
if ! command -v module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$module_init" ]]; then
      # shellcheck disable=SC1090
      source "$module_init"
      break
    fi
  done
fi

require_executable() {
  local executable="$1"
  local label="$2"
  if [[ -z "$executable" || ! -x "$executable" ]]; then
    echo "ERROR: The requested $label runtime is not executable: $executable" >&2
    exit 2
  fi
}

engine="$(printf '%s' "$engine" | tr '[:upper:]' '[:lower:]')"
case "$engine" in
  seurat)
    # Match CodeSpringLab's portable module setup before using R.
    module load EBModules 2>/dev/null || true
    module load R/4.3.2-gfbf-2023a 2>/dev/null || true
    runtime_executable="$(command -v Rscript || true)"
    require_executable "$runtime_executable" "Seurat Rscript"
    "$runtime_executable" -e 'for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork")) if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ", pkg)'
    "$runtime_executable" "$script_dir/scrna_pipeline_seurat.R" "$samples" "$out_dir" "$params" "$stage"
    ;;
  scanpy)
    module load EBModules 2>/dev/null || true
    module load Anaconda3/2021.05 2>/dev/null || true
    runtime_executable="$(command -v python3 || true)"
    require_executable "$runtime_executable" "Scanpy Python"
    "$runtime_executable" -c 'import anndata, igraph, leidenalg, numpy, pandas, scanpy, scipy; print("Scanpy runtime ready")'
    "$runtime_executable" "$script_dir/scrna_pipeline_scanpy.py" "$samples" "$out_dir" "$params" "$stage"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac
