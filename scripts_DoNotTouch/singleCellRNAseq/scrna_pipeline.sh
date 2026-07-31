#!/usr/bin/env bash
set -euo pipefail

# Arguments: engine samples.tsv output_dir params.tsv
engine="$1"
samples="$2"
out_dir="$3"
params="$4"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

runtime_executable="$(awk -F $'\t' '$1 == "runtime_executable" { print $2; exit }' "$params" 2>/dev/null | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"

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
    if [[ -z "$runtime_executable" ]]; then
      module load EBModules 2>/dev/null || true
      module load R/4.3.2-gfbf-2023a 2>/dev/null || true
      runtime_executable="$(command -v Rscript || true)"
    fi
    require_executable "$runtime_executable" "Seurat Rscript"
    "$runtime_executable" -e 'for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork")) if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ", pkg)'
    "$runtime_executable" "$script_dir/scrna_pipeline_seurat.R" "$samples" "$out_dir" "$params"
    ;;
  scanpy)
    if [[ -z "$runtime_executable" ]]; then
      module load EBModules 2>/dev/null || true
      module load Anaconda3/2021.05 2>/dev/null || true
      runtime_executable="$(command -v python3 || true)"
    fi
    require_executable "$runtime_executable" "Scanpy Python"
    "$runtime_executable" -c 'import anndata, igraph, leidenalg, numpy, pandas, scanpy, scipy; print("Scanpy runtime ready")'
    "$runtime_executable" "$script_dir/scrna_pipeline_scanpy.py" "$samples" "$out_dir" "$params"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac
