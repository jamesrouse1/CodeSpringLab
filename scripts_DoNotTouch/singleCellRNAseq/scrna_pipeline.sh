#!/usr/bin/env bash
set -euo pipefail

# Arguments: engine samples.tsv output_dir params.tsv
engine="$1"
samples="$2"
out_dir="$3"
params="$4"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${engine,,}" in
  seurat)
    # Match CodeSpringLab's portable module setup before using the configured
    # Python interpreter for Scanpy/scvi.
    module load EBModules 2>/dev/null || true
    module load R/4.3.2-gfbf-2023a 2>/dev/null || true
    Rscript "$script_dir/scrna_pipeline_seurat.R" "$samples" "$out_dir" "$params"
    ;;
  scanpy)
    module load EBModules 2>/dev/null || true
    module load Anaconda3/2021.05 2>/dev/null || true
    python_bin="${CSL_SCRNA_SCANPY_PYTHON:-python3}"
    "$python_bin" "$script_dir/scrna_pipeline_scanpy.py" "$samples" "$out_dir" "$params"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac
