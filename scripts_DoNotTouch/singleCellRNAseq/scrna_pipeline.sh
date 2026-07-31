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

scanpy_runtime() {
  # Scanpy is not part of many cluster-wide Anaconda modules. Keep one
  # managed, per-user environment outside either repository so every H5AD run
  # has the same complete runtime and no username is hard-coded.
  local allow_setup="${1:-false}"
  local runtime_root="${CSL_WEB_HOME:-$HOME/.codespringweb}/runtimes"
  local environment="${CSL_SCANPY_ENV:-$runtime_root/scanpy}"
  local python="$environment/bin/python"
  local conda_executable=""

  if [[ -x "$python" ]] && "$python" -c 'import anndata, harmonypy, igraph, leidenalg, numpy, pandas, scanpy, scipy, scrublet' >/dev/null 2>&1; then
    printf '%s\n' "$python"
    return 0
  fi

  if [[ "$allow_setup" != "true" ]]; then
    echo "ERROR: The managed Scanpy environment is not ready at $environment." >&2
    echo "In CodeSpringApp, use the 'Set up Scanpy environment' step for this H5AD project, then rerun the requested stage." >&2
    return 2
  fi

  conda_executable="$(command -v conda || true)"
  if [[ -z "$conda_executable" && -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    conda_executable="$CONDA_EXE"
  fi
  if [[ -z "$conda_executable" ]]; then
    echo "ERROR: Scanpy is not available and Conda could not be found to create the managed runtime." >&2
    echo "Load a supported Anaconda module, then rerun this stage." >&2
    return 2
  fi

  echo "INFO: Creating the one-time managed Scanpy runtime at $environment" >&2
  echo "INFO: This first H5AD job can take several minutes; later Scanpy jobs reuse it." >&2
  mkdir -p "$runtime_root"
  # This cluster's Conda mirror provides Scanpy and its compiled scientific
  # stack, but not Scrublet/HarmonyPy. Install the core reproducibly with
  # Conda, then use the environment's own pip for only those pure-Python
  # extensions. --force also repairs a partial environment left by a failed
  # first setup attempt at this dedicated CodeSpring-owned location.
  "$conda_executable" create -y --force -p "$environment" -c conda-forge \
    python=3.11 pip scanpy anndata python-igraph leidenalg >&2
  python="$environment/bin/python"
  if [[ -x "$python" ]]; then
    echo "INFO: Installing Scrublet and HarmonyPy with the managed environment's pip." >&2
    "$python" -m pip install --no-cache-dir scrublet harmonypy >&2
  fi
  if [[ ! -x "$python" ]] || ! "$python" -c 'import anndata, harmonypy, igraph, leidenalg, numpy, pandas, scanpy, scipy, scrublet' >/dev/null 2>&1; then
    echo "ERROR: The managed Scanpy runtime could not be initialized at $environment." >&2
    echo "Check the job log for the Conda error, then rerun this stage; no input data were modified." >&2
    return 2
  fi
  printf '%s\n' "$python"
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
    if [[ "$stage" == "setup-runtime" ]]; then
      runtime_executable="$(scanpy_runtime true)"
    else
      runtime_executable="$(scanpy_runtime)"
    fi
    require_executable "$runtime_executable" "managed Scanpy Python"
    echo "Scanpy runtime ready: $runtime_executable"
    if [[ "$stage" == "setup-runtime" ]]; then
      mkdir -p "$out_dir"
      printf 'ready\n' > "$out_dir/_SCANPY_RUNTIME_READY"
      echo "Managed Scanpy runtime setup completed."
      exit 0
    fi
    "$runtime_executable" "$script_dir/scrna_pipeline_scanpy.py" "$samples" "$out_dir" "$params" "$stage"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac
