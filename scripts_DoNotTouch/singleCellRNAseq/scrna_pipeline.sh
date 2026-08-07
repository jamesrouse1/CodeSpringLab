#!/usr/bin/env bash
set -euo pipefail

# Arguments: engine samples.tsv output_dir params.tsv [stage] [scanpy_container.sif]
engine="$1"
samples="$2"
out_dir="$3"
params="$4"
stage="${5:-all}"
scanpy_container_path="${6:-${CSL_SCANPY_SIF:-}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# This script is executed as a new shell by the SLURM wrapper. Shell functions
# such as `module` are not guaranteed to cross that boundary, even when the
# wrapper already initialized modules. Initialize them again here so a normal
# server run cannot fail with `module: command not found`.
if ! command -v module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$module_init" ]]; then
      # shellcheck disable=SC1090
      # The site module initializer references optional variables that are
      # legitimately absent in a clean SLURM environment.  Temporarily relax
      # nounset only while it defines the `module` shell function.
      set +u
      source "$module_init"
      set -u
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

scanpy_container() {
  # H5AD jobs run inside one immutable, tested SIF image. This avoids package
  # solving and makes results reproducible across users and projects.
  local default_container="$script_dir/containers/codespring-scanpy_1.0.0.sif"
  local container="${scanpy_container_path:-$default_container}"

  module load singularity/3.6.3 >/dev/null 2>&1 || true
  local singularity_executable="$(command -v singularity || true)"
  if [[ -z "$singularity_executable" ]]; then
    echo "ERROR: Singularity could not be loaded. This cluster requires singularity/3.6.3 for CodeSpringLab Scanpy jobs." >&2
    return 2
  fi
  if [[ ! -r "$container" ]]; then
    echo "ERROR: The shared Scanpy container is unavailable: $container" >&2
    echo "A CodeSpringLab maintainer must install the versioned SIF image once; individual users do not need to create Python environments." >&2
    return 2
  fi
  printf '%s\t%s\n' "$singularity_executable" "$container"
}

scanpy_bind_args() {
  # Bind project outputs, the manifest, the runner, and every parent input
  # directory explicitly. This works whether project data are in $HOME, /grid,
  # or another filesystem that Singularity does not auto-bind on a given node.
  local -a directories=("$out_dir" "$(dirname "$samples")" "$script_dir")
  local input_column
  input_column="$(awk -F '\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="input_path") {print i; exit}}' "$samples")"
  if [[ -n "$input_column" ]]; then
    while IFS= read -r input_path; do
      [[ -n "$input_path" ]] && directories+=("$(dirname "$input_path")")
    done < <(awk -F '\t' -v col="$input_column" 'NR>1 && $col != "" {print $col}' "$samples")
  fi
  local seen=$'\n'
  local directory
  for directory in "${directories[@]}"; do
    [[ -d "$directory" && "$seen" != *$'\n'"$directory"$'\n'* ]] || continue
    seen+="$directory"$'\n'
    printf '%s\n' "--bind=$directory"
  done
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
    runtime_info="$(scanpy_container)"
    singularity_executable="${runtime_info%%$'\t'*}"
    container="${runtime_info#*$'\t'}"
    if [[ -z "$singularity_executable" || -z "$container" || "$runtime_info" == "$singularity_executable" ]]; then
      echo "ERROR: Could not resolve the shared Scanpy container runtime." >&2
      exit 2
    fi
    bind_args=()
    while IFS= read -r bind_arg; do
      bind_args+=("$bind_arg")
    done < <(scanpy_bind_args)
    echo "Scanpy container ready: $(basename "$container")"
    SINGULARITYENV_MPLBACKEND=Agg \
      SINGULARITYENV_TMPDIR="${TMPDIR:-$out_dir/tmp}" \
      SINGULARITYENV_MPLCONFIGDIR="${MPLCONFIGDIR:-$out_dir/tmp/matplotlib}" \
      SINGULARITYENV_NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$out_dir/tmp/numba}" \
      SINGULARITYENV_XDG_CACHE_HOME="${XDG_CACHE_HOME:-$out_dir/tmp/cache}" \
      "$singularity_executable" exec --cleanenv "${bind_args[@]}" "$container" \
      python "$script_dir/scrna_pipeline_scanpy.py" "$samples" "$out_dir" "$params" "$stage"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac
