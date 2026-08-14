#!/usr/bin/env bash
#SBATCH --job-name=codespring_scrna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE

set -euo pipefail

# `--export=NONE` gives a clean SLURM environment, which can omit PATH on
# some clusters.  The module initializer needs basic system tools (including
# ps) before it can construct the scientific runtime environment.
export PATH="${PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"

# Propagate the SLURM allocation to numerical libraries used by Seurat and
# Scanpy. Stage-specific submissions can override the fallback allocation
# above without editing this wrapper.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export VECLIB_MAXIMUM_THREADS="$OMP_NUM_THREADS"
export NUMEXPR_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

# SLURM batch shells do not always inherit the interactive module function.
# Initialize it when available so the runner sees the intended R/Python stack.
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

runner="$1"
engine="$2"
samples="$3"
out_dir="$4"
params="$5"
stage="${6:-all}"
scanpy_container="${7:-}"

mkdir -p "$out_dir"
# Keep large R/Python temporary files in project storage, not a shared node
# /tmp that may be full or too small for a real single-cell object.
export TMPDIR="$out_dir/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
export MPLCONFIGDIR="$TMPDIR/matplotlib"
export NUMBA_CACHE_DIR="$TMPDIR/numba"
export XDG_CACHE_HOME="$TMPDIR/cache"
mkdir -p "$TMPDIR"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME"
if [[ "$stage" == "all" || "$stage" == "annotate" ]]; then
  rm -f "$out_dir/_COMPLETE"
fi
# A successfully completed upstream marker must never make stale downstream
# outputs look current. Keep the files recoverable for inspection, but clear
# only the exact completion markers whose inputs are being superseded.
case "$stage" in
  inspect) rm -f "$out_dir/_STAGE_QC_COMPLETE" "$out_dir/_STAGE_PREPROCESS_COMPLETE" "$out_dir/_STAGE_CLUSTER_COMPLETE" "$out_dir/_STAGE_ANNOTATE_COMPLETE" "$out_dir/_STAGE_SCORE_COMPLETE" "$out_dir/_STAGE_DIFFERENTIAL_COMPLETE" "$out_dir/_STAGE_PATHWAY_COMPLETE" "$out_dir/_COMPLETE" ;;
  qc) rm -f "$out_dir/_STAGE_PREPROCESS_COMPLETE" "$out_dir/_STAGE_CLUSTER_COMPLETE" "$out_dir/_STAGE_ANNOTATE_COMPLETE" "$out_dir/_STAGE_SCORE_COMPLETE" "$out_dir/_STAGE_DIFFERENTIAL_COMPLETE" "$out_dir/_STAGE_PATHWAY_COMPLETE" "$out_dir/_COMPLETE" ;;
  preprocess) rm -f "$out_dir/_STAGE_CLUSTER_COMPLETE" "$out_dir/_STAGE_ANNOTATE_COMPLETE" "$out_dir/_STAGE_SCORE_COMPLETE" "$out_dir/_STAGE_DIFFERENTIAL_COMPLETE" "$out_dir/_STAGE_PATHWAY_COMPLETE" "$out_dir/_COMPLETE" ;;
  cluster) rm -f "$out_dir/_STAGE_ANNOTATE_COMPLETE" "$out_dir/_STAGE_SCORE_COMPLETE" "$out_dir/_STAGE_DIFFERENTIAL_COMPLETE" "$out_dir/_STAGE_PATHWAY_COMPLETE" "$out_dir/_COMPLETE" ;;
  annotate) rm -f "$out_dir/_STAGE_DIFFERENTIAL_COMPLETE" "$out_dir/_STAGE_PATHWAY_COMPLETE" ;;
  differential)
    rm -f "$out_dir/_STAGE_PATHWAY_COMPLETE"
    rm -f "$out_dir/tables/pseudobulk_differential_expression.tsv" "$out_dir/tables/pseudobulk_normalized_counts.tsv" "$out_dir/tables/cell_level_differential_expression.tsv"
    ;;
esac
date '+%Y-%m-%dT%H:%M:%S%z' > "$out_dir/_RUN_STARTED_${stage}"
trap 'rm -f "$out_dir/_RUN_STARTED_${stage}"' EXIT

if [[ ! -x "$runner" ]]; then
  echo "ERROR: scRNA-seq runner is missing or not executable: $runner" >&2
  exit 2
fi
"$runner" "$engine" "$samples" "$out_dir" "$params" "$stage" "$scanpy_container"
# Preserve job-local temporary files on failure for diagnosis, but remove them
# after a successful run so a completed large analysis does not leave a second
# copy of its intermediate data in the project results directory.
rm -rf "$TMPDIR" || echo "WARNING: completed scRNA workflow could not remove temporary directory: $TMPDIR" >&2
