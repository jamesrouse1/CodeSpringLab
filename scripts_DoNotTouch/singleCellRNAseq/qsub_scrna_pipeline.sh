#!/usr/bin/env bash
#SBATCH --job-name=codespring_scrna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --export=NONE

set -euo pipefail

# SLURM batch shells do not always inherit the interactive module function.
# Initialize it when available so the runner sees the intended R/Python stack.
if ! command -v module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$module_init" ]]; then
      # shellcheck disable=SC1090
      source "$module_init"
      break
    fi
  done
fi

runner="$1"
engine="$2"
samples="$3"
out_dir="$4"
params="$5"

mkdir -p "$out_dir"
# Keep large R/Python temporary files in project storage, not a shared node
# /tmp that may be full or too small for a real single-cell object.
export TMPDIR="$out_dir/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"
rm -f "$out_dir/_COMPLETE"
date -Is > "$out_dir/_RUN_STARTED"
trap 'rm -f "$out_dir/_RUN_STARTED"' EXIT

if [[ ! -x "$runner" ]]; then
  echo "ERROR: scRNA-seq runner is missing or not executable: $runner" >&2
  exit 2
fi
"$runner" "$engine" "$samples" "$out_dir" "$params"
