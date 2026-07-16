#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
runner="${script_dir}/bowtie2_PE.sh"
[[ -s "$runner" ]] || { echo "ERROR: shared paired-end Bowtie2 runner is missing: $runner" >&2; exit 2; }

export PIPELINE_ASSAY_LABEL="ChIP-seq"
export BAMCOVERAGE_IGNORE_CHROMS="chrM"
exec bash "$runner" "$@"
