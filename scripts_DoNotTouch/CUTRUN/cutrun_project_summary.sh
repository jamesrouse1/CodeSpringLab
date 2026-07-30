#!/usr/bin/env bash
set -euo pipefail

r_script="$1"
data_dir="$2"
out_dir="$3"
kind="$4"

module load EBModules
module load R/4.3.2-gfbf-2023a

Rscript "$r_script" "$data_dir" "$out_dir" "$kind"
