#!/usr/bin/env bash
set -euo pipefail

r_script="$1"
sample_sheet="$2"
out_dir="$3"
reference_condition="$4"
min_replicates="$5"
genome="$6"
blacklist="${7:-none}"
comparison_condition="${8:-}"
cell_type="${9:-}"
mark="${10:-}"

module load EBModules
module load R/4.3.2-gfbf-2023a

Rscript "$r_script" "$sample_sheet" "$out_dir" "$reference_condition" "$min_replicates" "$genome" "$blacklist" "$comparison_condition" "$cell_type" "$mark"
