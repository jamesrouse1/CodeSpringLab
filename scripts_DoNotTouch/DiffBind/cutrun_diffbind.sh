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
minimum_peaks_per_sample="${11:-1}"
peak_source="${12:-legacy}"

module load EBModules
module load R/4.3.2-gfbf-2023a

# Keep gene annotation optional: the R step records a readable status file if
# HOMER is unavailable, but a successful DiffBind analysis is never discarded.
export PATH="$PATH:/grid/bsr/data/data/utama/tools/homer/bin"

Rscript "$r_script" "$sample_sheet" "$out_dir" "$reference_condition" "$min_replicates" "$genome" "$blacklist" "$comparison_condition" "$cell_type" "$mark" "$minimum_peaks_per_sample" "$peak_source"
