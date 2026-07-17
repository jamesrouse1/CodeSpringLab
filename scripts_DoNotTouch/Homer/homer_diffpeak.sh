#!/usr/bin/env bash
set -Eeuo pipefail

outdir="${1:?ERROR: output directory is required}"
input_tags="${2:?ERROR: reference tag directories are required}"
target_tags="${3:?ERROR: comparison tag directories are required}"
genome="${4:?ERROR: HOMER genome is required}"
reference="${5:?ERROR: reference condition is required}"
comparison="${6:?ERROR: comparison condition is required}"

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a
command -v getDifferentialPeaksReplicates.pl >/dev/null 2>&1 || { echo "ERROR: getDifferentialPeaksReplicates.pl was not found." >&2; exit 127; }
command -v annotatePeaks.pl >/dev/null 2>&1 || { echo "ERROR: annotatePeaks.pl was not found." >&2; exit 127; }

mkdir -p "$outdir"
prefix="DiffPeak_${comparison}_vs_${reference}(ref)"
result="${outdir}/${prefix}.txt"
annotation="${outdir}/${prefix}_annotated.txt"
result_tmp="${result}.tmp.$$"
annotation_tmp="${annotation}.tmp.$$"
cleanup() { rm -f "$result_tmp" "$annotation_tmp"; }
trap cleanup EXIT

getDifferentialPeaksReplicates.pl -t "$target_tags" -i "$input_tags" -DESeq2 -genome "$genome" -f 0.0001 -q 1 > "$result_tmp"
[[ -s "$result_tmp" ]] || { echo "ERROR: HOMER differential-peak output is empty." >&2; exit 1; }
annotatePeaks.pl "$result_tmp" "$genome" -raw > "$annotation_tmp"
[[ -s "$annotation_tmp" ]] || { echo "ERROR: HOMER differential-peak annotation is empty." >&2; exit 1; }
mv "$result_tmp" "$result"
mv "$annotation_tmp" "$annotation"

