#!/usr/bin/env bash
set -Eeuo pipefail

sample="${1:?ERROR: ATAC sample is required}"
outdir="${5:?ERROR: ATAC MACS2 output directory is required}"
if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
mkdir -p "$outdir"
tmp_dir="${outdir}/.macs2_tmp_${SLURM_JOB_ID:-$$}"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir" TMP="$tmp_dir" TEMP="$tmp_dir"
cleanup() { rm -rf "$tmp_dir"; }
trap cleanup EXIT

module load EBModules
module load MACS2/2.2.9.1-foss-2022b
command -v macs2 >/dev/null 2>&1 || { echo "ERROR: macs2 was not found after loading its module." >&2; exit 127; }
## Alternative local installation:
## /grid/bsr/home/utama/.local/bin/macs2

run_log="${outdir}/${sample}_macs2.log"
summary="${outdir}/${sample}_macs2_summary.txt"
complete_marker="${outdir}/${sample}_macs2_complete.txt"
rm -f "$complete_marker" "$run_log"
macs_status=0
macs2 callpeak --nomodel \
	-t "${2}" -q "${7}" \
	-f BED -g "${3}" \
	--shift -100 --extsize 200 --bdg --call-summits \
	-n "${1}" --keep-dup all \
	--outdir "${5}" 2> "$run_log" || macs_status=$?
cat "$run_log" >&2
if ((macs_status != 0)); then
    echo "ERROR: MACS2 peak calling failed for ${1}." >&2
    exit "$macs_status"
fi

fatal_pattern='Traceback \(most recent call last\):|Exception ignored in:|KeyError:|ValueError:|TypeError:|OSError:|No space left on device|Killed|Segmentation fault'
if grep -Eq "$fatal_pattern" "$run_log"; then
    echo "ERROR: MACS2 reported an internal peak-calling exception for ${1}. See ${run_log}." >&2
    exit 1
fi

peak_file="${5}/${1}_peaks.narrowPeak"
if [[ ! -f "$peak_file" ]]; then
    echo "ERROR: MACS2 did not create the expected peak file (an empty file is valid when zero peaks pass): $peak_file" >&2
    exit 1
fi
if [[ ! -s "${outdir}/${sample}_peaks.xls" ]]; then
    echo "ERROR: MACS2 did not create a non-empty peaks table for ${sample}." >&2
    exit 1
fi

#module load Anaconda3/2023.03-1
#conda activate deeptools

module load EBModules
module load deepTools/3.5.2-foss-2022a

### TSS was used for workshop instead of center
    
tss_bed="${6}"
signal_bw="${9}/${1}/${1}Aligned.sortedByCoord_removeDup.out.bw"
tss_matrix="${5}/${1}_TSS.gz"
if [[ -s "$tss_bed" && -s "$signal_bw" ]]; then
    if computeMatrix reference-point -p 4 \
        --referencePoint TSS \
        -b 1000 -a 1000 \
        -R "$tss_bed" \
        -S "$signal_bw" \
        --skipZeros \
        -o "$tss_matrix" \
        --outFileSortedRegions "${5}/${1}_genes_TSS.bed"; then
        plotHeatmap -m "$tss_matrix" -out "${5}/${1}_heatmap_TSS.png" || \
            echo "WARNING: TSS heatmap plotting failed for ${1}; MACS2 peaks are still valid." >&2
    else
        echo "WARNING: TSS matrix generation failed for ${1}; MACS2 peaks are still valid." >&2
    fi
else
    echo "WARNING: Skipping the optional TSS heatmap for ${1}; annotation BED or signal bigWig is missing." >&2
fi

#computeMatrix reference-point -p 4 \
#    --referencePoint center \
#    -b 1000 -a 1000 \
#    -R ${6} \
#    -S ${5}/${1}Aligned.sortedByCoord_removeDup.out.bw \
#    --skipZeros \
#    -o ${5}/${1}_peakCenter.gz \
#    --outFileSortedRegions ${5}/${1}_genes_peakCenter.bed
     
     
#plotHeatmap -m ${5}/${1}_peakCenter.gz \
#    -out ${5}/${1}_heatmap_peakCenter.png \

#conda deactivate

###### Homer Annotate #########

module load EBModules
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a

export PATH=$PATH:/grid/bsr/data/data/utama/tools/homer/bin/

annotation="${5}/${1}_peaks_annotated.txt"
annotation_tmp="${annotation}.tmp.$$"
rm -f "$annotation_tmp"
if [[ -s "$peak_file" ]]; then
    command -v annotatePeaks.pl >/dev/null 2>&1 || {
        echo "ERROR: HOMER annotatePeaks.pl is unavailable; ${1} MACS2 annotation was not completed." >&2
        exit 127
    }
    if ! annotatePeaks.pl "${5}/${1}_peaks.xls" "${8}" > "$annotation_tmp"; then
        echo "ERROR: HOMER peak annotation failed for ${1}." >&2
        rm -f "$annotation_tmp"
        exit 1
    fi
    if [[ ! -s "$annotation_tmp" ]] || ! head -n 1 "$annotation_tmp" | grep -q '^PeakID'; then
        echo "ERROR: HOMER produced an empty or invalid MACS2 annotation for ${1}." >&2
        rm -f "$annotation_tmp"
        exit 1
    fi
    mv -f "$annotation_tmp" "$annotation"
else
    printf 'Status\tNo MACS2 peaks passed the selected threshold\nSample\t%s\n' "${1}" > "$annotation"
fi

peak_count="$(wc -l < "$peak_file" | tr -d '[:space:]')"
caller_version="$(macs2 --version 2>&1 | head -n 1)"
{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tcomplete\n'
  printf 'caller\tmacs2\n'
  printf 'caller_version\t%s\n' "$caller_version"
  printf 'input_bed\t%s\n' "$2"
  printf 'genome_size\t%s\n' "$3"
  printf 'qvalue\t%s\n' "$7"
  printf 'peak_file\t%s\n' "$peak_file"
  printf 'peak_count\t%s\n' "$peak_count"
  printf 'run_log\t%s\n' "$run_log"
} > "$summary"
marker_tmp="${complete_marker}.tmp.$$"
{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tcomplete\n'
  printf 'peak_file\t%s\n' "$peak_file"
  printf 'peak_count\t%s\n' "$peak_count"
  printf 'summary\t%s\n' "$summary"
} > "$marker_tmp"
mv "$marker_tmp" "$complete_marker"

exit 0
