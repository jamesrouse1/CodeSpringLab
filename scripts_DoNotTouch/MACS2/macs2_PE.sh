

module load EBModules
module load MACS2/2.2.9.1-foss-2022b
## Alternative local installation:
## /grid/bsr/home/utama/.local/bin/macs2

macs2 callpeak --nomodel \
	-t "${2}" -q "${7}" \
	-f BED -g "${3}" \
	--shift -100 --extsize 200 --bdg --call-summits \
	-n "${1}" --keep-dup all \
	--outdir "${5}" || {
    echo "ERROR: MACS2 peak calling failed for ${1}." >&2
    return 1
}

peak_file="${5}/${1}_peaks.narrowPeak"
if [[ ! -s "$peak_file" ]]; then
    echo "ERROR: MACS2 did not create the expected peak file: $peak_file" >&2
    return 1
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

if command -v annotatePeaks.pl >/dev/null 2>&1; then
    annotatePeaks.pl "${5}/${1}_peaks.xls" "${8}" > "${5}/${1}_peaks_annotated.txt" || \
        echo "WARNING: HOMER peak annotation failed for ${1}; MACS2 peaks are still valid." >&2
else
    echo "WARNING: Skipping optional HOMER annotation for ${1}; annotatePeaks.pl is unavailable." >&2
fi

return 0
