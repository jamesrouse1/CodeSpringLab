
set -euo pipefail

module load EBModules
module load R/4.3.2-gfbf-2023a

Rscript "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"

outpath="$2"
refcond="$5"
compared="$6"
genome="$7"
prefix="DifferentialPeaks_${compared}_vs_${refcond}_ref"
bed_file="${outpath}/${prefix}.with_stats.bed"

if [[ ! -e "$bed_file" ]]; then
  echo "ERROR: DiffBind did not create the differential-peak BED file: ${bed_file}" >&2
  exit 1
fi

if [[ ! -s "$bed_file" ]]; then
  {
    printf 'Status\tNo differential peaks passed the DiffBind reporting threshold\n'
    printf 'Comparison\t%s vs %s\n' "$compared" "$refcond"
  } > "${outpath}/${prefix}_annotated_with_stats.txt"
  echo "No significant differential peaks; wrote ${outpath}/${prefix}_annotated_with_stats.txt"
  exit 0
fi

case "$genome" in
  mouse) homer_genome="mm39" ;;
  human) homer_genome="hg38" ;;
  *) echo "ERROR: unsupported Homer genome: ${genome}" >&2; exit 2 ;;
esac

module unload R/4.3.2-gfbf-2023a >/dev/null 2>&1 || true
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a
export PATH="$PATH:/grid/bsr/data/data/utama/tools/homer/bin"

annotatePeaks.pl "$bed_file" "$homer_genome" \
  > "${outpath}/${prefix}_annotated_with_stats.txt"

if [[ ! -s "${outpath}/${prefix}_annotated_with_stats.txt" ]]; then
  echo "ERROR: Homer annotation output is empty for ${prefix}" >&2
  exit 1
fi

echo "Wrote ${outpath}/${prefix}_annotated_with_stats.txt"
