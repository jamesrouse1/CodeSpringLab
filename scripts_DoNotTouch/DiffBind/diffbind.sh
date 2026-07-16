
set -euo pipefail

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }

outpath="$2"
mkdir -p "$outpath"
run_started="${outpath}/_RUN_STARTED"
complete_marker="${outpath}/_COMPLETE"
rm -f "$complete_marker"
printf 'status\trunning\nstarted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$run_started"
tmp_dir="${outpath}/.diffbind_tmp_${SLURM_JOB_ID:-$$}"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir" TMP="$tmp_dir" TEMP="$tmp_dir"
cleanup() { rm -rf "$tmp_dir"; }
trap cleanup EXIT

module load EBModules
module load R/4.3.2-gfbf-2023a

Rscript "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"

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
  rm -f "$run_started"
  printf 'status\tcomplete\ncompleted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$complete_marker"
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
rm -f "$run_started"
printf 'status\tcomplete\ncompleted_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "$complete_marker"
