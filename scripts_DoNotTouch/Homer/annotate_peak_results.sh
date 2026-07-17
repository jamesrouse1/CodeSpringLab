#!/usr/bin/env bash
set -Eeuo pipefail

data_dir="${1:?Project data directory is required}"
genome="${2:?HOMER genome name is required}"
gtf="${3:-}"

if [[ ! -d "$data_dir" ]]; then
  echo "ERROR: project data directory does not exist: $data_dir" >&2
  exit 1
fi
if [[ -n "$gtf" && ! -s "$gtf" ]]; then
  echo "ERROR: annotation GTF is missing or empty: $gtf" >&2
  exit 1
fi

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    if [[ -s "$module_init" ]]; then
      # shellcheck disable=SC1090
      source "$module_init"
      break
    fi
  done
fi
if type module >/dev/null 2>&1; then
  module load EBModules >/dev/null 2>&1 || true
  module load Anaconda3/2021.05 >/dev/null 2>&1 || true
  module load R/4.1.2-foss-2021a >/dev/null 2>&1 || true
fi

if [[ -n "${CSL_HOMER_BIN:-}" ]]; then
  export PATH="${CSL_HOMER_BIN}:$PATH"
fi
if ! command -v annotatePeaks.pl >/dev/null 2>&1 && [[ -x /grid/bsr/data/data/utama/tools/homer/bin/annotatePeaks.pl ]]; then
  export PATH="/grid/bsr/data/data/utama/tools/homer/bin:$PATH"
fi
if ! command -v annotatePeaks.pl >/dev/null 2>&1; then
  echo "ERROR: annotatePeaks.pl was not found. Set CSL_HOMER_BIN or install HOMER in the shared tools location." >&2
  exit 127
fi

annotation_root="$data_dir/peak_annotation"
summary="$annotation_root/peak_annotation_summary.tsv"
complete_marker="$annotation_root/_COMPLETE"
running_marker="$annotation_root/_RUN_STARTED"
tmp_root="$annotation_root/.tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "$annotation_root"
rm -rf "$tmp_root"
mkdir -p "$tmp_root"
rm -f "$complete_marker"
{
  printf 'status\trunning\n'
  printf 'job_id\t%s\n' "${SLURM_JOB_ID:-local}"
  printf 'started_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
} > "$running_marker"
cleanup() { rm -rf "$tmp_root"; }
trap cleanup EXIT

summary_tmp="$tmp_root/peak_annotation_summary.tsv"
printf 'result_type\tsample_or_comparison\tpeak_count\tsource_peak_file\tannotated_file\tstatus\n' > "$summary_tmp"
annotated_count=0
total_peaks=0

safe_label() {
  printf '%s' "$1" | tr -cs '[:alnum:]_.-' '_'
}

run_annotation() {
  local result_type="$1" label="$2" source="$3" prepared="$4" destination="$5"
  local peak_count output_tmp
  peak_count="$(wc -l < "$prepared" | tr -d '[:space:]')"
  if [[ "$peak_count" == "0" ]]; then
    echo "ERROR: no valid genomic intervals remained after preparing: $source" >&2
    return 1
  fi
  output_tmp="$tmp_root/$(safe_label "${result_type}_${label}").annotated.txt"
  if [[ -n "$gtf" ]]; then
    annotatePeaks.pl "$prepared" "$genome" -gtf "$gtf" > "$output_tmp"
  else
    annotatePeaks.pl "$prepared" "$genome" > "$output_tmp"
  fi
  if [[ ! -s "$output_tmp" ]]; then
    echo "ERROR: HOMER produced an empty annotation for: $source" >&2
    return 1
  fi
  mkdir -p "$(dirname "$destination")"
  mv -f "$output_tmp" "$destination"
  printf '%s\t%s\t%s\t%s\t%s\tcomplete\n' \
    "$result_type" "$label" "$peak_count" "$source" "$destination" >> "$summary_tmp"
  annotated_count=$((annotated_count + 1))
  total_peaks=$((total_peaks + peak_count))
}

prepare_bed_peak() {
  local source="$1" prepared="$2" format="$3"
  awk -v format="$format" 'BEGIN {OFS="\t"}
    NF >= 3 && $1 !~ /^#/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {
      name=(NF >= 4 && $4 != "" && $4 != ".") ? $4 : "peak_" NR;
      score=(NF >= 5 && $5 != "") ? $5 : 0;
      strand=(NF >= 6 && ($6 == "+" || $6 == "-")) ? $6 : "+";
      id=name "|score=" score;
      if (format == "narrow" && NF >= 10) id=id "|signalValue=" $7 "|pValue=" $8 "|qValue=" $9 "|summitOffset=" $10;
      else if (format == "broad" && NF >= 9) id=id "|signalValue=" $7 "|pValue=" $8 "|qValue=" $9;
      print $1, $2, $3, id, score, strand;
    }' "$source" > "$prepared"
}

prepare_diffbind_tsv() {
  local source="$1" prepared="$2"
  awk 'BEGIN {FS=OFS="\t"}
    NR == 1 {
      for (i=1; i<=NF; i++) {
        header[i]=$i;
        key=tolower($i); gsub(/[^a-z0-9]/, "", key);
        if (key == "seqnames" || key == "chr" || key == "chrom" || key == "chromosome") chr=i;
        else if (key == "start") start=i;
        else if (key == "end") end=i;
        else if (key == "fold" || key == "log2foldchange" || key == "log2fc") fold=i;
      }
      if (!chr || !start || !end) exit 42;
      next;
    }
    $start ~ /^[0-9]+$/ && $end ~ /^[0-9]+$/ && $end >= $start {
      bed_start=$start-1; if (bed_start < 0) bed_start=0;
      id=$chr ":" $start "-" $end;
      for (i=1; i<=NF; i++) {
        if (i == chr || i == start || i == end) continue;
        value=$i; gsub(/[|\r\n\t]/, "_", value);
        id=id "|" header[i] "=" value;
      }
      score=(fold && $fold != "") ? $fold : 0;
      print $chr, bed_start, $end, id, score, "+";
    }' "$source" > "$prepared"
}

while IFS= read -r -d '' source; do
  sample="$(basename "$(dirname "$source")")"
  prepared="$tmp_root/$(safe_label "macs2_${sample}").bed"
  if [[ "$source" == *.narrowPeak ]]; then format="narrow"; else format="broad"; fi
  prepare_bed_peak "$source" "$prepared" "$format"
  destination="$(dirname "$source")/${sample}_peaks_annotated.txt"
  run_annotation "MACS2" "$sample" "$source" "$prepared" "$destination"
done < <(find "$data_dir/macs2" -type f \( -name '*_peaks.narrowPeak' -o -name '*_peaks.broadPeak' \) -size +0c -print0 2>/dev/null)

while IFS= read -r -d '' source; do
  label="$(basename "$(dirname "$source")")"
  prepared="$tmp_root/$(safe_label "diffbind_${label}").bed"
  prepare_bed_peak "$source" "$prepared" "diffbind"
  destination="${source%.with_stats.bed}_annotated_with_stats.txt"
  run_annotation "Differential peaks" "$label" "$source" "$prepared" "$destination"
done < <(find "$data_dir/diffbind" -type f -name '*.with_stats.bed' -size +0c -print0 2>/dev/null)

while IFS= read -r -d '' source; do
  prepared_sibling="${source%.txt}.with_stats.bed"
  if [[ -s "$prepared_sibling" ]]; then
    continue
  fi
  label="$(basename "$(dirname "$source")")"
  prepared="$tmp_root/$(safe_label "diffbind_legacy_${label}").bed"
  prepare_diffbind_tsv "$source" "$prepared"
  destination="${source%.txt}_annotated_with_stats.txt"
  run_annotation "Differential peaks" "$label" "$source" "$prepared" "$destination"
done < <(find "$data_dir/diffbind" -type f -name 'DifferentialPeaks_*_ref.txt' -size +0c -print0 2>/dev/null)

while IFS= read -r -d '' source; do
  label="$(basename "$(dirname "$source")")"
  prepared="$tmp_root/$(safe_label "cutrun_all_${label}").bed"
  prepare_diffbind_tsv "$source" "$prepared"
  destination="${source%.tsv}_annotated_with_stats.txt"
  run_annotation "Differential peaks (all)" "$label" "$source" "$prepared" "$destination"
done < <(find "$data_dir/cutrun_diffbind" -type f -name 'all_differential_peaks.tsv' -size +0c -print0 2>/dev/null)

while IFS= read -r -d '' source; do
  label="$(basename "$(dirname "$source")")"
  prepared="$tmp_root/$(safe_label "cutrun_significant_${label}").bed"
  prepare_bed_peak "$source" "$prepared" "diffbind"
  destination="${source%.bed}_annotated_with_stats.txt"
  run_annotation "Differential peaks (significant)" "$label" "$source" "$prepared" "$destination"
done < <(find "$data_dir/cutrun_diffbind" -type f -name 'significant_differential_peaks.bed' -size +0c -print0 2>/dev/null)

if ((annotated_count == 0)); then
  echo "ERROR: no non-empty MACS2 or differential peak files were found under: $data_dir" >&2
  exit 1
fi

mv -f "$summary_tmp" "$summary"
{
  printf 'status\tcomplete\n'
  printf 'annotated_files\t%s\n' "$annotated_count"
  printf 'total_peaks\t%s\n' "$total_peaks"
  printf 'genome\t%s\n' "$genome"
  printf 'gtf\t%s\n' "$gtf"
  printf 'completed_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
} > "$complete_marker"
rm -f "$running_marker"

echo "Annotated $annotated_count peak files containing $total_peaks total intervals."
echo "Summary: $summary"
