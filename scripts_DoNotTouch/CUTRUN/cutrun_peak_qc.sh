#!/usr/bin/env bash
set -euo pipefail

seacr_dir="$1"
bowtie2_dir="$2"
out_dir="$3"
project_name="$4"

module load EBModules
module load BEDTools/2.30.0-GCC-10.3.0

mkdir -p "$out_dir"
peak_list="${out_dir}/seacr_peak_files.txt"
bam_list="${out_dir}/bowtie2_bam_files.txt"
summary_out="${out_dir}/cutrun_peak_qc_summary.txt"

find "$seacr_dir" -mindepth 2 -type f \( -name "*.stringent.bed" -o -name "*.relaxed.bed" \) | sort > "$peak_list"
find "$bowtie2_dir" -mindepth 2 -type f -name "*Aligned.sortedByCoord.out.bam" | sort > "$bam_list"

peak_count="$(wc -l < "$peak_list")"
bam_count="$(wc -l < "$bam_list")"
if [[ "$peak_count" -eq 0 ]]; then
  echo "ERROR: No SEACR peak BED files found under $seacr_dir" >&2
  exit 2
fi

while IFS= read -r peak_file; do
  sample="$(basename "$(dirname "$peak_file")")"
  awk -v s="$sample" 'BEGIN{OFS="\t"} NF>=3 {print $1,$2,$3,s}' "$peak_file"
done < "$peak_list" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -c 4 -o distinct,count \
  > "${out_dir}/seacr_consensus_peaks.bed"

if [[ -s "${out_dir}/seacr_consensus_peaks.bed" && "$bam_count" -gt 0 ]]; then
  mapfile -t bams < "$bam_list"
  {
    printf "chrom\tstart\tend\tsamples\tpeak_support"
    for bam in "${bams[@]}"; do
      printf "\t%s" "$(basename "$(dirname "$bam")")"
    done
    printf "\n"
    bedtools multicov -bams "${bams[@]}" -bed "${out_dir}/seacr_consensus_peaks.bed"
  } > "${out_dir}/seacr_consensus_peak_counts.tsv"
fi

{
  echo -e "project\t${project_name}"
  echo -e "seacr_dir\t${seacr_dir}"
  echo -e "bowtie2_dir\t${bowtie2_dir}"
  echo -e "peak_files\t${peak_count}"
  echo -e "bam_files\t${bam_count}"
  echo -e "consensus_peaks\t$(wc -l < "${out_dir}/seacr_consensus_peaks.bed")"
  echo -e "consensus_bed\t${out_dir}/seacr_consensus_peaks.bed"
  echo -e "consensus_counts\t${out_dir}/seacr_consensus_peak_counts.tsv"
} > "$summary_out"

frip_out="${out_dir}/seacr_frip_summary.tsv"
echo -e "sample\tfrip\tfragments_in_peaks\ttotal_fragments\tpeak_file" > "$frip_out"
find "$seacr_dir" -mindepth 2 -type f -name "*_seacr_summary.txt" | sort | while IFS= read -r summary; do
  sample="$(basename "$(dirname "$summary")")"
  frip="$(awk -F'\t' '$1=="frip"{print $2}' "$summary")"
  in_peaks="$(awk -F'\t' '$1=="fragments_in_peaks"{print $2}' "$summary")"
  total="$(awk -F'\t' '$1=="total_fragments"{print $2}' "$summary")"
  peak="$(awk -F'\t' '$1=="target_bedgraph"{print $2}' "$summary")"
  echo -e "${sample}\t${frip:-NA}\t${in_peaks:-NA}\t${total:-NA}\t${peak:-NA}"
done >> "$frip_out"
