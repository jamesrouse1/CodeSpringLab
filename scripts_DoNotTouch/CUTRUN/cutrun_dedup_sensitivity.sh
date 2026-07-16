#!/usr/bin/env bash
set -euo pipefail

dedup_bam="$1"
chrom_sizes="$2"
alignment_summary="$3"
control_bedgraph="$4"
seacr_norm="${5:-non}"
stringency="${6:-stringent}"
out_prefix="$7"
project_name="$8"
max_fragment="${9:-1000}"
remove_mito="${10:-y}"
seacr_runner="${11}"

for required in "$dedup_bam" "$chrom_sizes" "$alignment_summary" "$seacr_runner"; do
  [[ -s "$required" ]] || { echo "ERROR: required sensitivity-analysis input is missing or empty: $required" >&2; exit 2; }
done
if [[ "$control_bedgraph" != "none" && ! -s "$control_bedgraph" ]]; then
  echo "ERROR: matched control bedGraph is missing or empty: $control_bedgraph" >&2
  exit 2
fi

module load EBModules
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

out_dir="$(dirname "$out_prefix")"
sample="$(basename "$out_prefix")"
mkdir -p "$out_dir"
tmp_dir="${out_dir}/tmp_dedup_sensitivity_${SLURM_JOB_ID:-$$}"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"
export TEMP="$tmp_dir"
export TMP="$tmp_dir"
cleanup() { rm -rf -- "$tmp_dir"; }
trap cleanup EXIT

metric() {
  awk -F '\t' -v key="$2" '$1==key {print $2; exit}' "$1"
}

normalization_mode="$(metric "$alignment_summary" normalization_mode | tr '[:upper:]' '[:lower:]')"
spikein_scale_factor="$(metric "$alignment_summary" spikein_scale_factor)"
if [[ "$normalization_mode" != "spikein" && "$normalization_mode" != "cpm" && "$normalization_mode" != "none" ]]; then
  echo "ERROR: alignment summary has unsupported normalization_mode: ${normalization_mode:-missing}" >&2
  exit 2
fi
if [[ "$normalization_mode" == "spikein" ]] && ! awk -v x="$spikein_scale_factor" 'BEGIN{exit !(x+0>0)}'; then
  echo "ERROR: alignment summary does not contain a positive spikein_scale_factor." >&2
  exit 2
fi

samtools quickcheck -v "$dedup_bam"
name_sorted_bam="${tmp_dir}/${sample}.name_sorted.bam"
tmp_fragments="${tmp_dir}/${sample}_fragments.bed"
samtools sort -n -@ 8 -T "${tmp_dir}/name_sort" -o "$name_sorted_bam" "$dedup_bam"
bedtools bamtobed -bedpe -i "$name_sorted_bam" \
  | awk -v maxfrag="$max_fragment" -v rmmito="$remove_mito" 'BEGIN{OFS="\t"} $1==$4 && $6>$2 {frag=$6-$2; if (frag<=maxfrag && frag>0) {if (rmmito=="y" && ($1=="chrM" || $1=="MT")) next; print $1,$2,$6}}' \
  | sort -T "$tmp_dir" -S 4G -k1,1 -k2,2n > "$tmp_fragments"
[[ -s "$tmp_fragments" ]] || { echo "ERROR: deduplicated fragment generation produced an empty BED for $sample" >&2; exit 2; }
mv "$tmp_fragments" "${out_prefix}_fragments.bed"

bedtools genomecov -bg -i "${out_prefix}_fragments.bed" -g "$chrom_sizes" > "${out_prefix}_fragments.raw.bedgraph"
normalized_bedgraph="${out_prefix}_fragments.raw.bedgraph"
scale_factor="1"
if [[ "$normalization_mode" == "spikein" ]]; then
  scale_factor="$spikein_scale_factor"
  awk -v sf="$scale_factor" 'BEGIN{OFS="\t"} {$4=$4*sf; print}' "${out_prefix}_fragments.raw.bedgraph" > "${out_prefix}_fragments.spikein.bedgraph"
  normalized_bedgraph="${out_prefix}_fragments.spikein.bedgraph"
elif [[ "$normalization_mode" == "cpm" ]]; then
  fragment_count="$(wc -l < "${out_prefix}_fragments.bed")"
  scale_factor="$(awk -v c="$fragment_count" 'BEGIN{if(c>0) printf "%.10f",1000000/c; else print "0"}')"
  awk -v sf="$scale_factor" 'BEGIN{OFS="\t"} {$4=$4*sf; print}' "${out_prefix}_fragments.raw.bedgraph" > "${out_prefix}_fragments.CPM.bedgraph"
  normalized_bedgraph="${out_prefix}_fragments.CPM.bedgraph"
fi

seacr_target_bedgraph="$normalized_bedgraph"
if [[ "$seacr_norm" == "norm" ]]; then
  seacr_target_bedgraph="${out_prefix}_fragments.raw.bedgraph"
fi

module load deepTools/3.5.2-foss-2022a
bigwig="${out_prefix}_fragments.${normalization_mode}.bw"
bamCoverage -b "$dedup_bam" --normalizeUsing None --scaleFactor "$scale_factor" \
  --binSize 10 --ignoreForNormalization chrM MT --outFileFormat bigwig --outFileName "$bigwig"

source "$seacr_runner" "$seacr_target_bedgraph" "$control_bedgraph" "$seacr_norm" "$stringency" "$out_prefix" "$project_name" "${out_prefix}_fragments.bed"

{
  echo -e "analysis\tdeduplicated_target_sensitivity"
  echo -e "source_bam\t${dedup_bam}"
  echo -e "source_alignment_summary\t${alignment_summary}"
  echo -e "normalization_mode\t${normalization_mode}"
  echo -e "scale_factor\t${scale_factor}"
  echo -e "normalized_bedgraph\t${normalized_bedgraph}"
  echo -e "seacr_target_bedgraph\t${seacr_target_bedgraph}"
  echo -e "bigwig\t${bigwig}"
} >> "${out_prefix}_seacr_summary.txt"
