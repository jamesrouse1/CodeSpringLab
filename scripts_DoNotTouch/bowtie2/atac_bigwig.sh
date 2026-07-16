#!/bin/bash

set -euo pipefail

prefix="$1"
reference_index="$2"
sample="$(basename "$prefix")"
sample_dir="$(dirname "$prefix")"
input_bam="${prefix}Aligned.sortedByCoord.out.bam"
dedup_bam="${prefix}Aligned.sortedByCoord_removeDup.out.bam"
dedup_bai="${dedup_bam}.bai"
bigwig="${prefix}Aligned.sortedByCoord_removeDup.out.bw"
alignment_summary="${prefix}_alignment_summary.txt"
repair_summary="${prefix}_postprocess_summary.txt"
tmp_bigwig="${prefix}.bigwig.${SLURM_JOB_ID:-$$}.bw"
tmp_bai="${prefix}.bigwig.${SLURM_JOB_ID:-$$}.bam.bai"
tmp_summary="${prefix}.bigwig.${SLURM_JOB_ID:-$$}.summary.txt"
tmp_dir="${sample_dir}/tmp_bamcoverage_${SLURM_JOB_ID:-$$}"

cleanup() {
  rm -f "$tmp_bigwig" "$tmp_bai" "$tmp_summary"
  rm -rf -- "$tmp_dir"
}
trap cleanup EXIT

[[ -s "$dedup_bam" ]] || { echo "ERROR: deduplicated BAM is missing or empty: $dedup_bam" >&2; exit 1; }

mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"
export TEMP="$tmp_dir"
export TMP="$tmp_dir"

module load EBModules
module load deepTools/3.5.2-foss-2022a
module load SAMtools/1.14-GCC-10.3.0

samtools quickcheck -v "$dedup_bam"
if [[ ! -s "$dedup_bai" ]] || [[ "$dedup_bai" -ot "$dedup_bam" ]]; then
  samtools index "$dedup_bam" "$tmp_bai"
  [[ -s "$tmp_bai" ]] || { echo "ERROR: BAM index is empty for $sample" >&2; exit 1; }
  mv "$tmp_bai" "$dedup_bai"
fi

bamCoverage \
  -b "$dedup_bam" \
  -o "$tmp_bigwig" \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --numberOfProcessors 2 \
  --binSize 10 \
  --extendReads \
  --ignoreForNormalization chrM chrX

[[ -s "$tmp_bigwig" ]] || { echo "ERROR: direct CPM bigWig is empty for $sample" >&2; exit 1; }
mv "$tmp_bigwig" "$bigwig"

mapped_reads=""
if [[ -s "$input_bam" ]]; then
  samtools quickcheck -v "$input_bam"
  mapped_reads="$(samtools view -c -F 4 "$input_bam")"
fi
deduplicated_reads="$(samtools view -c -F 4 "$dedup_bam")"

{
  printf 'sample\t%s\n' "$sample"
  printf 'reference_index\t%s\n' "$reference_index"
  printf 'mapped_reads\t%s\n' "$mapped_reads"
  printf 'deduplicated_reads\t%s\n' "$deduplicated_reads"
  printf 'bigwig_normalization\tCPM\n'
  printf 'bigwig\t%s\n' "$bigwig"
} > "$tmp_summary"
mv "$tmp_summary" "$alignment_summary"

{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tbigwig_repaired\n'
  printf 'deduplicated_bam\t%s\n' "$dedup_bam"
  printf 'bigwig_normalization\tCPM\n'
  printf 'bigwig\t%s\n' "$bigwig"
} > "$repair_summary"

echo "CPM bigWig repair complete: $sample"
