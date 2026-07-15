#!/bin/bash

set -euo pipefail

prefix="$1"
effective_genome_size="$2"
chrom_sizes="$3"
reference_index="$4"
sample="$(basename "$prefix")"
sample_dir="$(dirname "$prefix")"

input_bam="${prefix}Aligned.sortedByCoord.out.bam"
dedup_bam="${prefix}Aligned.sortedByCoord_removeDup.out.bam"
dedup_bai="${dedup_bam}.bai"
input_bed="${prefix}Aligned.sortedByCoord.out.bed"
dedup_bed="${prefix}Aligned.sortedByCoord_removeDup.out.bed"
bigwig="${prefix}Aligned.sortedByCoord_removeDup.out.bw"
dup_metrics="${prefix}_markedDup_metrics.txt"
insert_metrics="${prefix}_insert_size_metrics.txt"
insert_jpg="${prefix}_insert_size_histogram.jpg"
alignment_summary="${prefix}_alignment_summary.txt"
repair_summary="${prefix}_postprocess_summary.txt"

tmp_root="${prefix}.postprocess.${SLURM_JOB_ID:-$$}"
tmp_bam="${tmp_root}.dedup.bam"
tmp_bai="${tmp_root}.dedup.bam.bai"
tmp_input_bed="${tmp_root}.aligned.bed"
tmp_dedup_bed="${tmp_root}.dedup.bed"
tmp_dup_metrics="${tmp_root}.duplicate_metrics.txt"
tmp_insert_metrics="${tmp_root}.insert_metrics.txt"
tmp_insert_pdf="${tmp_root}.insert_histogram.pdf"
tmp_insert_jpg_prefix="${tmp_root}.insert_histogram"
tmp_bedgraph="${tmp_root}.bedgraph"
tmp_filtered_bedgraph="${tmp_root}.filtered.bedgraph"
tmp_bigwig="${tmp_root}.bw"

cleanup() {
  rm -f "$tmp_bam" "$tmp_bai" "$tmp_input_bed" "$tmp_dedup_bed" \
    "$tmp_dup_metrics" "$tmp_insert_metrics" "$tmp_insert_pdf" \
    "${tmp_insert_jpg_prefix}-1.jpg" "$tmp_bedgraph" \
    "$tmp_filtered_bedgraph" "$tmp_bigwig"
}
trap cleanup EXIT

mkdir -p "$sample_dir"
cd "$sample_dir"

if [[ ! -s "$input_bam" ]]; then
  echo "ERROR: aligned BAM is missing or empty: $input_bam" >&2
  exit 1
fi
if [[ ! -s "$chrom_sizes" ]]; then
  echo "ERROR: chromosome sizes file is missing or empty: $chrom_sizes" >&2
  exit 1
fi

module load EBModules
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

samtools quickcheck -v "$input_bam"

java -Xmx80g -Djava.io.tmpdir="$sample_dir" -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  REMOVE_DUPLICATES=true \
  I="$input_bam" \
  O="$tmp_bam" \
  M="$tmp_dup_metrics"

[[ -s "$tmp_bam" ]] || { echo "ERROR: Picard produced an empty deduplicated BAM for $sample" >&2; exit 1; }
samtools quickcheck -v "$tmp_bam"
samtools index -b "$tmp_bam" "$tmp_bai"

java -Xmx80g -Djava.io.tmpdir="$sample_dir" -jar "$EBROOTPICARD/picard.jar" CollectInsertSizeMetrics \
  I="$tmp_bam" \
  O="$tmp_insert_metrics" \
  H="$tmp_insert_pdf" \
  M=0.5

pdftoppm -jpeg "$tmp_insert_pdf" "$tmp_insert_jpg_prefix"
[[ -s "${tmp_insert_jpg_prefix}-1.jpg" ]] || { echo "ERROR: insert-size image was not created for $sample" >&2; exit 1; }

bedtools bamtobed -i "$input_bam" > "$tmp_input_bed"
bedtools bamtobed -i "$tmp_bam" > "$tmp_dedup_bed"
[[ -s "$tmp_dedup_bed" ]] || { echo "ERROR: deduplicated BED is empty for $sample" >&2; exit 1; }

module load deepTools/3.5.2-foss-2022a

bamCoverage -b "$tmp_bam" \
  --normalizeUsing RPGC \
  --effectiveGenomeSize "$effective_genome_size" \
  --binSize 10 \
  --extendReads \
  --ignoreForNormalization chrX chrM \
  --outFileFormat bedgraph \
  --outFileName "$tmp_bedgraph"

awk 'BEGIN{OFS="\t"}
     NR==FNR {valid[$1]=1; next}
     ($1 in valid) && $1 != "chrM" && $1 != "chrUn" &&
       $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ &&
       $4 ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ {
         print $1, $2, $3, $4
       }' "$chrom_sizes" "$tmp_bedgraph" > "$tmp_filtered_bedgraph"

[[ -s "$tmp_filtered_bedgraph" ]] || { echo "ERROR: filtered bedGraph is empty for $sample" >&2; exit 1; }

/grid/bsr/data/data/utama/tools/bin/x86_64/bedGraphToBigWig \
  "$tmp_filtered_bedgraph" "$chrom_sizes" "$tmp_bigwig"
[[ -s "$tmp_bigwig" ]] || { echo "ERROR: bigWig is empty for $sample" >&2; exit 1; }

mapped_reads="$(samtools view -c -F 4 "$input_bam")"
deduplicated_reads="$(samtools view -c -F 4 "$tmp_bam")"

# Publish the complete repaired output set only after every command succeeds.
mv "$tmp_bam" "$dedup_bam"
mv "$tmp_bai" "$dedup_bai"
mv "$tmp_input_bed" "$input_bed"
mv "$tmp_dedup_bed" "$dedup_bed"
mv "$tmp_dup_metrics" "$dup_metrics"
mv "$tmp_insert_metrics" "$insert_metrics"
mv "${tmp_insert_jpg_prefix}-1.jpg" "$insert_jpg"
mv "$tmp_bigwig" "$bigwig"

{
  printf 'sample\t%s\n' "$sample"
  printf 'reference_index\t%s\n' "$reference_index"
  printf 'mapped_reads\t%s\n' "$mapped_reads"
  printf 'deduplicated_reads\t%s\n' "$deduplicated_reads"
  printf 'effective_genome_size\t%s\n' "$effective_genome_size"
  printf 'bigwig\t%s\n' "$bigwig"
} > "$alignment_summary"

{
  printf 'sample\t%s\n' "$sample"
  printf 'status\tcomplete\n'
  printf 'input_bam\t%s\n' "$input_bam"
  printf 'deduplicated_bam\t%s\n' "$dedup_bam"
  printf 'deduplicated_bed\t%s\n' "$dedup_bed"
  printf 'bigwig\t%s\n' "$bigwig"
} > "$repair_summary"

echo "Post-alignment repair complete: $sample"
