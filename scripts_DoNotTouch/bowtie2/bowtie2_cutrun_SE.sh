#!/usr/bin/env bash
set -euo pipefail

out_prefix="$1"
genome_index="$2"
read1="$3"
chromsize="$4"
project_name="$5"
mapq="${6:-30}"
max_fragment="${7:-1000}"
dedup_mode="${8:-keepdup}"
remove_mito="${9:-y}"

module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

out_dir="$(dirname "$out_prefix")"
sample="$(basename "$out_prefix")"
mkdir -p "$out_dir"

bowtie2 --very-sensitive --threads 8 \
  --no-unal --end-to-end --phred33 \
  -x "$genome_index" \
  -U "$read1" 2> "${out_prefix}Log.final.out" \
  | samtools view -h -q "$mapq" -bS - \
  | samtools sort -o "${out_prefix}Aligned.sortedByCoord.out.bam" -

samtools index -b "${out_prefix}Aligned.sortedByCoord.out.bam"
samtools idxstats "${out_prefix}Aligned.sortedByCoord.out.bam" > "${out_prefix}_chr_counts.txt"

java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  REMOVE_DUPLICATES=true \
  I="${out_prefix}Aligned.sortedByCoord.out.bam" \
  O="${out_prefix}Aligned.sortedByCoord_removeDup.out.bam" \
  M="${out_prefix}_markedDup_metrics.txt" \
  VALIDATION_STRINGENCY=SILENT
samtools index -b "${out_prefix}Aligned.sortedByCoord_removeDup.out.bam"

bam_for_signal="${out_prefix}Aligned.sortedByCoord.out.bam"
if [[ "$dedup_mode" == "dedup" ]]; then
  bam_for_signal="${out_prefix}Aligned.sortedByCoord_removeDup.out.bam"
fi

bedtools bamtobed -i "${out_prefix}Aligned.sortedByCoord.out.bam" > "${out_prefix}Aligned.sortedByCoord.out.bed"
bedtools bamtobed -i "${out_prefix}Aligned.sortedByCoord_removeDup.out.bam" > "${out_prefix}Aligned.sortedByCoord_removeDup.out.bed"
bedtools bamtobed -i "$bam_for_signal" \
  | awk -v rmmito="$remove_mito" 'BEGIN{OFS="\t"} {if (rmmito=="y" && ($1=="chrM" || $1=="MT")) next; print $1,$2,$3}' \
  | sort -k1,1 -k2,2n > "${out_prefix}_fragments.bed"

bedtools genomecov -bg -i "${out_prefix}_fragments.bed" -g "$chromsize" > "${out_prefix}_fragments.raw.bedgraph"

module load deepTools/3.5.2-foss-2022a
bamCoverage -b "$bam_for_signal" \
  --normalizeUsing CPM \
  --binSize 10 \
  --ignoreForNormalization chrM MT \
  --outFileFormat bigwig \
  --outFileName "${out_prefix}_fragments.CPM.bw"

mapped_reads="$(samtools view -c -F 4 "${out_prefix}Aligned.sortedByCoord.out.bam")"
dedup_reads="$(samtools view -c -F 4 "${out_prefix}Aligned.sortedByCoord_removeDup.out.bam")"
fragment_count="$(wc -l < "${out_prefix}_fragments.bed")"
{
  echo -e "sample\t${sample}"
  echo -e "read1\t${read1}"
  echo -e "mapq\t${mapq}"
  echo -e "dedup_mode\t${dedup_mode}"
  echo -e "mapped_reads\t${mapped_reads}"
  echo -e "deduplicated_reads\t${dedup_reads}"
  echo -e "fragments_used_for_signal\t${fragment_count}"
} > "${out_prefix}_alignment_summary.txt"
