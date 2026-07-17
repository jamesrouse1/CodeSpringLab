#!/usr/bin/env bash
set -Eeuo pipefail

prefix="${1:?ERROR: output prefix is required}"
index="${2:?ERROR: Bowtie2 index is required}"
read1="${3:?ERROR: read FASTQ is required}"
effective_genome_size="${5:?ERROR: effective genome size is required}"
sample_dir="$(dirname "$prefix")"
sample="$(basename "$prefix")"
job_key="${SLURM_JOB_ID:-$$}"
tmp_dir="${sample_dir}/.chip_bowtie2_tmp_${job_key}"
bowtie_log="${prefix}Log.final.out"
raw_bam="${prefix}Aligned.sortedByCoord.out.bam"
dedup_bam="${prefix}Aligned.sortedByCoord_removeDup.out.bam"
bigwig="${prefix}Aligned.sortedByCoord_removeDup.out.bw"
tmp_bigwig="${bigwig}.${job_key}.tmp"
current_stage="initialization"

cleanup() {
  rm -rf -- "$tmp_dir"
  rm -f -- "$tmp_bigwig"
}
failure() {
  rc=$?
  echo "ERROR: ChIP-seq Bowtie2 failed during ${current_stage} for ${sample} (exit ${rc})." >&2
  [[ -s "$bowtie_log" ]] && tail -n 80 "$bowtie_log" >&2 || true
  return "$rc"
}
trap failure ERR
trap cleanup EXIT
mkdir -p "$sample_dir" "$tmp_dir"
export TMPDIR="$tmp_dir" TMP="$tmp_dir" TEMP="$tmp_dir"

if ! type module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
    [[ -s "$module_init" ]] && source "$module_init" && break
  done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }

current_stage="loading alignment modules"
module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

current_stage="aligning single-end reads"
bowtie2 --very-sensitive --threads 8 --no-unal --end-to-end --phred33 -x "$index" -U "$read1" 2> "$bowtie_log" |
  samtools view -h -q 30 -F 4 - |
  awk 'BEGIN{OFS="\t"} /^@/ || ($3 ~ /^chr([0-9]+|X|Y)$/) {print}' |
  samtools view -b -o "${tmp_dir}/filtered.bam" -
samtools sort -@ 8 -T "${tmp_dir}/sort" -o "$raw_bam" "${tmp_dir}/filtered.bam"
samtools quickcheck -v "$raw_bam"
samtools index -b "$raw_bam" "${raw_bam}.bai"
samtools idxstats "$raw_bam" > "${prefix}_chr_counts.txt"

current_stage="removing duplicates"
java -Xmx80g -Djava.io.tmpdir="$tmp_dir" -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  REMOVE_DUPLICATES=true I="$raw_bam" O="$dedup_bam" M="${prefix}_markedDup_metrics.txt"
[[ -s "$dedup_bam" ]] || { echo "ERROR: Picard produced an empty deduplicated BAM for ${sample}." >&2; exit 1; }
samtools quickcheck -v "$dedup_bam"
samtools index -b "$dedup_bam" "${dedup_bam}.bai"

current_stage="creating BED files"
bedtools bamtobed -i "$raw_bam" > "${prefix}Aligned.sortedByCoord.out.bed"
bedtools bamtobed -i "$dedup_bam" > "${prefix}Aligned.sortedByCoord_removeDup.out.bed"

current_stage="creating CPM bigWig"
module load deepTools/3.5.2-foss-2022a
bamCoverage -b "$dedup_bam" -o "$tmp_bigwig" --outFileFormat bigwig --normalizeUsing CPM \
  --numberOfProcessors 2 --binSize 10 --extendReads 200 --ignoreForNormalization chrM
[[ -s "$tmp_bigwig" ]] || { echo "ERROR: bamCoverage produced an empty bigWig for ${sample}." >&2; exit 1; }
mv "$tmp_bigwig" "$bigwig"

current_stage="writing alignment summary"
mapped_reads="$(samtools view -c -F 4 "$raw_bam")"
deduplicated_reads="$(samtools view -c -F 4 "$dedup_bam")"
{
  printf 'sample\t%s\n' "$sample"
  printf 'reference_index\t%s\n' "$index"
  printf 'mapped_reads\t%s\n' "$mapped_reads"
  printf 'deduplicated_reads\t%s\n' "$deduplicated_reads"
  printf 'effective_genome_size\t%s\n' "$effective_genome_size"
  printf 'bigwig_normalization\tCPM\n'
  printf 'bigwig\t%s\n' "$bigwig"
} > "${prefix}_alignment_summary.txt"

current_stage="complete"
trap - ERR
