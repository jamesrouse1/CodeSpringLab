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
normalization_mode="${10:-CPM}"
spikein_index="${11:-none}"
spikein_name="${12:-spikein}"
spikein_min_reads="${13:-1000}"

# Match the tested Radutama CodeSpringLab Bowtie2/Picard module stack.
module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

PICARD_JAR="${EBROOTPICARD:-}/picard.jar"
if [[ ! -s "$PICARD_JAR" ]]; then
  echo "ERROR: Picard 2.21.6 loaded, but picard.jar was not found at: $PICARD_JAR" >&2
  exit 127
fi

out_dir="$(dirname "$out_prefix")"
sample="$(basename "$out_prefix")"
mkdir -p "$out_dir"
tmp_dir="${out_dir}/tmp_cutrun_${SLURM_JOB_ID:-$$}"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"
export TEMP="$tmp_dir"
export TMP="$tmp_dir"
cleanup() {
  rm -rf -- "$tmp_dir"
}
trap cleanup EXIT
normalization_mode="$(echo "$normalization_mode" | tr '[:upper:]' '[:lower:]')"
if [[ "$normalization_mode" != "cpm" && "$normalization_mode" != "spikein" && "$normalization_mode" != "none" ]]; then
  echo "ERROR: normalization_mode must be one of CPM, spikein, or none. Got: $normalization_mode" >&2
  exit 2
fi
if [[ "$normalization_mode" == "spikein" && ( "$spikein_index" == "none" || ! -s "${spikein_index}.1.bt2" ) ]]; then
  echo "ERROR: spike-in normalization was requested but no readable Bowtie2 spike-in index was provided: $spikein_index" >&2
  exit 2
fi

bowtie2 --very-sensitive --threads 8 \
  --no-unal --end-to-end --phred33 \
  -x "$genome_index" \
  -U "$read1" 2> "${out_prefix}Log.final.out" \
  | samtools view -h -q "$mapq" -bS - \
  | samtools sort -T "${tmp_dir}/target_coordinate" -o "${out_prefix}Aligned.sortedByCoord.out.bam" -

samtools index -b "${out_prefix}Aligned.sortedByCoord.out.bam"
samtools idxstats "${out_prefix}Aligned.sortedByCoord.out.bam" > "${out_prefix}_chr_counts.txt"

java -Djava.io.tmpdir="$tmp_dir" -jar "$PICARD_JAR" MarkDuplicates \
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
  | sort -T "$tmp_dir" -S 4G -k1,1 -k2,2n > "${out_prefix}_fragments.bed"

bedtools genomecov -bg -i "${out_prefix}_fragments.bed" -g "$chromsize" > "${out_prefix}_fragments.raw.bedgraph"

# Always write CPM from this exact fragment coverage, including for spike-in
# runs, so SEACR and visualization can use one common pre-normalized track.
fragment_count_for_scale="$(wc -l < "${out_prefix}_fragments.bed")"
cpm_scale_factor="$(awk -v c="$fragment_count_for_scale" 'BEGIN{if (c>0) printf "%.10f", 1000000/c; else print "0"}')"
awk -v sf="$cpm_scale_factor" 'BEGIN{OFS="\t"} {$4=$4*sf; print}' "${out_prefix}_fragments.raw.bedgraph" > "${out_prefix}_fragments.CPM.bedgraph"

spikein_mapped_reads="0"
spikein_scale_factor="NA"
normalized_bedgraph="${out_prefix}_fragments.CPM.bedgraph"
if [[ "$normalization_mode" == "spikein" ]]; then
  bowtie2 --very-sensitive --threads 8 \
    --no-unal --end-to-end --phred33 \
    -x "$spikein_index" \
    -U "$read1" 2> "${out_prefix}_${spikein_name}Log.final.out" \
    | samtools view -h -bS - \
    | samtools sort -T "${tmp_dir}/spikein_coordinate" -o "${out_prefix}_${spikein_name}.bam" -
  samtools index -b "${out_prefix}_${spikein_name}.bam"
  # Use the same MAPQ threshold as the primary-genome signal when deriving
  # the spike-in scale factor; low-MAPQ E. coli alignments are not reliable
  # quantitative spike-in observations.
  spikein_mapped_reads="$(samtools view -c -F 4 -q "$mapq" "${out_prefix}_${spikein_name}.bam")"
  if [[ "$spikein_mapped_reads" -lt "$spikein_min_reads" ]]; then
    echo "WARNING: spike-in mapped reads (${spikein_mapped_reads}) are below the requested minimum (${spikein_min_reads}). Spike-in-normalized bedGraph will still be written, but interpretation should be cautious." >&2
  fi
  spikein_scale_factor="$(awk -v c="$spikein_mapped_reads" 'BEGIN{if (c>0) printf "%.10f", 10000/c; else print "0"}')"
  awk -v sf="$spikein_scale_factor" 'BEGIN{OFS="\t"} {$4=$4*sf; print}' "${out_prefix}_fragments.raw.bedgraph" > "${out_prefix}_fragments.spikein.bedgraph"
  normalized_bedgraph="${out_prefix}_fragments.spikein.bedgraph"
elif [[ "$normalization_mode" == "none" ]]; then
  normalized_bedgraph="${out_prefix}_fragments.raw.bedgraph"
fi

module load deepTools/3.5.2-foss-2022a
 bamCoverage -b "$bam_for_signal" \
  --normalizeUsing CPM \
  --binSize 10 \
  --ignoreForNormalization chrM MT \
  --outFileFormat bigwig \
  --outFileName "${out_prefix}_fragments.CPM.bw"
if [[ "$normalization_mode" == "none" ]]; then
  bamCoverage -b "$bam_for_signal" \
    --normalizeUsing None \
    --binSize 10 \
    --ignoreForNormalization chrM MT \
    --outFileFormat bigwig \
    --outFileName "${out_prefix}_fragments.raw.bw"
elif [[ "$normalization_mode" == "spikein" ]]; then
  bamCoverage -b "$bam_for_signal" \
    --normalizeUsing None \
    --scaleFactor "$spikein_scale_factor" \
    --binSize 10 \
    --ignoreForNormalization chrM MT \
    --outFileFormat bigwig \
    --outFileName "${out_prefix}_fragments.spikein.bw"
fi

mapped_reads="$(samtools view -c -F 4 "${out_prefix}Aligned.sortedByCoord.out.bam")"
dedup_reads="$(samtools view -c -F 4 "${out_prefix}Aligned.sortedByCoord_removeDup.out.bam")"
fragment_count="$(wc -l < "${out_prefix}_fragments.bed")"
duplicate_fraction="$(awk 'BEGIN{v="NA"} /^PERCENT_DUPLICATION/ {getline; split($0,a,"\t"); v=a[9]} END{print v}' "${out_prefix}_markedDup_metrics.txt")"
awk '{print $3-$2}' "${out_prefix}_fragments.bed" > "${out_prefix}_fragment_lengths.txt"
awk 'BEGIN{OFS="\t"; print "metric","value"} {n++; s+=$1; if (n==1 || $1<min) min=$1; if ($1>max) max=$1} END{if(n>0){print "fragments",n; print "mean_fragment_length",s/n; print "min_fragment_length",min; print "max_fragment_length",max}}' "${out_prefix}_fragment_lengths.txt" > "${out_prefix}_fragment_length_summary.txt"
{
  echo -e "sample\t${sample}"
  echo -e "read1\t${read1}"
  echo -e "target_genome_index\t${genome_index}"
  echo -e "chrom_sizes\t${chromsize}"
  echo -e "mapq\t${mapq}"
  echo -e "dedup_mode\t${dedup_mode}"
  echo -e "normalization_mode\t${normalization_mode}"
  echo -e "normalized_bedgraph\t${normalized_bedgraph}"
  echo -e "raw_bedgraph\t${out_prefix}_fragments.raw.bedgraph"
  echo -e "cpm_bedgraph\t${out_prefix}_fragments.CPM.bedgraph"
  echo -e "cpm_scale_factor\t${cpm_scale_factor}"
  echo -e "spikein_name\t${spikein_name}"
  echo -e "spikein_index\t${spikein_index}"
  echo -e "spikein_mapq\t${mapq}"
  echo -e "spikein_mapped_reads\t${spikein_mapped_reads}"
  echo -e "spikein_scale_factor\t${spikein_scale_factor}"
  echo -e "duplicate_fraction\t${duplicate_fraction}"
  echo -e "mapped_reads\t${mapped_reads}"
  echo -e "deduplicated_reads\t${dedup_reads}"
  echo -e "fragments_used_for_signal\t${fragment_count}"
} > "${out_prefix}_alignment_summary.txt"
