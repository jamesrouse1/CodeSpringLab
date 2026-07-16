set -Eeuo pipefail

prefix="${1:?ERROR: output prefix is required}"
assay_label="${PIPELINE_ASSAY_LABEL:-ATAC}"
coverage_ignore_chroms="${BAMCOVERAGE_IGNORE_CHROMS:-chrM chrX}"
bowtie_log="${prefix}Log.final.out"
current_stage="initialization"
sample_dir="$(dirname "${prefix}")"
bamcoverage_tmp_dir="${sample_dir}/tmp_bamcoverage_${SLURM_JOB_ID:-$$}"
bigwig_tmp="${prefix}Aligned.sortedByCoord_removeDup.out.bw.${SLURM_JOB_ID:-$$}.tmp"

cleanup() {
	rm -f -- "${bigwig_tmp}"
	rm -rf -- "${bamcoverage_tmp_dir}"
}

report_stage() {
	current_stage="$1"
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${assay_label} Bowtie2: ${current_stage} ($(basename "${prefix}"))"
}

report_failure() {
	local rc=$?
	echo "ERROR: ${assay_label} Bowtie2 failed during ${current_stage} for $(basename "${prefix}") (exit ${rc})." >&2
	if [[ -s "$bowtie_log" ]]; then
		echo "Last 80 lines of ${bowtie_log}:" >&2
		tail -n 80 "$bowtie_log" >&2 || true
	fi
	return "$rc"
}

trap report_failure ERR
trap cleanup EXIT

report_stage "loading modules"
if ! type module >/dev/null 2>&1; then
	for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /cm/local/apps/environment-modules/current/init/bash; do
		[[ -s "$module_init" ]] && source "$module_init" && break
	done
fi
type module >/dev/null 2>&1 || { echo "ERROR: cluster module command is unavailable." >&2; exit 127; }
module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

report_stage "aligning paired-end reads"
bowtie2  --very-sensitive -k 10 -X 1000 --dovetail --threads 8 \
		-x $2 \
		-1 $3 -2 $4 > ${1}_temp1.sam 2> ${1}Log.final.out

report_stage "filtering and sorting alignments"
samtools view -h -q 30 -f 0x2 -o ${1}_temp2.sam ${1}_temp1.sam
rm ${1}_temp1.sam

samtools sort -o ${1}_temp2.sortedByCoord.sam ${1}_temp2.sam
samtools idxstats ${1}_temp2.sortedByCoord.sam > ${1}_chr_counts.txt
rm ${1}_temp2.sortedByCoord.sam

samtools sort -n -o ${1}_temp3.sam ${1}_temp2.sam
rm ${1}_temp2.sam

awk '( ($3 ~ /^chr/ && $3 != "chrM" && $3 != "chrUn") || (/^@/) )' ${1}_temp3.sam > ${1}Aligned.sortedByName.out.sam
rm ${1}_temp3.sam

samtools view -h -bS -o ${1}Aligned.sortedByName.out.bam ${1}Aligned.sortedByName.out.sam
rm ${1}Aligned.sortedByName.out.sam

samtools sort -o ${1}Aligned.sortedByCoord.out.bam ${1}Aligned.sortedByName.out.bam
rm ${1}Aligned.sortedByName.out.bam

samtools index -b ${1}Aligned.sortedByCoord.out.bam ${1}Aligned.sortedByCoord.out.bam.bai
samtools quickcheck -v ${1}Aligned.sortedByCoord.out.bam

#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#REMOVE_DUPLICATES=true \
#I=${1}Aligned.sortedByName.out.bam \
#O=${1}Aligned.sortedByName_removeDup.out.bam \
#M=${1}_markedDup_metrics.txt

report_stage "removing duplicates with Picard"
java -Xmx80g -Djava.io.tmpdir="$(dirname "${1}")" -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
REMOVE_DUPLICATES=true \
I=${1}Aligned.sortedByCoord.out.bam \
O=${1}Aligned.sortedByCoord_removeDup.out.bam \
M=${1}_markedDup_metrics.txt

if [[ ! -s "${1}Aligned.sortedByCoord_removeDup.out.bam" ]]; then
  echo "ERROR: Picard produced an empty deduplicated BAM for $(basename "${1}")" >&2
  exit 1
fi
samtools quickcheck -v ${1}Aligned.sortedByCoord_removeDup.out.bam

report_stage "collecting insert-size metrics"
java -Xmx80g -Djava.io.tmpdir="$(dirname "${1}")" -jar "$EBROOTPICARD/picard.jar" CollectInsertSizeMetrics \
I=${1}Aligned.sortedByCoord_removeDup.out.bam \
O=${1}_insert_size_metrics.txt \
H=${1}_insert_size_histogram.pdf \
M=0.5

pdftoppm -jpeg ${1}_insert_size_histogram.pdf ${1}_insert_size_histogram
mv ${1}_insert_size_histogram-1.jpg ${1}_insert_size_histogram.jpg
rm ${1}_insert_size_histogram.pdf

samtools index -b ${1}Aligned.sortedByCoord_removeDup.out.bam ${1}Aligned.sortedByCoord_removeDup.out.bam.bai

#bedtools bamtobed -i ${1}Aligned.sortedByName.out.bam > ${1}Aligned.sortedByName.out.bed
#bedtools bamtobed -i ${1}Aligned.sortedByName_removeDup.out.bam > ${1}Aligned.sortedByName_removeDup.out.bed
report_stage "creating BED files"
bedtools bamtobed -i ${1}Aligned.sortedByCoord.out.bam > ${1}Aligned.sortedByCoord.out.bed
bedtools bamtobed -i ${1}Aligned.sortedByCoord_removeDup.out.bam > ${1}Aligned.sortedByCoord_removeDup.out.bed

#####################################

report_stage "loading deepTools"
module load EBModules
module load deepTools/3.5.2-foss-2022a

#--effectiveGenomeSize 2913022398 \ # Mice:2150570000; GRCh38:2913022398
### For male mice chrX should be ignored

report_stage "creating CPM bigWig"
mkdir -p "${bamcoverage_tmp_dir}"
export TMPDIR="${bamcoverage_tmp_dir}"
export TEMP="${bamcoverage_tmp_dir}"
export TMP="${bamcoverage_tmp_dir}"
rm -f "${bigwig_tmp}"
coverage_args=(
  -b "${1}Aligned.sortedByCoord_removeDup.out.bam"
  -o "${bigwig_tmp}"
  --outFileFormat bigwig
  --normalizeUsing CPM
  --numberOfProcessors 2
  --binSize 10
  --extendReads
)
if [[ -n "${coverage_ignore_chroms// }" ]]; then
  read -r -a ignore_chroms <<< "$coverage_ignore_chroms"
  coverage_args+=(--ignoreForNormalization "${ignore_chroms[@]}")
fi
if ! bamCoverage "${coverage_args[@]}"; then
  rm -f "${bigwig_tmp}"
  echo "ERROR: direct CPM bigWig generation failed for $(basename "${1}")" >&2
  exit 1
fi

if [[ ! -s "${bigwig_tmp}" ]]; then
  echo "ERROR: bigWig conversion produced an empty file for $(basename "${1}")" >&2
  exit 1
fi
mv "${bigwig_tmp}" "${1}Aligned.sortedByCoord_removeDup.out.bw"

report_stage "writing alignment summary"
mapped_reads="$(samtools view -c -F 4 ${1}Aligned.sortedByCoord.out.bam)"
deduplicated_reads="$(samtools view -c -F 4 ${1}Aligned.sortedByCoord_removeDup.out.bam)"
{
  echo -e "sample\t$(basename "${1}")"
  echo -e "reference_index\t${2}"
  echo -e "mapped_reads\t${mapped_reads}"
  echo -e "deduplicated_reads\t${deduplicated_reads}"
  echo -e "effective_genome_size\t${5}"
  echo -e "bigwig_normalization\tCPM"
  echo -e "bigwig\t${1}Aligned.sortedByCoord_removeDup.out.bw"
} > ${1}_alignment_summary.txt

report_stage "complete"
trap - ERR
