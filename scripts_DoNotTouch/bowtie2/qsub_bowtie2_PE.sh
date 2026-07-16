#!/usr/bin/env bash
#SBATCH --job-name=bowtie2
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

runner="${8:-}"
if [[ -z "$runner" || ! -s "$runner" ]]; then
  runner="${SLURM_SUBMIT_DIR:-$PWD}/../scripts_DoNotTouch/bowtie2/bowtie2_PE.sh"
fi
[[ -s "$runner" ]] || { echo "ERROR: Bowtie2 runner not found: $runner" >&2; exit 2; }
exec bash "$runner" "$1" "$2" "$3" "$4" "$5" "$6" "$7"

# Ori script with grid qsub
#qsub -l mem_free=50G -pe threads 8 -cwd -o ../../csl_results/${7}/log/output_bowtie2.txt -e ${1}Log.final.out -V ../scripts_DoNotTouch/bowtie2/bowtie2_PE.sh $1 $2 $3 $4 $5 $6 $7
