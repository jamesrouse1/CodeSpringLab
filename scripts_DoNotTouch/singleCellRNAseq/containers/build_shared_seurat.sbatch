#!/usr/bin/env bash
#SBATCH --job-name=build_csl_seurat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=codespring_seurat_build_%j.out
#SBATCH --error=codespring_seurat_build_%j.err

set -Eeuo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
definition="$script_dir/codespring-seurat_1.0.0.def"
target_dir="${CSL_SHARED_CONTAINER_DIR:-/grid/bsr/data/data/bsr_readable_data/containers/seurat}"
target="$target_dir/codespring-seurat_1.0.0.sif"
building="$target_dir/.codespring-seurat_1.0.0.${SLURM_JOB_ID:-manual}.building.sif"

if [[ -e "$target" ]]; then
  echo "ERROR: Refusing to replace existing versioned image: $target" >&2
  exit 2
fi

module purge >/dev/null 2>&1 || true
module load EBModules
module load singularity/3.6.3
mkdir -p "$target_dir"

cleanup() {
  rm -f "$building"
}
trap cleanup EXIT

if grep -q "^${USER}:" /etc/subuid 2>/dev/null && grep -q "^${USER}:" /etc/subgid 2>/dev/null; then
  echo "Building with the cluster's local fakeroot configuration."
  singularity build --fakeroot "$building" "$definition"
else
  echo "No local fakeroot mapping is configured for ${USER}; using the configured Sylabs remote builder."
  if ! singularity remote get-login-password >/dev/null 2>&1; then
    echo "ERROR: Singularity remote-build credentials are unavailable. Run 'singularity remote login' once, then resubmit this job." >&2
    exit 2
  fi
  singularity build --remote "$building" "$definition"
fi
singularity test "$building"
singularity exec --cleanenv "$building" Rscript -e 'cat("Seurat", as.character(packageVersion("Seurat")), "\n"); cat("DESeq2", as.character(packageVersion("DESeq2")), "\n"); cat("fgsea", as.character(packageVersion("fgsea")), "\n")'
chmod 444 "$building"
mv "$building" "$target"
trap - EXIT
sha256sum "$target" > "$target.sha256"
chmod 444 "$target.sha256"

echo "Installed shared Seurat image: $target"
echo "Recorded image checksum: $target.sha256"
