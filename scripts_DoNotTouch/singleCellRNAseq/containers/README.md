# Shared single-cell containers

CodeSpringApp runs every H5AD/Scanpy job inside the immutable
`codespring-scanpy_1.0.0.sif` image. The image is deliberately **not** stored
in Git: it is a large binary. A CodeSpringLab maintainer builds it once and
places it in a shared, group-readable location.

## Deploy once per cluster

Choose a shared location readable by every CodeSpringLab user, then build the
image from a node where Singularity image building is permitted:

```bash
module load singularity/3.6.3
singularity build --fakeroot /shared/codespringlab/containers/codespring-scanpy_1.0.0.sif \
  scripts_DoNotTouch/singleCellRNAseq/containers/codespring-scanpy_1.0.0.def
singularity test /shared/codespringlab/containers/codespring-scanpy_1.0.0.sif
chmod 444 /shared/codespringlab/containers/codespring-scanpy_1.0.0.sif
```

The Seurat runtime is built with a versioned batch job because compiling its R
and Bioconductor dependencies can take a while:

```bash
sbatch scripts_DoNotTouch/singleCellRNAseq/containers/build_shared_seurat.sbatch
```

By default this installs the completed, tested image at:

```text
/grid/bsr/data/data/bsr_readable_data/containers/seurat/codespring-seurat_1.0.0.sif
```

The job builds to a temporary filename, runs the definition's `%test` section,
prints the installed Seurat, DESeq2, and fgsea versions, and only then moves the
read-only image to its final name. It also records a SHA-256 checksum and
refuses to overwrite an existing version. It uses local fakeroot when the HPC
account has subordinate-user mappings; otherwise it uses the configured Sylabs
remote builder. A remote build requires a one-time `singularity remote login`.

If the cluster does not permit `--fakeroot`, build the image in a CI runner or
ask the cluster administrator to build this single SIF from the supplied
definition file. Do not put the SIF in a user's home directory.

## Connect CodeSpringApp to the shared image

Set `CSL_SCANPY_SIF` once in the shared app launcher/service configuration:

```bash
export CSL_SCANPY_SIF=/shared/codespringlab/containers/codespring-scanpy_1.0.0.sif
export CSL_SEURAT_SIF=/shared/codespringlab/containers/codespring-seurat_1.0.0.sif
```

The application then passes that immutable path to every SLURM job. Individual
users do not install Conda, Python, Scanpy, or any packages. For a small local
installation, the same image may instead be placed beside this README as
`codespring-scanpy_1.0.0.sif`; it is detected automatically.

## Updating deliberately

Create a new image filename and version when packages change (for example,
`codespring-scanpy_1.1.0.sif`). Keep old images available for reproducibility,
then change the shared `CSL_SCANPY_SIF` setting only after testing the new
image with a small H5AD project.
