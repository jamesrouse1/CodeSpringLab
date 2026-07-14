# CodeSpringLab CUT&RUN

This module is the first CodeSpringLab implementation for CUT&RUN/CUT&Tag-style paired-end chromatin profiling.

## Current workflow

1. Create or provide `design_matrix.txt`
2. Optional FASTQ copy into `csl_results/<project>/data/fastq`
3. FastQC
4. Adapter trimming with cutadapt
5. Bowtie2 alignment and CUT&RUN fragment files
6. Peak calling with SEACR or MACS2
7. Peak QC with consensus SEACR peaks, consensus peak counts, and FRiP summaries

## Design matrix

The design matrix should contain at least:

```text
sample	include	target	condition	replicate	control_sample	filename
```

- `sample`: short sample name used for output folders.
- `include`: `TRUE` or `FALSE`.
- `target`: antibody/mark, for example `H3K27ac`, `H3K4me3`, `Creb`, `pCreb`, `IgG`.
- `condition`: treatment or phenotype.
- `replicate`: biological replicate.
- `control_sample`: matching IgG/control sample name for peak calling.
- `filename`: R1 FASTQ filename, or `R1,R2` pair.

`bulkCUTRUNseq.create_design_matrix_from_fastqs()` can draft this table from filenames and leave metadata editable.

## Defaults

- paired-end Bowtie2 alignment
- MAPQ >= 30
- max fragment length 1000 bp
- mitochondrial reads removed from fragment signal
- target duplicate reads kept by default
- IgG/control duplicate reads deduplicated by default
- CPM-normalized fragment bedGraph and bigWig tracks are written by default
- optional spike-in alignment/normalization can be enabled with a spike-in Bowtie2 index prefix
- SEACR uses the selected normalized paired-fragment bedGraph when available
- MACS2 uses `BAMPE` and `--keep-dup all`

These defaults follow common CUT&RUN conventions, where identical fragment starts can reflect real antibody-tethered cleavage rather than PCR duplicates, but IgG controls often have higher nonspecific duplication.

Spike-in normalization is off by default. If E. coli or another spike-in genome was included experimentally and indexed with Bowtie2, pass `normalization_mode="spikein"` and `spikein_index_path="/path/to/index/prefix"` to `bowtie2_Prep()`. If spike-in reads are very low, use CPM or no normalization instead.

## Best-practice references used for this first pass

- nf-core/cutandrun output documentation: preprocessing, Bowtie2 alignment, quality filtering, duplicate strategy, fragment QC, SEACR/MACS2 peak calling, FRiP, and fragment-length QC.
- Fred Hutch SEACR documentation: SEACR is designed for sparse CUT&RUN/chromatin profiling data and expects paired-end fragment bedGraph input.
- Meers, Tenenbaum, and Henikoff 2019: SEACR method paper for sparse enrichment analysis in CUT&RUN.

## SEACR setup

The server does not currently expose a `seacr` module. Install the SEACR script once from the CodeSpringLab folder:

```bash
bash ../scripts_DoNotTouch/SEACR/download_seacr.sh
```

Then use the `seacr_Prep()` and `seacr_RunPeakCalling()` functions.

## Minimal notebook-style run

```python
import sys
sys.path.append("../scripts_DoNotTouch/")
import make_config as mc
mc.config("cutrun")
import bulkCUTRUNseq as csl

read_path_original, read_path_destination, scriptpath_copy, genome, pairing, inpath_design = csl.filetransfer_Prep()

# Optional if reads need to be copied
# jobid = csl.filetransfer_Copy(read_path_original, scriptpath_copy)

adapter, adapter2, minlen, read1, read2, trim1, trim2, out_cutadapt, script_cutadapt = csl.cutadapt_Prep(read_path_destination, pairing)
jobid = csl.cutadapt_RunTrimming(adapter, adapter2, minlen, read1, read2, trim1, trim2, script_cutadapt)

bt = csl.bowtie2_Prep(genome=genome, pairing=pairing, inpath_design=inpath_design)
jobid = csl.bowtie2_RunAlignment(*bt)

script_seacr, seacr_table = csl.seacr_Prep(inpath_design=inpath_design)
jobid = csl.seacr_RunPeakCalling(script_seacr, seacr_table)

peakqc = csl.peakqc_Prep()
jobid = csl.peakqc_Run(*peakqc)
```
