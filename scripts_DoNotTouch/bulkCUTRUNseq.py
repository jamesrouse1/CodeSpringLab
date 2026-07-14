import config_store
import importlib
import os
import re
import shutil
import types
import pandas as pd
from IPython.display import Image

DEFAULT_SPIKEIN_INDEX = "/grid/bsr/data/data/utama/genome/ecoli_k12/bowtie2_index/ecoli_k12_mg1655"
DEFAULT_SPIKEIN_NAME = "ecoli"

def _load_config_module():
    try:
        cfg = importlib.import_module("config")
        return importlib.reload(cfg)
    except ImportError:
        return types.SimpleNamespace(
            project_name="example_dataset",
            parameters_exist="n",
            results_directory="../../csl_results/"
        )


config = _load_config_module()

ANALYSIS_TYPE = "cutrun"
project_name = getattr(config, "project_name", "example_dataset")
param = getattr(config, "parameters_exist", "n")
res_dir = getattr(config, "results_directory", "../../csl_results/")
CONFIG_KEYS = config_store.CONFIG_KEYS


REFERENCE_DEFAULTS = {
    "mouse": {
        "genome_version": "mouse_gencodeM39",
        "genome_index_path": "/grid/bsr/data/data/utama/genome/mouse_gencodeM39/bowtie2_index/GRCm39_gencodeM39",
        "chromsize": "/grid/bsr/data/data/utama/genome/mouse_gencodeM39/GRCm39.chrom.sizes",
        "spikein_index_path": DEFAULT_SPIKEIN_INDEX,
        "spikein_name": DEFAULT_SPIKEIN_NAME,
        "effective_genome_size": "2654621783",
        "macs2_genome_size": "mm",
        "homer_species": "mm39"
    },
    "human": {
        "genome_version": "human_gencode50",
        "genome_index_path": "/grid/bsr/data/data/utama/genome/human_gencode50/bowtie2_index/GRCh38_gencode50",
        "chromsize": "/grid/bsr/data/data/utama/genome/human_gencode50/GRCh38.chrom.sizes",
        "spikein_index_path": DEFAULT_SPIKEIN_INDEX,
        "spikein_name": DEFAULT_SPIKEIN_NAME,
        "effective_genome_size": "2913022398",
        "macs2_genome_size": "hs",
        "homer_species": "hg38"
    }
}


def _analysis_type():
    return config_store.infer_analysis_type(ANALYSIS_TYPE)


def _config_value(key, default=""):
    value = getattr(config, key, default)
    if default != "" and len(str(value).strip()) == 0:
        return default
    return value


def _as_config_dir(path):
    if path is None:
        return ""
    path = os.path.expanduser(str(path).strip())
    if len(path) == 0:
        return ""
    return path.rstrip("/")


def _with_slash(path):
    path = _as_config_dir(path)
    if len(path) == 0:
        return path
    return path + "/"


def _config_snapshot():
    values = {}
    for key in CONFIG_KEYS:
        if hasattr(config, key):
            values[key] = getattr(config, key)
    values.setdefault("analysis_type", _analysis_type())
    values.setdefault("project_name", project_name)
    values.setdefault("parameters_exist", param)
    values.setdefault("results_directory", res_dir)
    return values


def _save_config_updates(**updates):
    global project_name
    global param
    global res_dir

    values = _config_snapshot()
    for key, value in updates.items():
        if value is None:
            continue
        if isinstance(value, str):
            value = os.path.expanduser(value.strip())
        values[key] = value
    config_store.save_config_values(values, analysis_type=_analysis_type())
    importlib.invalidate_caches()
    cfg = _load_config_module()
    project_name = getattr(cfg, "project_name", project_name)
    param = getattr(cfg, "parameters_exist", param)
    res_dir = getattr(cfg, "results_directory", res_dir)
    return cfg


def _saved_or_prompt(key, prompt, default="", example=None):
    value = _config_value(key, default)
    if value is not None and len(str(value).strip()) > 0:
        value = str(value).strip()
        print("Using saved " + key + ": " + value)
        return value
    print("========================================")
    print(prompt)
    if example is not None:
        print("\033[91m" + "Example:" + "\x1b[0m")
        print("\033[94m" + example + "\x1b[0m")
    value = os.path.expanduser(input().strip())
    _save_config_updates(**{key: value})
    return value


def _parse_bool(value, default=False):
    if value is None:
        return default
    value = str(value).strip().lower()
    if len(value) == 0:
        return default
    return value in ["y", "yes", "true", "t", "1"]


def _project_dirs():
    data_dir = os.path.join(res_dir, project_name, "data")
    log_dir = os.path.join(res_dir, project_name, "log")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    return data_dir, log_dir


def _read_design(inpath_design):
    path = config_store.design_matrix_path(inpath_design)
    if not os.path.exists(path):
        raise FileNotFoundError("Could not find design_matrix.txt: " + path)
    design = pd.read_table(path)
    if design.shape[1] < 2:
        raise ValueError("design_matrix.txt needs at least sample and filename columns.")
    if "sample" not in design.columns:
        design = design.rename(columns={design.columns[0]: "sample"})
    if "filename" not in design.columns:
        design = design.rename(columns={design.columns[-1]: "filename"})
    return design


def _included_design(design):
    if "include" not in design.columns:
        return design.copy()
    return design[design["include"].astype(str).str.upper() != "FALSE"].copy()


def _column_or_default(design, column, default=""):
    if column in design.columns:
        return design[column].astype(str)
    return pd.Series([default] * len(design), index=design.index)


def _sample_prefix_from_filename(filename):
    name = os.path.basename(str(filename).split(",")[0])
    name = re.sub(r"(_S\d+)?_R[12](_001)?\.f(ast)?q\.gz$", "", name)
    name = re.sub(r"_R[12]\.f(ast)?q\.gz$", "", name)
    return name


def _infer_target(sample):
    sample = str(sample)
    parts = re.split(r"[_-]+", sample)
    known = ["IgG", "H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3", "CREB", "Creb", "pCREB", "pCreb"]
    for token in parts:
        for k in known:
            if token.lower() == k.lower():
                return k
    return ""


def _infer_condition(sample):
    sample = str(sample)
    match = re.search(r"[-_](AA|Veh|Vehicle|Control|Ctrl|Treat|Treated)\d*$", sample, flags=re.IGNORECASE)
    if match:
        value = re.sub(r"\d+$", "", match.group(1))
        if value.lower() == "veh":
            return "Veh"
        return value
    return ""


def create_design_matrix_from_fastqs(raw_fastq_dir, out_dir=None):
    """Create an editable CUT&RUN design matrix from paired FASTQ names.

    The output intentionally leaves most metadata editable. CodeSpringApp can
    later expose this as a table with include TRUE/FALSE, target, condition,
    replicate, and IgG/control assignment.
    """
    raw_fastq_dir = _as_config_dir(raw_fastq_dir)
    out_dir = _as_config_dir(out_dir or os.path.join(res_dir, project_name, "data", "manifest"))
    os.makedirs(out_dir, exist_ok=True)

    files = sorted([f for f in os.listdir(raw_fastq_dir) if re.search(r"\.f(ast)?q\.gz$", f)])
    r1_files = [f for f in files if re.search(r"_R1(_001)?\.f(ast)?q\.gz$", f)]
    rows = []
    for r1 in r1_files:
        prefix = _sample_prefix_from_filename(r1)
        r2 = re.sub(r"_R1", "_R2", r1)
        filename = r1 + "," + r2 if r2 in files else r1
        rows.append({
            "sample": prefix,
            "include": "TRUE",
            "target": _infer_target(prefix),
            "condition": _infer_condition(prefix),
            "replicate": "",
            "control_sample": "",
            "filename": filename
        })
    design = pd.DataFrame(rows)
    design_path = os.path.join(out_dir, "design_matrix.txt")
    design.to_csv(design_path, sep="\t", index=False)
    _save_config_updates(inpath_design=out_dir, read_path_original=raw_fastq_dir, pairing="y")
    print("Design matrix written to:")
    print(design_path)
    return design


def cutrun_reference(genome=None):
    genome = str(genome or _config_value("genome", "mouse")).strip().lower()
    if genome not in REFERENCE_DEFAULTS:
        raise ValueError("Supported genomes are: " + ", ".join(sorted(REFERENCE_DEFAULTS)))
    return REFERENCE_DEFAULTS[genome]


def filetransfer_Prep():
    global param
    data_dir, _ = _project_dirs()
    read_path_destination = os.path.join(data_dir, "fastq")
    os.makedirs(read_path_destination, exist_ok=True)

    if param == "n":
        print("==================================")
        print("Here's the list of available species:")
        genome_list = pd.Series(["human", "mouse"])
        print("Index")
        print(genome_list)
        print("==================================")
        print("Specify the index to the species/genome:(e.g 1)")
        genome = genome_list[int(input())]

        print("==================================")
        print("Copy the path to your original CUT&RUN FASTQ folder:")
        read_path_original = os.path.expanduser(input().strip())

        print("==================================")
        print("Do you want to copy fastq files into this project folder? (y/n)")
        copyfastq = input().strip().lower()
        if copyfastq == "n":
            read_path_destination = read_path_original

        print("==================================")
        print("Paste a design_matrix.txt folder, or leave blank to auto-create one from FASTQ filenames:")
        inpath_design = os.path.expanduser(input().strip())
        if len(inpath_design) == 0:
            create_design_matrix_from_fastqs(read_path_original, os.path.join(data_dir, "manifest"))
            inpath_design = os.path.join(data_dir, "manifest")

        design = _read_design(inpath_design)
        pairing = "y" if design["filename"].astype(str).str.contains("_R2|,").any() else "n"
        scriptpath_copy = "../scripts_DoNotTouch/fastq/qsub_copy.sh"
        scriptpath_listdir = "../scripts_DoNotTouch/fastq/qsub_listdir.sh"

        _save_config_updates(
            project_name=project_name,
            parameters_exist="y",
            results_directory=res_dir,
            read_path_original=read_path_original,
            read_path_destination=read_path_destination,
            genome=genome,
            pairing=pairing,
            inpath_design=inpath_design,
            scriptpath_listdir=scriptpath_listdir,
            scriptpath_copy=scriptpath_copy
        )
    else:
        read_path_original = _config_value("read_path_original")
        read_path_destination = _config_value("read_path_destination")
        genome = _config_value("genome")
        pairing = _config_value("pairing")
        inpath_design = _config_value("inpath_design")
        scriptpath_copy = _config_value("scriptpath_copy", "../scripts_DoNotTouch/fastq/qsub_copy.sh")
        print("Using saved CUT&RUN project setup.")

    return _with_slash(read_path_original), _with_slash(read_path_destination), scriptpath_copy, genome, pairing, _with_slash(inpath_design)


def filetransfer_Copy(read_path_original, scriptpath_copy):
    if os.path.exists(os.path.join(res_dir, project_name, "data", "fastq")):
        shutil.rmtree(os.path.join(res_dir, project_name, "data", "fastq"))
    os.makedirs(os.path.join(res_dir, project_name, "log"), exist_ok=True)
    stderr = "-e " + os.path.join(res_dir, project_name, "log", "error_copyFastq.txt")
    stdout = "-o " + os.path.join(res_dir, project_name, "log", "output_copyFastq.txt")
    dest = os.path.join(res_dir, project_name, "data", "fastq")
    command = "sbatch " + stderr + " " + stdout + " " + scriptpath_copy + " " + read_path_original + " " + dest + " " + project_name
    job = os.popen(command).read().splitlines()
    print(job[0])
    return [job[0].split(" ")[3]]


def fastqc_Prep(directory=None, use_trimmed=None):
    data_dir, _ = _project_dirs()
    directory = _as_config_dir(directory or _config_value("read_path_destination", os.path.join(data_dir, "fastq")))
    if use_trimmed is None:
        use_trimmed = os.path.isdir(os.path.join(data_dir, "cutadapt"))
    if use_trimmed:
        directory = os.path.join(data_dir, "cutadapt")
        folder_fastqc = "fastqc_cutadapt"
    else:
        folder_fastqc = "fastqc"

    readlist = pd.Series([f for f in os.listdir(directory) if re.search(r"\.f(ast)?q\.gz$", f)])
    outdir_fastqc = os.path.join(data_dir, folder_fastqc)
    os.makedirs(outdir_fastqc, exist_ok=True)
    scriptpath_fastqc = "../scripts_DoNotTouch/FastQC/qsub_fastqc.sh"
    print("FastQC results will be stored in " + outdir_fastqc + "/")
    return readlist, _with_slash(directory), _with_slash(outdir_fastqc), scriptpath_fastqc


def fastqc_RunQC(readlist, outdir_fastqc, read_path_destination, scriptpath_fastqc):
    jobid = []
    for file in readlist:
        stderr = "-e " + os.path.join(res_dir, project_name, "log", "error_fastQC.txt")
        stdout = "-o " + os.path.join(res_dir, project_name, "log", "output_fastQC.txt")
        command = "sbatch " + stderr + " " + stdout + " " + scriptpath_fastqc + " " + read_path_destination + file + " " + outdir_fastqc + "/. " + project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        jobid.append(job[0].split(" ")[3])
    return jobid


def cutadapt_Prep(directory=None, pairing=None, adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", minlen="20"):
    data_dir, _ = _project_dirs()
    directory = _as_config_dir(directory or _config_value("read_path_destination", os.path.join(data_dir, "fastq")))
    pairing = str(pairing or _config_value("pairing", "y")).strip().lower()

    design = _included_design(_read_design(_config_value("inpath_design")))
    prefixes = pd.Series(design["sample"].astype(str).values)
    filenames = design["filename"].astype(str).str.split(",").str[0].map(_sample_prefix_from_filename)

    outdir_cutadapt = os.path.join(data_dir, "cutadapt")
    os.makedirs(outdir_cutadapt, exist_ok=True)
    if pairing == "y":
        scriptpath_cutadapt = "../scripts_DoNotTouch/cutadapt_PE/qsub_cutadapt_PE.sh"
        read1_list = directory + "/" + filenames + "_R1_001.fastq.gz"
        read2_list = directory + "/" + filenames + "_R2_001.fastq.gz"
        trimmed1_list = outdir_cutadapt + "/" + prefixes + "_R1_001.fastq.gz"
        trimmed2_list = outdir_cutadapt + "/" + prefixes + "_R2_001.fastq.gz"
    else:
        scriptpath_cutadapt = "../scripts_DoNotTouch/cutadapt_SE/qsub_cutadapt_SE.sh"
        read1_list = directory + "/" + filenames + "_R1_001.fastq.gz"
        read2_list = directory + "/" + filenames + "_R2_001.fastq.gz"
        trimmed1_list = outdir_cutadapt + "/" + prefixes + "_R1_001.fastq.gz"
        trimmed2_list = outdir_cutadapt + "/" + prefixes + "_R2_001.fastq.gz"

    _save_config_updates(read_path_destination=directory, pairing=pairing)
    print("Trimmed reads will be stored in " + outdir_cutadapt + "/")
    return adapter, adapter2, str(minlen), read1_list, read2_list, trimmed1_list, trimmed2_list, _with_slash(outdir_cutadapt), scriptpath_cutadapt


def cutadapt_RunTrimming(adapter, adapter2, minlen, read1_list, read2_list, trimmed1_list, trimmed2_list, scriptpath_cutadapt):
    jobid = []
    for i in range(len(read1_list)):
        stderr = "-e " + os.path.join(res_dir, project_name, "log", "error_cutadapt.txt")
        stdout = "-o " + os.path.join(res_dir, project_name, "log", "output_cutadapt.txt")
        command = "sbatch " + stderr + " " + stdout + " " + scriptpath_cutadapt + " " + str(minlen) + " " + adapter + " " + adapter2 + " " + trimmed1_list[i] + " " + trimmed2_list[i] + " " + read1_list[i] + " " + read2_list[i] + " " + project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        jobid.append(job[0].split(" ")[3])
    return jobid


def bowtie2_Prep(genome=None, pairing=None, read_dir=None, inpath_design=None, use_trimmed=None, mapq=None, max_fragment_length=None,
                 dedup_target_reads=None, dedup_control_reads=None, remove_mitochondrial_reads=None,
                 normalization_mode=None, spikein_index_path=None, spikein_name=None, spikein_min_reads=None):
    data_dir, _ = _project_dirs()
    genome = str(genome or _config_value("genome", "mouse")).strip().lower()
    pairing = str(pairing or _config_value("pairing", "y")).strip().lower()
    if use_trimmed is None:
        use_trimmed = os.path.isdir(os.path.join(data_dir, "cutadapt"))
    read_dir = _as_config_dir(read_dir or (_config_value("read_path_destination", os.path.join(data_dir, "fastq"))))
    if use_trimmed:
        read_dir = os.path.join(data_dir, "cutadapt")
    inpath_design = _as_config_dir(inpath_design or _config_value("inpath_design"))
    mapq = str(mapq or _config_value("minimum_alignment_q_score", "30"))
    max_fragment_length = str(max_fragment_length or _config_value("max_fragment_length", "1000"))
    dedup_target_reads = str(dedup_target_reads if dedup_target_reads is not None else _config_value("dedup_target_reads", "n")).lower()
    dedup_control_reads = str(dedup_control_reads if dedup_control_reads is not None else _config_value("dedup_control_reads", "y")).lower()
    remove_mitochondrial_reads = str(remove_mitochondrial_reads if remove_mitochondrial_reads is not None else _config_value("remove_mitochondrial_reads", "y")).lower()
    normalization_mode = str(normalization_mode or _config_value("normalization_mode", _config_value("normalisation_mode", "CPM"))).strip()
    spikein_index_path = str(spikein_index_path or _config_value("spikein_index_path", _config_value("spikein_genome", DEFAULT_SPIKEIN_INDEX))).strip()
    spikein_name = str(spikein_name or _config_value("spikein_name", DEFAULT_SPIKEIN_NAME)).strip() or DEFAULT_SPIKEIN_NAME
    spikein_min_reads = str(spikein_min_reads or _config_value("spikein_min_reads", "1000")).strip() or "1000"

    ref = cutrun_reference(genome)
    design = _included_design(_read_design(inpath_design))
    design["is_control"] = _column_or_default(design, "target").str.lower().isin(["igg", "input", "control"])
    design["dedup_mode"] = design["is_control"].map(lambda x: "dedup" if x and _parse_bool(dedup_control_reads, True) else "keepdup")
    if _parse_bool(dedup_target_reads, False):
        design.loc[~design["is_control"], "dedup_mode"] = "dedup"

    prefixes = pd.Series(design["sample"].astype(str).values)
    filenames = design["filename"].astype(str).str.split(",").str[0].map(_sample_prefix_from_filename)
    out_dir = os.path.join(data_dir, "bowtie2")
    for prefix in prefixes:
        os.makedirs(os.path.join(out_dir, prefix), exist_ok=True)

    read1_list = read_dir + "/" + filenames + "_R1_001.fastq.gz"
    read2_list = read_dir + "/" + filenames + "_R2_001.fastq.gz"
    out_prefix_list = out_dir + "/" + prefixes + "/" + prefixes
    dedup_mode_list = pd.Series(design["dedup_mode"].values)
    scriptpath_bowtie2 = "../scripts_DoNotTouch/bowtie2/qsub_bowtie2_cutrun_PE.sh" if pairing == "y" else "../scripts_DoNotTouch/bowtie2/qsub_bowtie2_cutrun_SE.sh"

    _save_config_updates(
        genome=genome,
        pairing=pairing,
        read_path_destination=read_dir,
        inpath_design=inpath_design,
        out_dir_bowtie2=out_dir,
        minimum_alignment_q_score=mapq,
        max_fragment_length=max_fragment_length,
        dedup_target_reads=dedup_target_reads,
        dedup_control_reads=dedup_control_reads,
        remove_mitochondrial_reads=remove_mitochondrial_reads,
        normalization_mode=normalization_mode,
        normalisation_mode=normalization_mode,
        spikein_index_path=spikein_index_path,
        spikein_genome=spikein_index_path,
        spikein_name=spikein_name,
        spikein_min_reads=spikein_min_reads
    )
    print("Bowtie2 CUT&RUN alignment results will be stored in " + out_dir + "/")
    return ref["genome_index_path"], read1_list, read2_list, out_prefix_list, out_dir, ref["chromsize"], scriptpath_bowtie2, mapq, max_fragment_length, dedup_mode_list, remove_mitochondrial_reads, normalization_mode, spikein_index_path, spikein_name, spikein_min_reads


def bowtie2_RunAlignment(genome_index_path, read1_list, read2_list, out_prefix_list, out_dir, chromsize, scriptpath_bowtie2, mapq, max_fragment_length, dedup_mode_list, remove_mitochondrial_reads, normalization_mode="CPM", spikein_index_path=DEFAULT_SPIKEIN_INDEX, spikein_name=DEFAULT_SPIKEIN_NAME, spikein_min_reads="1000"):
    jobid = []
    for i in range(len(out_prefix_list)):
        stderr = "-e " + out_prefix_list[i] + "_cutrun_bowtie2.err"
        stdout = "-o " + out_prefix_list[i] + "_cutrun_bowtie2.out"
        command = "sbatch " + stderr + " " + stdout + " " + scriptpath_bowtie2 + " " + out_prefix_list[i] + " " + genome_index_path + " " + read1_list[i] + " " + read2_list[i] + " " + chromsize + " " + project_name + " " + str(mapq) + " " + str(max_fragment_length) + " " + str(dedup_mode_list[i]) + " " + str(remove_mitochondrial_reads) + " " + str(normalization_mode) + " " + str(spikein_index_path) + " " + str(spikein_name) + " " + str(spikein_min_reads)
        job = os.popen(command).read().splitlines()
        print(job[0])
        jobid.append(job[0].split(" ")[3])
    return jobid


def bowtie2_Summary(directory=None):
    directory = _as_config_dir(directory or _config_value("out_dir_bowtie2", os.path.join(res_dir, project_name, "data", "bowtie2")))
    rows = []
    for sample in sorted(os.listdir(directory)):
        log_file = os.path.join(directory, sample, sample + "_alignment_summary.txt")
        if os.path.exists(log_file):
            vals = {}
            with open(log_file) as handle:
                for line in handle:
                    if "\t" in line:
                        key, value = line.rstrip("\n").split("\t", 1)
                        vals[key] = value
            vals["sample"] = sample
            rows.append(vals)
    summary = pd.DataFrame(rows)
    if len(summary) > 0:
        out = os.path.join(res_dir, project_name, "data", "bowtie2_summary")
        os.makedirs(out, exist_ok=True)
        summary.to_csv(os.path.join(out, "cutrun_alignment_summary.txt"), sep="\t", index=False)
    return summary


def _cutrun_normalized_bedgraph(out_dir_bowtie2, sample):
    sample_dir = os.path.join(out_dir_bowtie2, sample)
    summary_path = os.path.join(sample_dir, sample + "_alignment_summary.txt")
    if os.path.exists(summary_path):
        with open(summary_path) as handle:
            for line in handle:
                if line.startswith("normalized_bedgraph\t"):
                    path = line.rstrip("\n").split("\t", 1)[1]
                    if os.path.exists(path):
                        return path
    for suffix in ["_fragments.spikein.bedgraph", "_fragments.CPM.bedgraph", "_fragments.raw.bedgraph"]:
        path = os.path.join(sample_dir, sample + suffix)
        if os.path.exists(path):
            return path
    return os.path.join(sample_dir, sample + "_fragments.raw.bedgraph")


def seacr_Prep(inpath_design=None, out_dir_bowtie2=None, control_column=None, seacr_norm=None, seacr_stringency=None):
    data_dir, _ = _project_dirs()
    inpath_design = _as_config_dir(inpath_design or _config_value("inpath_design"))
    out_dir_bowtie2 = _as_config_dir(out_dir_bowtie2 or _config_value("out_dir_bowtie2", os.path.join(data_dir, "bowtie2")))
    control_column = control_column or _config_value("igg_control_column", "control_sample")
    seacr_norm = seacr_norm or _config_value("seacr_norm", "norm")
    seacr_stringency = seacr_stringency or _config_value("seacr_stringency", "stringent")
    design = _included_design(_read_design(inpath_design))

    seacr_dir = os.path.join(data_dir, "seacr")
    os.makedirs(seacr_dir, exist_ok=True)
    scriptpath_seacr = "../scripts_DoNotTouch/SEACR/qsub_seacr_cutrun.sh"
    rows = []
    for _, row in design.iterrows():
        target = str(row.get("target", "")).lower()
        if target in ["igg", "input", "control"]:
            continue
        sample = str(row["sample"])
        control = str(row.get(control_column, "")).strip() if control_column in row.index else ""
        target_bdg = _cutrun_normalized_bedgraph(out_dir_bowtie2, sample)
        control_bdg = _cutrun_normalized_bedgraph(out_dir_bowtie2, control) if control else "none"
        target_fragments = os.path.join(out_dir_bowtie2, sample, sample + "_fragments.bed")
        rows.append({
            "sample": sample,
            "target_bdg": target_bdg,
            "control_bdg": control_bdg,
            "target_fragments": target_fragments,
            "out_prefix": os.path.join(seacr_dir, sample, sample),
            "norm": seacr_norm,
            "stringency": seacr_stringency
        })
        os.makedirs(os.path.join(seacr_dir, sample), exist_ok=True)
    seacr_table = pd.DataFrame(rows)
    _save_config_updates(out_peak_seacr=seacr_dir, seacr_norm=seacr_norm, seacr_stringency=seacr_stringency, igg_control_column=control_column)
    print("SEACR peaks will be stored in " + seacr_dir + "/")
    return scriptpath_seacr, seacr_table


def seacr_RunPeakCalling(scriptpath_seacr, seacr_table):
    jobid = []
    for _, row in seacr_table.iterrows():
        stderr = "-e " + row["out_prefix"] + "_seacr.err"
        stdout = "-o " + row["out_prefix"] + "_seacr.out"
        command = "sbatch " + stderr + " " + stdout + " " + scriptpath_seacr + " " + row["target_bdg"] + " " + row["control_bdg"] + " " + row["norm"] + " " + row["stringency"] + " " + row["out_prefix"] + " " + project_name + " " + row.get("target_fragments", "none")
        job = os.popen(command).read().splitlines()
        print(job[0])
        jobid.append(job[0].split(" ")[3])
    return jobid


def peakqc_Prep(out_peak_seacr=None, out_dir_bowtie2=None):
    data_dir, _ = _project_dirs()
    out_peak_seacr = _as_config_dir(out_peak_seacr or _config_value("out_peak_seacr", os.path.join(data_dir, "seacr")))
    out_dir_bowtie2 = _as_config_dir(out_dir_bowtie2 or _config_value("out_dir_bowtie2", os.path.join(data_dir, "bowtie2")))
    out_dir_peakqc = os.path.join(data_dir, "cutrun_peak_qc")
    os.makedirs(out_dir_peakqc, exist_ok=True)
    scriptpath_peakqc = "../scripts_DoNotTouch/CUTRUN/qsub_cutrun_peak_qc.sh"
    _save_config_updates(out_peak_seacr=out_peak_seacr, out_dir_bowtie2=out_dir_bowtie2, out_dir_peakqc=out_dir_peakqc)
    print("CUT&RUN peak QC will be stored in " + out_dir_peakqc + "/")
    return scriptpath_peakqc, out_peak_seacr, out_dir_bowtie2, out_dir_peakqc


def peakqc_Run(scriptpath_peakqc, out_peak_seacr, out_dir_bowtie2, out_dir_peakqc):
    stderr = "-e " + os.path.join(out_dir_peakqc, "error_cutrun_peak_qc.txt")
    stdout = "-o " + os.path.join(out_dir_peakqc, "output_cutrun_peak_qc.txt")
    command = "sbatch " + stderr + " " + stdout + " " + scriptpath_peakqc + " " + out_peak_seacr + " " + out_dir_bowtie2 + " " + out_dir_peakqc + " " + project_name
    job = os.popen(command).read().splitlines()
    if len(job):
        print(job[0])
        return [job[0].split(" ")[3]]
    return []


def macs2_Prep(genome=None, inpath_design=None, out_dir_bowtie2=None, control_column=None, qval=None, broad=None):
    data_dir, _ = _project_dirs()
    genome = str(genome or _config_value("genome", "mouse")).strip().lower()
    ref = cutrun_reference(genome)
    inpath_design = _as_config_dir(inpath_design or _config_value("inpath_design"))
    out_dir_bowtie2 = _as_config_dir(out_dir_bowtie2 or _config_value("out_dir_bowtie2", os.path.join(data_dir, "bowtie2")))
    control_column = control_column or _config_value("igg_control_column", "control_sample")
    qval = str(qval or _config_value("qval", "0.01"))
    broad = "broad" if _parse_bool(broad, False) else "narrow"
    design = _included_design(_read_design(inpath_design))

    macs2_dir = os.path.join(data_dir, "macs2")
    os.makedirs(macs2_dir, exist_ok=True)
    scriptpath_macs2 = "../scripts_DoNotTouch/MACS2/qsub_macs2_cutrun_PE.sh"
    rows = []
    for _, row in design.iterrows():
        target = str(row.get("target", "")).lower()
        if target in ["igg", "input", "control"]:
            continue
        sample = str(row["sample"])
        control = str(row.get(control_column, "")).strip() if control_column in row.index else ""
        rows.append({
            "sample": sample,
            "target_bam": os.path.join(out_dir_bowtie2, sample, sample + "Aligned.sortedByCoord.out.bam"),
            "control_bam": os.path.join(out_dir_bowtie2, control, control + "Aligned.sortedByCoord.out.bam") if control else "none",
            "outdir": os.path.join(macs2_dir, sample),
            "genome_size": ref["macs2_genome_size"],
            "qval": qval,
            "peak_type": broad
        })
        os.makedirs(os.path.join(macs2_dir, sample), exist_ok=True)
    macs_table = pd.DataFrame(rows)
    _save_config_updates(out_peak_macs2=macs2_dir, qval=qval, peakcaller="macs2")
    print("MACS2 peaks will be stored in " + macs2_dir + "/")
    return scriptpath_macs2, macs_table


def macs2_RunPeakCalling(scriptpath_macs2, macs_table):
    jobid = []
    for _, row in macs_table.iterrows():
        stderr = "-e " + os.path.join(row["outdir"], row["sample"] + "_macs2.err")
        stdout = "-o " + os.path.join(row["outdir"], row["sample"] + "_macs2.out")
        command = "sbatch " + stderr + " " + stdout + " " + scriptpath_macs2 + " " + row["sample"] + " " + row["target_bam"] + " " + row["control_bam"] + " " + row["genome_size"] + " " + row["qval"] + " " + row["peak_type"] + " " + row["outdir"] + " " + project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        jobid.append(job[0].split(" ")[3])
    return jobid


def peakcaller_Recommendation():
    print("CUT&RUN peak-calling defaults:")
    print("- SEACR is recommended for sparse CUT&RUN peak calls, especially histone marks and low background data.")
    print("- MACS2 BAMPE is useful as a familiar fallback and often works well for TF-like/narrow targets.")
    print("- Use an IgG/control bedgraph or BAM when available. If no control exists, SEACR can run with a numeric threshold, but interpret cautiously.")
    print("- Default duplicate handling keeps target duplicates and deduplicates IgG/control reads only; this can be changed.")


def ListDir(directory):
    print("Here's the list of contents:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    print(dirlist)
    return dirlist
