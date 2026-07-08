import importlib.util
import inspect
import os
import re


CONFIG_KEYS = [
    "analysis_type", "project_name", "parameters_exist", "results_directory",
    "read_path_original", "read_path_destination", "genome", "genome_version", "pairing",
    "inpath_design", "scriptpath_listdir", "scriptpath_copy",
    "feature", "featurecounts_feature", "outpath_counts", "outpath_deseq2",
    "refcond", "compared", "redundant", "geneset",
    "visualizer_data_dir",
    "out_dir_star", "out_dir_kallisto", "out_dir_featurecounts", "out_dir_rsem",
    "out_dir_bowtie2", "out_peak_macs2", "outpath_diffbind", "outpath_homer",
    "tracks_dir", "qval", "removeDup"
]

ANALYSIS_LABELS = {
    "rna": "RNA-seq",
    "atac": "ATAC-seq",
    "chip": "ChIP-seq"
}

REQUIRED_SETUP_KEYS = [
    "read_path_original", "read_path_destination", "genome", "pairing",
    "inpath_design", "scriptpath_listdir", "scriptpath_copy"
]

def setup_is_complete(values):
    return all(len(str(values.get(key, "")).strip()) > 0 for key in REQUIRED_SETUP_KEYS)

def with_slash(path):
    path = os.path.expanduser(str(path).strip())
    if len(path) == 0:
        return path
    return path.rstrip("/")+"/"

def design_matrix_path(path):
    path = os.path.expanduser(str(path).strip())
    if len(path) == 0:
        return ""
    if path.endswith("design_matrix.txt"):
        return path
    return os.path.join(path.rstrip("/"), "design_matrix.txt")

def infer_pairing_from_design(inpath_design):
    design_path = design_matrix_path(inpath_design)
    if not os.path.exists(design_path):
        return ""
    try:
        with open(design_path) as handle:
            header = handle.readline()
            rows = [line.rstrip("\n").split("\t") for line in handle if len(line.strip()) > 0]
        if not rows:
            return ""
        filenames = [row[-1] for row in rows if len(row) > 0]
        joined = " ".join(filenames)
        return "y" if ("_R2_" in joined or "_R2." in joined or ",") else "n"
    except Exception:
        return ""

def first_existing_dir(paths):
    for path in paths:
        if os.path.isdir(path):
            return path
    return ""

def first_existing_design_dir(paths):
    for path in paths:
        if os.path.exists(design_matrix_path(path)):
            return os.path.dirname(design_matrix_path(path))
    return ""

def infer_standard_project_values(values, analysis_type=None):
    values = dict(values)
    analysis_type = analysis_type or values.get("analysis_type") or infer_analysis_type("rna")
    project_name = str(values.get("project_name", "")).strip()
    if len(project_name) == 0:
        return values

    res_dir = with_slash(values.get("results_directory", "../../csl_results/"))
    values["results_directory"] = res_dir
    data_dir = os.path.join(res_dir, project_name, "data")

    if os.path.isdir(data_dir) and len(str(values.get("visualizer_data_dir", "")).strip()) == 0:
        values["visualizer_data_dir"] = data_dir

    fastq_dir = os.path.join(data_dir, "fastq")
    if len(str(values.get("read_path_destination", "")).strip()) == 0 and os.path.isdir(fastq_dir):
        values["read_path_destination"] = fastq_dir
    if len(str(values.get("read_path_original", "")).strip()) == 0 and len(str(values.get("read_path_destination", "")).strip()) > 0:
        values["read_path_original"] = values["read_path_destination"]

    design_dir = first_existing_design_dir([
        os.path.join(data_dir, "manifest"),
        os.path.join(data_dir, "design_matrix"),
        os.path.join(data_dir, "design")
    ])
    if len(str(values.get("inpath_design", "")).strip()) == 0 and len(design_dir) > 0:
        values["inpath_design"] = design_dir

    if len(str(values.get("pairing", "")).strip()) == 0 and len(str(values.get("inpath_design", "")).strip()) > 0:
        inferred_pairing = infer_pairing_from_design(values["inpath_design"])
        if len(inferred_pairing) > 0:
            values["pairing"] = inferred_pairing

    values.setdefault("scriptpath_listdir", "../scripts_DoNotTouch/fastq/qsub_listdir.sh")
    values.setdefault("scriptpath_copy", "../scripts_DoNotTouch/fastq/qsub_copy.sh")
    if len(str(values.get("scriptpath_listdir", "")).strip()) == 0:
        values["scriptpath_listdir"] = "../scripts_DoNotTouch/fastq/qsub_listdir.sh"
    if len(str(values.get("scriptpath_copy", "")).strip()) == 0:
        values["scriptpath_copy"] = "../scripts_DoNotTouch/fastq/qsub_copy.sh"

    if analysis_type == "rna":
        standard_dirs = {
            "out_dir_star": os.path.join(data_dir, "star"),
            "out_dir_kallisto": os.path.join(data_dir, "kallisto"),
            "out_dir_featurecounts": os.path.join(data_dir, "featurecounts"),
            "out_dir_rsem": os.path.join(data_dir, "rsem"),
            "outpath_counts": os.path.join(data_dir, "counts"),
            "outpath_deseq2": os.path.join(data_dir, "deseq2")
        }
        for key, path in standard_dirs.items():
            if len(str(values.get(key, "")).strip()) == 0 and os.path.isdir(path):
                values[key] = path
    elif analysis_type == "atac":
        standard_dirs = {
            "out_dir_bowtie2": os.path.join(data_dir, "bowtie2"),
            "out_peak_macs2": os.path.join(data_dir, "macs2"),
            "outpath_diffbind": os.path.join(data_dir, "diffbind"),
            "outpath_homer": os.path.join(data_dir, "homer"),
            "tracks_dir": os.path.join(data_dir, "tracks")
        }
        for key, path in standard_dirs.items():
            if len(str(values.get(key, "")).strip()) == 0 and os.path.isdir(path):
                values[key] = path

    values["parameters_exist"] = "y" if setup_is_complete(values) else "n"
    return values

def project_data_dir_exists(values):
    project_name = str(values.get("project_name", "")).strip()
    if len(project_name) == 0:
        return False
    res_dir = with_slash(values.get("results_directory", "../../csl_results/"))
    return os.path.isdir(os.path.join(res_dir, project_name, "data"))

def scripts_dir():
    return os.path.dirname(os.path.abspath(__file__))

def config_path():
    return os.path.join(scripts_dir(), "config.py")

def project_config_dir(analysis_type):
    path = os.path.join(scripts_dir(), "project_configs", analysis_type)
    os.makedirs(path, exist_ok=True)
    return path

def safe_project_name(project_name):
    project_name = str(project_name).strip()
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", project_name)

def project_config_path(analysis_type, project_name):
    return os.path.join(project_config_dir(analysis_type), safe_project_name(project_name)+".py")

def read_values(path=None):
    path = path or config_path()
    values = {}
    if os.path.exists(path):
        spec = importlib.util.spec_from_file_location("_codespring_config_values", path)
        cfg = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cfg)
        for key, value in vars(cfg).items():
            if not key.startswith("_"):
                values[key] = value
    return values

def write_values(values, path=None):
    path = path or config_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ordered = [key for key in CONFIG_KEYS if key in values]
    ordered += sorted([key for key in values if key not in ordered and not key.startswith("__")])
    with open(path, "w") as conf:
        for key in ordered:
            conf.write(key+"="+repr(str(values[key]))+"\n")

def infer_analysis_type(default=None):
    candidates = [os.getcwd()]
    for frame in inspect.stack()[1:]:
        candidates.append(frame.filename)
    joined = " ".join(candidates).lower()
    if "bulkrnaseq" in joined or "rnaseq" in joined or "rna_seq" in joined:
        return "rna"
    if "bulkatacseq" in joined or "atac" in joined:
        return "atac"
    if "bulkchipseq" in joined or "chip" in joined:
        return "chip"
    return default or "rna"

def save_config_values(values, analysis_type=None):
    analysis_type = analysis_type or values.get("analysis_type") or infer_analysis_type("rna")
    active = read_values(config_path())
    merged = dict(active)
    merged.update({key: value for key, value in values.items() if value is not None})
    merged["analysis_type"] = analysis_type

    project_name = str(merged.get("project_name", "")).strip()
    if project_name:
        merged["last_"+analysis_type+"_project_name"] = project_name
        merged[analysis_type+"_project_name"] = project_name

    write_values(merged, config_path())
    if project_name:
        project_values = {key: merged[key] for key in CONFIG_KEYS if key in merged}
        project_values["analysis_type"] = analysis_type
        project_values["project_name"] = project_name
        write_values(project_values, project_config_path(analysis_type, project_name))
    return merged

def load_project_values(analysis_type, project_name):
    return read_values(project_config_path(analysis_type, project_name))

def activate_project(analysis_type, project_name, results_directory=None):
    active = read_values(config_path())
    project_values = load_project_values(analysis_type, project_name)
    if not project_values:
        if str(active.get("project_name", "")).strip() == str(project_name).strip():
            project_values = dict(active)
        else:
            project_values = {}
    project_values.update({
        "analysis_type": analysis_type,
        "project_name": project_name,
        "results_directory": results_directory or project_values.get("results_directory") or active.get("results_directory", "../../csl_results/")
    })
    project_values = infer_standard_project_values(project_values, analysis_type)
    return save_config_values(project_values, analysis_type)
