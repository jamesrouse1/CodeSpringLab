import importlib.util
import inspect
import os
import re


CONFIG_KEYS = [
    "analysis_type", "project_name", "parameters_exist", "results_directory",
    "read_path_original", "read_path_destination", "genome", "pairing",
    "inpath_design", "scriptpath_listdir", "scriptpath_copy",
    "feature", "outpath_counts", "outpath_deseq2",
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

def activate_project(analysis_type, project_name):
    active = read_values(config_path())
    project_values = load_project_values(analysis_type, project_name)
    if not project_values:
        project_values = {
            "analysis_type": analysis_type,
            "project_name": project_name,
            "parameters_exist": "y",
            "results_directory": active.get("results_directory", "../../csl_results/")
        }
    project_values["parameters_exist"] = "y"
    return save_config_values(project_values, analysis_type)
