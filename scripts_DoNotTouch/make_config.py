import importlib.util
import os


CONFIG_KEYS = [
    "project_name", "parameters_exist", "results_directory",
    "read_path_original", "read_path_destination", "genome", "pairing",
    "inpath_design", "scriptpath_listdir", "scriptpath_copy",
    "feature", "outpath_counts", "outpath_deseq2",
    "refcond", "compared", "redundant", "geneset",
    "visualizer_data_dir"
]

def _config_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.py")

def _read_config_values():
    path = _config_path()
    values = {}
    if os.path.exists(path):
        spec = importlib.util.spec_from_file_location("_codespring_config", path)
        cfg = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cfg)
        for key in CONFIG_KEYS:
            if hasattr(cfg, key):
                values[key] = getattr(cfg, key)
    return values

def _write_config_values(values):
    ordered = [key for key in CONFIG_KEYS if key in values]
    ordered += sorted([key for key in values if key not in ordered])
    with open(_config_path(), "w") as conf:
        for key in ordered:
            conf.write(key+"="+repr(str(values[key]))+"\n")

def config():
    
    print("Do you want to use your most recent project name, genome, design matrix, and reads folders:(e.g y/n)")
    print("\033[91m"+"If this is your first time, type"+"\033[94m"+" n"+"\x1b[0m")
    param = input().strip().lower()
    values = _read_config_values()

    if param == 'n':
        print("========================================")
        print("Provide any unique project name:")
        print("\033[91m"+"If you want to use our example dataset, type"+"\033[94m"+" example_dataset"+"\x1b[0m")
        project_name = input().strip()

        print("========================================")
        print("Provide a path/location where to store your analysis results:(Sometimes your home has limited space)")
        print("\033[91m"+"If not sure, leave it blank and press enter. Results will be stored in the same location"+"\033[94m")
        res_dir = input().strip()
        if len(res_dir) == 0:
            res_dir = "../../csl_results/"
        else:
            res_dir = os.path.expanduser(res_dir).rstrip("/")+"/csl_results/"

        values["project_name"] = project_name
        values["parameters_exist"] = param
        values["results_directory"] = res_dir
        _write_config_values(values)
    else:
        values.setdefault("project_name", "example_dataset")
        values.setdefault("results_directory", "../../csl_results/")
        values["parameters_exist"] = param
        _write_config_values(values)
        print("Using saved project: "+str(values["project_name"]))
        print("Using saved results directory: "+str(values["results_directory"]))
