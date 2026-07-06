import os
import config_store


def _new_project_values(analysis_type, active_values):
    label = config_store.ANALYSIS_LABELS.get(analysis_type, analysis_type)
    keep = {key: value for key, value in active_values.items() if key.startswith("last_") or key.endswith("_project_name")}

    print("========================================")
    print("Provide a new "+label+" project name:")
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

    for key in config_store.CONFIG_KEYS:
        if key not in ["analysis_type", "project_name", "parameters_exist", "results_directory"]:
            keep[key] = ""

    keep.update({
        "analysis_type": analysis_type,
        "project_name": project_name,
        "parameters_exist": "n",
        "results_directory": res_dir
    })
    return keep

def config(analysis_type=None):
    analysis_type = analysis_type or config_store.infer_analysis_type("rna")
    label = config_store.ANALYSIS_LABELS.get(analysis_type, analysis_type)
    active_values = config_store.read_values()
    last_key = "last_"+analysis_type+"_project_name"
    last_project = active_values.get(last_key, active_values.get(analysis_type+"_project_name", ""))

    print("Do you want to use your most recent "+label+" project? (e.g y/n)")
    if len(str(last_project).strip()) > 0:
        print("Most recent "+label+" project: "+"\033[94m"+str(last_project)+"\x1b[0m")
    else:
        print("\033[91m"+"No saved "+label+" project was found. Type"+"\033[94m"+" n"+"\x1b[0m"+" to start one.")
    print("\033[91m"+"If this is your first time for this analysis type, type"+"\033[94m"+" n"+"\x1b[0m")

    param = input().strip().lower()
    if param.startswith("y") and len(str(last_project).strip()) > 0:
        values = config_store.activate_project(analysis_type, last_project)
        print("Using saved "+label+" project: "+str(values.get("project_name", last_project)))
        print("Using saved results directory: "+str(values.get("results_directory", "../../csl_results/")))
    else:
        values = _new_project_values(analysis_type, active_values)
        config_store.save_config_values(values, analysis_type)
