import os
import config_store


def _results_directory_from_prompt():
    print("========================================")
    print("Provide a path/location where analysis results are stored or should be stored:(Sometimes your home has limited space)")
    print("\033[91m"+"If not sure, leave it blank and press enter. Results will be stored in the same location"+"\033[94m")
    res_dir = input().strip()
    if len(res_dir) == 0:
        return "../../csl_results/"
    if os.path.basename(os.path.expanduser(res_dir).rstrip("/")) == "csl_results":
        return os.path.expanduser(res_dir).rstrip("/")+"/"
    return os.path.expanduser(res_dir).rstrip("/")+"/csl_results/"

def _new_project_values(analysis_type, active_values):
    label = config_store.ANALYSIS_LABELS.get(analysis_type, analysis_type)
    keep = {key: value for key, value in active_values.items() if key.startswith("last_") or key.endswith("_project_name")}

    print("========================================")
    print("Provide a "+label+" project name:")
    print("\033[91m"+"If you want to use our example dataset, type"+"\033[94m"+" example_dataset"+"\x1b[0m")
    print("\033[91m"+"If this is an older existing project, type its project name here too."+"\x1b[0m")
    project_name = input().strip()
    res_dir = _results_directory_from_prompt()

    existing_values = config_store.load_project_values(analysis_type, project_name)
    candidate_values = dict(existing_values)
    candidate_values.update({
        "analysis_type": analysis_type,
        "project_name": project_name,
        "results_directory": res_dir
    })
    candidate_values = config_store.infer_standard_project_values(candidate_values, analysis_type)

    if existing_values or config_store.project_data_dir_exists(candidate_values):
        print("========================================")
        print("Using existing "+label+" project: "+str(project_name))
        print("Using results directory: "+str(candidate_values.get("results_directory", res_dir)))
        if candidate_values.get("parameters_exist", "n") == "n":
            if len(str(candidate_values.get("genome", "")).strip()) == 0:
                print("No saved genome was found for this existing project.")
                print("Specify genome index: 0 human, 1 mouse")
                genome_list = ["human", "mouse"]
                genome_index = int(input())
                candidate_values["genome"] = genome_list[genome_index]
            if len(str(candidate_values.get("inpath_design", "")).strip()) == 0:
                print("No design_matrix.txt was found under this project.")
                print("Paste the path to design_matrix.txt or its folder:")
                design_path = input().strip()
                if design_path.endswith("design_matrix.txt"):
                    candidate_values["inpath_design"] = os.path.dirname(os.path.expanduser(design_path))
                else:
                    candidate_values["inpath_design"] = os.path.expanduser(design_path)
            candidate_values = config_store.infer_standard_project_values(candidate_values, analysis_type)
        return candidate_values

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
    print("\033[91m"+"If this is your first time, or you want to use a different older project, type"+"\033[94m"+" n"+"\x1b[0m")

    param = input().strip().lower()
    if param.startswith("y") and len(str(last_project).strip()) > 0:
        values = config_store.activate_project(analysis_type, last_project)
        print("Using saved "+label+" project: "+str(values.get("project_name", last_project)))
        print("Using saved results directory: "+str(values.get("results_directory", "../../csl_results/")))
        if values.get("parameters_exist", "n") == "n":
            print("Saved setup is missing genome, design matrix, or reads paths.")
            print("Please answer the setup prompts once; they will be saved for this project.")
    else:
        values = _new_project_values(analysis_type, active_values)
        config_store.save_config_values(values, analysis_type)
