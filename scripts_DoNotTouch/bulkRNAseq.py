
import config
import config_store
import pandas as pd
import os
import time
import re
import importlib
import shlex
import socket
import subprocess
from IPython.display import IFrame,clear_output,HTML,Image
import shutil
import gzip
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame
import gseapy as gp
from gseapy import GSEA,dotplot,heatmap
import imgkit

project_name=getattr(config, 'project_name', 'example_dataset')
param=getattr(config, 'parameters_exist', 'n')
res_dir=getattr(config, 'results_directory', '../../csl_results/')


CONFIG_KEYS = config_store.CONFIG_KEYS

def _analysis_type():
    return config_store.infer_analysis_type('rna')

def _config_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.py")

def _config_value(key, default=""):
    return getattr(config, key, default)

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
    return path+"/"

def _design_matrix_path(path):
    path = _as_config_dir(path)
    if len(path) == 0:
        return ""
    if path.endswith("design_matrix.txt"):
        return path
    return os.path.join(path, "design_matrix.txt")

def _config_snapshot():
    values = {}
    for key in CONFIG_KEYS:
        if hasattr(config, key):
            values[key] = getattr(config, key)
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

    for required_key in ["project_name", "parameters_exist", "results_directory"]:
        values.setdefault(required_key, _config_value(required_key, ""))

    ordered_keys = [key for key in CONFIG_KEYS if key in values]
    ordered_keys += sorted([key for key in values if key not in ordered_keys])

    config_store.save_config_values(values, analysis_type=_analysis_type())

    importlib.invalidate_caches()
    cfg = importlib.reload(config)
    project_name = getattr(cfg, "project_name", project_name)
    param = getattr(cfg, "parameters_exist", param)
    res_dir = getattr(cfg, "results_directory", res_dir)
    return cfg

def _saved_or_prompt(key, prompt, default="", example=None, normalize_dir=False, allow_blank=False):
    value = _config_value(key, default)
    if value is not None and len(str(value).strip()) > 0:
        value = str(value).strip()
        print("Using saved "+key+": "+value)
        return _with_slash(value) if normalize_dir else value

    print("========================================")
    print(prompt)
    if example is not None:
        print("\033[91m"+"If you want to use our example dataset, use:"+"\x1b[0m")
        print("\033[94m"+example+"\x1b[0m")
    value = input().strip()
    if len(value) == 0 and not allow_blank:
        value = default
    value = os.path.expanduser(value)
    _save_config_updates(**{key: value})
    return _with_slash(value) if normalize_dir else value

def _default_data_dir():
    return os.path.abspath(os.path.expanduser(os.path.join(res_dir, project_name, "data")))


def Tree():
    
    global res_dir
    
    cmd = "tree "+res_dir
    print(os.popen(cmd).read())

def DeleteJobs(jobid):
    
    for i in range(len(jobid)):
        command = "scancel "+jobid[i]
        #command = "qdel "+jobid[i]
        jobdel = os.popen(command).read().splitlines()
        print('Deleting job ID: '+jobid[i])
        while True:
            command = "squeue -u $USER"
            #command = "qstat"
            del_status = os.popen(command).read().splitlines()
            if len(del_status) == 1:
            #if del_status == []:
                break
            else:
                del_stat = pd.DataFrame(del_status)[0][1:].str.split(expand=True)[0].str.contains(jobid[i]).any()
                if del_stat == False:
                    break
            time.sleep(5)
            
def DeleteOneJob(jobid):

    command = "scancel "+jobid[0]
    #command = "qdel "+jobid[0]
    jobdel=os.popen(command).read().splitlines()
    print('Deleting job ID: '+jobid[0])
    while True:
        command = "squeue -u $USER"
        #command = "qstat"
        del_status=os.popen(command).read().splitlines()
        if len(del_status) == 1:
        #if del_status == []:
            break
        else:
            del_stat = pd.DataFrame(del_status)[0][1:].str.split(expand=True)[0].str.contains(jobid[0]).any()
            if del_stat == False:
                break
        time.sleep(5)
            
def Qstat(jobid):
        
    print("==================================")
    starttime = time.time()
    while True:
        command = "squeue -u $USER"
        #command = "qstat"
        running_status = os.popen(command).read().splitlines()
    
        clear_output(wait=True)

        if len(running_status) == 1:
        #if running_status == []:
            print("All jobs done !" )
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "done"
        elif pd.DataFrame(running_status)[0][1:].str.split(expand=True)[0].str.contains('|'.join(jobid)).any():
        #elif pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid)).any():

            all_stat = pd.DataFrame(running_status)[0][1:].str.split(expand=True)[4]            
            #all_stat = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[4]
            id_sub = pd.DataFrame(running_status)[0][1:].str.split(expand=True)[0].str.contains('|'.join(jobid))
            #id_sub = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid))
            id_stat = all_stat[id_sub]
            status_stat = id_stat.value_counts().rename(index={'R':'Number of jobs still running = ',
                                                'PD':'Number of jobs waiting in line = ',
                                                'CG' : 'Number of jobs still finishing = ',
                                                'CD' : 'Number of jobs completed = ',
                                                'F' : 'Number of jobs failed = ',
                                                'CA' : 'Number of jobs cancelled = ',
                                                'TO' : 'Number of jobs timeout = ',
                                                'CF' : 'Number of jobs configuring = '
                                                })
            #status_stat = id_stat.value_counts().rename(index={'r':'Number of jobs still running = ',
            #                                        'qw':'Number of jobs waiting in line = ',
            #                                        't':'Number of jobs about to run = ',
            #                                        'Rt':'Number of jobs about to run = ',
            #                                        'Rr':'Number of jobs about to run = ',
            #                                        'Eqw':'Number of jobs cannot run (server error) = '
            #                                       })
            print(status_stat)
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "undone"
        else:
            print("All jobs done !" )
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "done"
    
        if runstat=="done":
            break
        time.sleep(10)

def ListDir(directory):
        
    print("Here's the list of contents:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    print(dirlist)
    
    return dirlist


def _maybe_launch_results_explorer_from_setup():

    global project_name
    global res_dir

    print("==================================")
    print("Have you already run the analysis and only want to launch the RNA-Seq Results Explorer? (y/n)")
    jump = input().strip().lower()
    if not jump.startswith("y"):
        return None

    default_data_dir = _config_value("visualizer_data_dir", _default_data_dir())
    default_data_dir = os.path.abspath(os.path.expanduser(default_data_dir))
    default_design = _design_matrix_path(_config_value("inpath_design", ""))
    if len(default_design) > 0:
        default_design = os.path.abspath(os.path.expanduser(default_design))

    data_dir = default_data_dir
    print("Using saved/default visualizer data folder:")
    print(data_dir)

    if len(default_design) > 0:
        design_matrix = default_design
        print("Using saved design matrix path:")
        print(design_matrix)
    else:
        example_design = "../scripts_DoNotTouch/test/manifest/"
        print("==================================")
        print("Copy the path to design_matrix.txt used for the visualizer:")
        print("\033[91m"+"If you want to use our example dataset, use:"+"\x1b[0m")
        print("\033[94m"+example_design+"\x1b[0m")
        design_matrix = input().strip()
        design_matrix = _design_matrix_path(design_matrix)
        design_matrix = os.path.abspath(os.path.expanduser(design_matrix)) if len(design_matrix) > 0 else design_matrix

    _save_config_updates(visualizer_data_dir=data_dir, inpath_design=os.path.dirname(design_matrix) if len(design_matrix) > 0 else "")

    shiny_dir, outpath_shiny, shiny_config_path = shiny_Prep(data_dir=data_dir, inpath_design=design_matrix)
    _ = shiny_OutsideOneLiner(shiny_dir, shiny_config_path, start_server_here=True)
    print("Results Explorer launch commands printed above. You can skip the analysis cells below unless you need to rerun them.")

    design_dir = os.path.dirname(design_matrix) if len(design_matrix) > 0 else getattr(config, "inpath_design", "")
    read_path_destination = getattr(config, "read_path_destination", os.path.join(data_dir, "fastq"))
    read_path_original = getattr(config, "read_path_original", read_path_destination)
    genome = getattr(config, "genome", "")
    pairing = getattr(config, "pairing", "")
    scriptpath_copy = getattr(config, "scriptpath_copy", "../scripts_DoNotTouch/fastq/qsub_copy.sh")
    return read_path_original, read_path_destination, scriptpath_copy, genome, pairing, design_dir


def filetransfer_Prep():
        
    global project_name
    global param
    global res_dir

    jump_results = _maybe_launch_results_explorer_from_setup()
    if jump_results is not None:
        read_path_original,read_path_destination,scriptpath_copy,genome,pairing,inpath_design = jump_results
        return _with_slash(read_path_original),_with_slash(read_path_destination),scriptpath_copy,genome,pairing,_with_slash(inpath_design)

    os.makedirs(res_dir+project_name+"/data/",exist_ok=True)
    #os.makedirs(res_dir+project_name+"/data/manifest/",exist_ok=True)
    os.makedirs(res_dir+project_name+"/log/",exist_ok=True)
    
    if os.path.exists(res_dir+project_name+"/log/output_listFastq.txt") & os.path.exists(res_dir+project_name+"/log/error_listFastq.txt"):
        os.remove(res_dir+project_name+"/log/output_listFastq.txt")
        os.remove(res_dir+project_name+"/log/error_listFastq.txt")
    if param == "n":
        read_path_destination = res_dir+project_name+"/data/fastq/"
        #print("Read files will be copied to "+res_dir+project_name+"/data/fastq/")
        print("==================================")
        print("Here's the list of available genomes:")
        genome_list = pd.Series(['human','mouse'])
        print("Index")
        print(genome_list)
        print("==================================")
        print("Specify the index to the genome:(e.g 0)")
        print("\033[91m"+"If you want to use our example dataset, type the number"+"\033[94m"+" 1"+"\x1b[0m")
        genome_index = int(input())
        genome = genome_list[genome_index]
        print("==================================")
        print("Copy the path to your original read files folder:")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/fastq/")
        read_path_original = input()
        read_path_original = os.path.expanduser(read_path_original)
        print("==================================")
        print("Do you want to copy your fastq files to your home folder? Only recommended if the files are not in your lab/home folder yet or they will be removed in the near future:(y/n)")
        copyfastq = input()
        if copyfastq == 'n':
            read_path_destination = read_path_original
            os.makedirs(res_dir+project_name+"/data/fastq/",exist_ok=True)
        else:
            print("Read files will be copied to "+res_dir+project_name+"/data/fastq/")
        print("==================================")
        print("Copy the path to design matrix folder (If it's in your home folder, type tilde sign ~):")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/manifest/")
        inpath_design = input()
        inpath_design = os.path.expanduser(inpath_design)
        print("==================================")
        print("You'll be working with "+"\033[91m"+project_name+"\x1b[0m"+" folder")
        print("Re-running any cell will overwrite exisiting outputs in folder "+"\033[91m"+project_name+"\x1b[0m")
        print("If you don't want to overwrite, please re-run this cell and specify different unique"+"\033[91m project_name")
    
        scriptpath_listdir = "../scripts_DoNotTouch/fastq/qsub_listdir.sh"
        scriptpath_copy = "../scripts_DoNotTouch/fastq/qsub_copy.sh"
    
        des=pd.read_table(inpath_design+"/design_matrix.txt")
        filename=des.iloc[:,len(des.columns)-1]
        if filename.str.contains('_R2_').any() ==True:
            pairing = 'y'
        else:
            pairing = 'n'

        _save_config_updates(project_name=project_name,
                             parameters_exist=param,
                             results_directory=res_dir,
                             read_path_original=read_path_original,
                             read_path_destination=read_path_destination,
                             genome=genome,
                             pairing=pairing,
                             inpath_design=inpath_design,
                             scriptpath_listdir=scriptpath_listdir,
                             scriptpath_copy=scriptpath_copy,
                             visualizer_data_dir=res_dir+project_name+"/data")

    else:
    
        required = ["read_path_original", "read_path_destination", "genome", "pairing", "inpath_design", "scriptpath_listdir", "scriptpath_copy"]
        missing = [key for key in required if len(str(_config_value(key, "")).strip()) == 0]
        if missing:
            print("Saved project setup is missing: "+", ".join(missing))
            print("Please answer the setup prompts once; the values will be saved for the next tools.")
            param = "n"
            return filetransfer_Prep()

        read_path_original=_config_value("read_path_original")
        read_path_destination=_config_value("read_path_destination")
        genome=_config_value("genome")
        pairing=_config_value("pairing")
        inpath_design=_config_value("inpath_design")
        scriptpath_listdir=_config_value("scriptpath_listdir")
        scriptpath_copy=_config_value("scriptpath_copy")
        print("Using saved project setup from config.py.")

    #command = "sbatch "+scriptpath_listdir+" "+read_path_original+" "+project_name
    ##command = "source "+scriptpath_listdir+" "+read_path_original+" "+project_name
    ##joblist=os.popen(command).read().splitlines()
    ##print(joblist)
    
    #return read_path_original,read_path_destination,scriptpath_copy,scriptpath_listdir,genome,pairing,inpath_design
    return _with_slash(read_path_original),_with_slash(read_path_destination),scriptpath_copy,genome,pairing,_with_slash(inpath_design)

def filetransfer_ListDir(read_path_original):
    
    global project_name
    global res_dir
    
    dirfileset = ['empty']
    
    if os.path.exists(res_dir+project_name+"/log/output_copyFastq.txt") & os.path.exists(res_dir+project_name+"/log/error_copyFastq.txt"):
        os.remove(res_dir+project_name+"/log/output_copyFastq.txt")
        os.remove(res_dir+project_name+"/log/error_copyFastq.txt")
    
    try:
        print("Here's the list of files in the original folder:")
        print("Index")
        listori = pd.read_csv(res_dir+project_name+"/log/output_listFastq.txt",header=None)
        listori = listori[listori[0].str.endswith('fastq.gz')]
        print(listori)

        dirfileset = read_path_original + listori[listori[0].str.endswith('fastq.gz')]
    except OSError:
        print("Access to view files in original directory is still pending")
    
    return dirfileset[0]

def filetransfer_PrepDirect():
    
    read_path_destination = _saved_or_prompt("read_path_destination",
                                             "Specify the path to fastq folder used for QC:",
                                             default=res_dir+project_name+"/data/fastq/")
    _save_config_updates(read_path_destination=read_path_destination)
    return _with_slash(read_path_destination)

def filetransfer_Copy(read_path_original,scriptpath_copy):
    
    global project_name
    global res_dir
    
    if os.path.exists(res_dir+project_name+"/data/fastq") :
        shutil.rmtree(res_dir+project_name+"/data/fastq")
    
    jobid = []

    stderr = "-e "+res_dir+project_name+"/log/error_copyFastq.txt"
    stdout = "-o "+res_dir+project_name+"/log/output_copyFastq.txt"
    command = "sbatch "+stderr+" "+stdout+" "+scriptpath_copy+" "+read_path_original+" "+res_dir+project_name+"/data/fastq"+" "+project_name
    #command = "source "+scriptpath_copy+" "+read_path_original+" "+res_dir+project_name+"/data/fastq"+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job[0])
    #print(job[1])
    jobid.append(job[0].split(' ')[3])
    #jobid.append(job[1].split(' ')[2])
    
    return jobid

def filetransfer_ListDest(directory):
    
    global project_name
    global res_dir
    
    print("Here's the list of contents:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    print(dirlist)
    
    dirfileset = directory + dirlist
    
    rmhidden = [shutil.rmtree(f) for f in os.listdir(res_dir+project_name+"/data/fastq") if f.startswith(".")]
    
    return dirfileset

def filetransfer_Convert(directory,inpath_design):
    
    des=pd.read_table(inpath_design+"/design_matrix.txt")
    filename=des.iloc[:,len(des.columns)-1]

    dirlist = pd.Series(os.listdir(directory))    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    dirlist.index = range(len(dirlist))
    
    for i in range(len(dirlist)):
        for j in range(len(filename)):
            if dirlist[i] in filename[j]:
                #newname=des.iloc[j,0]+re.sub(r'^.*?_R', '_R', dirlist[i])
                #os.rename(directory+dirlist[i],directory+newname)
                if dirlist.str.contains("_R1_")[i]:
                    newname=des.iloc[j,0]+re.sub(r'^.*?_R1', '_R1', dirlist[i])
                    os.rename(directory+dirlist[i],directory+newname)
                elif dirlist.str.contains("_R2_")[i]:
                    newname=des.iloc[j,0]+re.sub(r'^.*?_R2', '_R2', dirlist[i])
                    os.rename(directory+dirlist[i],directory+newname)
     
    print("Here's the list of name-converted read files:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    print(dirlist)
    
    dirfileset = directory + dirlist
    
    return dirfileset

def fastqc_Prep(directory):
    
    global project_name
    global res_dir
    
    if os.path.exists(res_dir+project_name+"/log/output_fastQC.txt") & os.path.exists(res_dir+project_name+"/log/error_fastQC.txt"):
        os.remove(res_dir+project_name+"/log/output_fastQC.txt")
        os.remove(res_dir+project_name+"/log/error_fastQC.txt")
    
    folder_fastqc = "fastqc"
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        directory = res_dir+project_name+"/data/cutadapt/"
        folder_fastqc = "fastqc_cutadapt"
    print("========================================")
    
    readlist = pd.Series(os.listdir(directory))
    readlist = readlist[readlist.str.endswith('fastq.gz')]
    
    outdir_fastqc = res_dir+project_name+"/data/"+folder_fastqc+"/"
    os.makedirs(outdir_fastqc,exist_ok=True)
    
    print("FastQC results will be stored in "+res_dir+project_name+"/data/"+folder_fastqc+"/")
    
    scriptpath_fastqc = '../scripts_DoNotTouch/FastQC/qsub_fastqc.sh'

    return readlist,directory,outdir_fastqc,scriptpath_fastqc

def fastqc_PrepDirect():
    
    read_path_destination = _saved_or_prompt("read_path_destination",
                                             "Specify the path to fastq folder used for QC:",
                                             default=res_dir+project_name+"/data/fastq/")
    _save_config_updates(read_path_destination=read_path_destination)
    return _with_slash(read_path_destination)

def fastqc_RunQC(readlist,outdir_fastqc,read_path_destination,scriptpath_fastqc):
            
    global project_name   
    global res_dir
    
    jobid = []
    for file in readlist:
        stderr = "-e "+res_dir+project_name+"/log/error_fastQC.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_fastQC.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name
        #command = "source "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid
    
def fastqc_ListDir(outdir_fastqc):
    
    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    
    DataFrame(dirlist).to_html(outdir_fastqc+'/QC_list.html')
    imgkit.from_file(outdir_fastqc+'/QC_list.html', outdir_fastqc+'/QC_list.png')
    
    #dir_html = IFrame(outdir_fastqc+'/QC_list.html', width=1000, height=800)
    #dir_html = HTML(open(outdir_fastqc+'/QC_list.html').read())
    dir_html = Image(outdir_fastqc+'/QC_list.png')
    os.remove(outdir_fastqc+'/QC_list.png')
    
    return dir_html

def fastqc_Visualization(outdir_fastqc):

    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    #print(dirlist)
    
    print("Specify index to visualize HTML file:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist[index_files]
    options = {
    'format': 'png',
    'width': 1500,  # adjust as needed
    'height': 8000  # large enough to capture full content
    }
    imgkit.from_file(outdir_fastqc+file, outdir_fastqc+file+'.png',options=options)
    
    #qc = IFrame(outdir_fastqc+file, width=1000, height=800)
    #qc = HTML(open(outdir_fastqc+file).read())
    qc = Image(outdir_fastqc+file+'.png')
    os.remove(outdir_fastqc+file+'.png')
    
    return qc

def cutadapt_Prep(directory,pairing):
    
    global project_name
    global res_dir
    
    if os.path.exists(res_dir+project_name+"/log/output_cutadapt.txt") & os.path.exists(res_dir+project_name+"/log/error_cutadapt.txt"):
        os.remove(res_dir+project_name+"/log/output_cutadapt.txt")
        os.remove(res_dir+project_name+"/log/error_cutadapt.txt")
    
    readlist = pd.Series(os.listdir(directory))
    readlist = readlist[readlist.str.endswith('fastq.gz')]
    
    prefix = readlist.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    outdir_cutadapt = res_dir+project_name+"/data/cutadapt/"
    os.makedirs(outdir_cutadapt,exist_ok=True)
    
    print("Trimmed reads results will be stored in "+res_dir+project_name+"/data/cutadapt/")
    
    if pairing == "y":
        scriptpath_cutadapt = '../scripts_DoNotTouch/cutadapt_PE/qsub_cutadapt_PE.sh'
        read1_list = directory+'/'+prefix+'_R1_001.fastq.gz'
        read2_list = directory +'/'+prefix+'_R2_001.fastq.gz'
        trimmed1_list = outdir_cutadapt+'/'+prefix+'_R1_001.fastq.gz'
        trimmed2_list = outdir_cutadapt+'/'+prefix+'_R2_001.fastq.gz'
    else:
        scriptpath_cutadapt = '../scripts_DoNotTouch/cutadapt_SE/qsub_cutadapt_SE.sh'
        read1_list = directory+'/'+prefix+'_R1_001.fastq.gz'
        read2_list = directory +'/'+prefix+'_R2_001.fastq.gz'
        trimmed1_list = outdir_cutadapt+'/'+prefix+'_R1_001.fastq.gz'
        trimmed2_list = outdir_cutadapt+'/'+prefix+'_R2_001.fastq.gz'
        
    print("========================================")
    print("Specify the adapter for R1/read1:")
    print("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA for Illumina Universal TruSeq RNA")
    print("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC for Nextera Transposase ATAC")
    print("TGGAATTCTCGG for Illumina Small RNA 3' ")
    print("GATCGTCGGACT for Illumina Small RNA 5' ")
    adapter = input()
    print("========================================")
    print("Specify the adapter for R2/read2:")
    print("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT for Illumina Universal TruSeq RNA")
    print("CTGTCTCTTATACACATCTGACGCTGCCGACGA for Nextera Transposase ATAC")
    print("TGGAATTCTCGG for Illumina Small RNA 3' ")
    print("GATCGTCGGACT for Illumina Small RNA 5' ")
    adapter2 = input()
    print("========================================")
    print("Specify minimum length after trimming (default 20):")
    minlen = input()

    return adapter,adapter2,minlen,read1_list,read2_list,trimmed1_list,trimmed2_list,outdir_cutadapt,scriptpath_cutadapt

def cutadapt_PrepDirect():
    
    read_path_destination = _saved_or_prompt("read_path_destination",
                                             "Specify the path to fastq folder used for adapter trimming:",
                                             default=res_dir+project_name+"/data/fastq/")
    _save_config_updates(read_path_destination=read_path_destination)
    return _with_slash(read_path_destination)

def cutadapt_RunTrimming(adapter,adapter2,minlen,read1_list,read2_list,trimmed1_list,trimmed2_list,scriptpath_cutadapt):
            
    global project_name
    global res_dir
    
    jobid = []
    for i in range(len(read1_list)):
        stderr = "-e "+res_dir+project_name+"/log/error_cutadapt.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_cutadapt.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        #command = "source "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def star_Prep(genome,pairing,read_dir,inpath_design):
        
    global project_name
    global res_dir

    if os.path.exists(res_dir+project_name+"/log/output_star.txt") & os.path.exists(res_dir+project_name+"/log/error_star.txt"):
        os.remove(res_dir+project_name+"/log/output_star.txt")
        os.remove(res_dir+project_name+"/log/error_star.txt")
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = res_dir+project_name+"/data/cutadapt/"
    print("========================================")
    
    des = pd.read_table(inpath_design+"/design_matrix.txt")
    prefix = des.iloc[:,0]
    prefix = pd.Series(prefix)
    prefix_filename = des.iloc[:,len(des.columns)-1]   
    prefix_filename = prefix_filename.str.replace('_R1_001.fastq.gz','',regex=False)
    prefix_filename = [x.split(",")[0] if "," in x else x for x in prefix_filename]
    prefix_filename = pd.Series(prefix_filename)
    
    #prefix = pd.Series(os.listdir(read_dir))
    #prefix = prefix[prefix.str.endswith('fastq.gz')]
    #prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir = res_dir+project_name+"/data/star/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir+prefix[i],exist_ok=True)

    print("STAR alignment results will be stored in "+res_dir+project_name+"/data/star/")
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_starindex"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_rel42_all_starindex"
    
    #read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    #read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    #out_prefix_list = out_dir+prefix+'/'+prefix
    read1_list = read_dir+'/'+prefix_filename+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix_filename+'_R2_001.fastq.gz'
    out_prefix_list = out_dir+prefix+'/'+prefix
    
    if pairing == 'y':
        scriptpath_star = '../scripts_DoNotTouch/STAR/qsub_star_PE.sh'
    else:
        scriptpath_star = '../scripts_DoNotTouch/STAR/qsub_star_SE.sh'

    _save_config_updates(genome=genome, pairing=pairing, inpath_design=inpath_design, read_path_destination=read_dir, out_dir_star=out_dir)
    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_star

def star_PrepDirect():
    
    genome = _saved_or_prompt("genome", "Specify genome:(e.g human, mouse, etc)")
    pairing = _saved_or_prompt("pairing", "Are the reads paired-end:(e.g y/n)")
    read_path_destination = _saved_or_prompt("read_path_destination",
                                             "Specify the path to fastq folder used for alignment:",
                                             default=res_dir+project_name+"/data/fastq/")
    _save_config_updates(genome=genome, pairing=pairing, read_path_destination=read_path_destination)
    return genome,pairing,_with_slash(read_path_destination)

def star_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_star):
        
    global project_name
    global res_dir
    
    jobid = []
    for i in range(len(out_prefix_list)):
        stderr = "-e "+res_dir+project_name+"/log/error_star.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_star.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_star+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        #command = "source "+scriptpath_star+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])
        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def star_ListDir(directory):
    
    global project_name
    global res_dir
    
    starlogdir = res_dir+project_name+"/data/star_summary/"
    os.makedirs(starlogdir,exist_ok=True)
    
    dirlist = pd.Series(os.listdir(directory))
    
    i = 0

    for file in dirlist:
        i += 1
        
        logfile = directory+file+"/"+file+"Log.final.out"

        log_df = pd.read_table(logfile,comment='#',header=None,sep="\t",index_col=[0])
        log_raw = log_df.rename({log_df.columns[0]:file},axis='columns')
        if i == 1 :
            log_matrix = log_raw
        else :
            log_matrix = pd.concat([log_matrix,log_raw],axis=1)

    log_matrix.to_csv(starlogdir+'summary_matrix.txt',sep='\t')        
    
    return log_matrix

def kallisto_Prep(genome,pairing,read_dir,inpath_design):
        
    global project_name
    global res_dir

    if os.path.exists(res_dir+project_name+"/log/output_kallisto.txt") & os.path.exists(res_dir+project_name+"/log/error_kallisto.txt"):
        os.remove(res_dir+project_name+"/log/output_kallisto.txt")
        os.remove(res_dir+project_name+"/log/error_kallisto.txt")
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = res_dir+project_name+"/data/cutadapt/"
    print("========================================")

    des = pd.read_table(inpath_design+"/design_matrix.txt")
    prefix = des.iloc[:,0]
    prefix = pd.Series(prefix)
    prefix_filename = des.iloc[:,len(des.columns)-1]   
    prefix_filename = prefix_filename.str.replace('_R1_001.fastq.gz','',regex=False)
    prefix_filename = [x.split(",")[0] if "," in x else x for x in prefix_filename]
    prefix_filename = pd.Series(prefix_filename)

    #prefix = pd.Series(os.listdir(read_dir))
    #prefix = prefix[prefix.str.endswith('fastq.gz')]
    #prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir_kal = res_dir+project_name+"/data/kallisto/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir_kal+prefix[i],exist_ok=True)

    print("Kallisto pseudo-alignment results will be stored in "+res_dir+project_name+"/data/kallisto/")
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.transcripts.idx"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v45.transcripts.idx"
    
    #read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    #read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    #out_prefix_list = out_dir_kal+prefix+'/'
    read1_list = read_dir+'/'+prefix_filename+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix_filename+'_R2_001.fastq.gz'
    out_prefix_list = out_dir_kal+prefix+'/'
    
    if pairing == 'y':
        scriptpath_kallisto = '../scripts_DoNotTouch/Kallisto/qsub_kallisto_PE.sh'
    else:
        scriptpath_kallisto = '../scripts_DoNotTouch/Kallisto/qsub_kallisto_SE.sh'

    _save_config_updates(genome=genome, pairing=pairing, inpath_design=inpath_design, read_path_destination=read_dir, out_dir_kallisto=out_dir_kal, feature='gene_id')
    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir_kal,scriptpath_kallisto

def kallisto_PrepDirect():
    
    genome = _saved_or_prompt("genome", "Specify genome:(e.g human, mouse, etc)")
    pairing = _saved_or_prompt("pairing", "Are the reads paired-end:(e.g y/n)")
    read_path_destination = _saved_or_prompt("read_path_destination",
                                             "Specify the path to fastq folder used for alignment:",
                                             default=res_dir+project_name+"/data/fastq/")
    _save_config_updates(genome=genome, pairing=pairing, read_path_destination=read_path_destination)
    return genome,pairing,_with_slash(read_path_destination)

def kallisto_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir_kal,scriptpath_kallisto):
        
    global project_name
    global res_dir
    
    jobid = []
    for i in range(len(out_prefix_list)):
        stderr = "-e "+res_dir+project_name+"/log/error_kallisto.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_kallisto.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_kallisto+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        #command = "source "+scriptpath_kallisto+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])
        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def featurecounts_ListDir(prefix,count_prefix_list):
    
    global project_name
    global res_dir
    
    outpath_counts = res_dir+project_name+"/data/counts/"
    os.makedirs(outpath_counts,exist_ok=True)
    print("Featurecounts summary matrix is stored in "+res_dir+project_name+"/data/counts/")
    
    i = 0

    for file in count_prefix_list:
        i += 1
        
        logfile = file+"_counts.txt.summary"

        log_df = pd.read_table(logfile,comment='#',header=None,sep="\t",index_col=[0])
        log_raw = log_df.rename({log_df.columns[0]:prefix[i-1]},axis='columns')
        if i == 1 :
            log_matrix = log_raw
        else :
            log_matrix = pd.concat([log_matrix,log_raw],axis=1)

    log_matrix.to_csv(outpath_counts+'featurecounts_summary.txt',sep='\t')        
    
    return log_matrix

def featurecounts_Prep(genome,out_dir,pairing):
    
    global project_name
    global res_dir
    
    if os.path.exists(res_dir+project_name+"/log/output_featurecounts.txt") & os.path.exists(res_dir+project_name+"/log/error_featurecounts.txt"):
        os.remove(res_dir+project_name+"/log/output_featurecounts.txt")
        os.remove(res_dir+project_name+"/log/error_featurecounts.txt")
    
    prefix = pd.Series(os.listdir(out_dir))
    
    count_dir = res_dir+project_name+"/data/featurecounts/"
   
    for i in range(len(prefix)):
        os.makedirs(count_dir+prefix[i],exist_ok=True)
    
    if genome == 'mouse':
        GTF = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf"
        strandBED = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed"
    elif genome == 'human':
        GTF = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
        strandBED = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed"
    
    if pairing == "y":
        scriptpath_featurecounts = '../scripts_DoNotTouch/featureCounts/qsub_featurecounts_PE.sh'
    else:
        scriptpath_featurecounts = '../scripts_DoNotTouch/featureCounts/qsub_featurecounts_SE.sh'

    feature = _config_value("feature", "")
    if len(str(feature).strip()) > 0:
        print("Using saved feature: "+str(feature))
    else:
        print("Specify the genomic feature to quantify (e.g gene_name, gene_id, etc):")
        feature = input()
    
    count_prefix_list = count_dir+prefix+'/'+prefix
    bam_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord.out.bam'

    _save_config_updates(genome=genome, pairing=pairing, out_dir_star=out_dir, out_dir_featurecounts=count_dir, feature=feature)
    return scriptpath_featurecounts,GTF,bam_list,count_prefix_list,prefix,feature,strandBED

def featurecounts_PrepDirect():
    
    genome = _saved_or_prompt("genome", "Specify genome:(e.g human, mouse, etc)")
    pairing = _saved_or_prompt("pairing", "Are the reads paired-end:(e.g y/n)")
    out_dir = _saved_or_prompt("out_dir_star",
                               "Specify the path to alignment folder used for quantification:",
                               default=res_dir+project_name+"/data/star/")
    _save_config_updates(genome=genome, pairing=pairing, out_dir_star=out_dir)
    return genome,pairing,_with_slash(out_dir)

def featurecounts_RunQuantification(scriptpath_featurecounts,GTF,bam_list,count_prefix_list,feature,strandBED):
     
    global project_name
    global res_dir
    
    jobid = []
    for i in range(len(bam_list)):
        stderr = "-e "+res_dir+project_name+"/log/error_featurecounts.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_featurecounts.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_featurecounts+" "+bam_list[i]+" "+GTF+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+project_name
        #command = "source "+scriptpath_featurecounts+" "+bam_list[i]+" "+GTF+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])
        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def featurecounts_CreateCountMatrix():

    global project_name
    global res_dir
    
    inpath_counts = res_dir+project_name+"/data/featurecounts/"
    outpath_counts = res_dir+project_name+"/data/counts/"
    print("Count matrix is stored in "+res_dir+project_name+"/data/counts/")

    filelist=sorted([f for f in os.listdir(inpath_counts) if not f.startswith('.')])

    i = 0

    for file in filelist:
        i += 1
        count_df = pd.read_table(os.path.join(inpath_counts,file,file+'_counts.txt'),comment='#',header=[0],index_col=[0])
        count_raw = count_df.drop(['Chr','Start','End','Strand','Length'],axis=1)
        count_raw = count_raw.rename({count_raw.columns[0]:file},axis='columns')
        if i == 1 :
            count_matrix = count_raw
        else :
            count_matrix = pd.concat([count_matrix,count_raw],axis=1)

    count_matrix.columns=count_matrix.columns.str.rstrip('_counts.txt')
    count_matrix.to_csv(outpath_counts+'count_matrix.txt',sep='\t')

    _save_config_updates(outpath_counts=outpath_counts)
    return outpath_counts,count_matrix

def rsem_Prep(genome,out_dir,pairing):
    
    global project_name
    global res_dir
    
    if os.path.exists(res_dir+project_name+"/log/output_rsem.txt") & os.path.exists(res_dir+project_name+"/log/error_rsem.txt"):
        os.remove(res_dir+project_name+"/log/output_rsem.txt")
        os.remove(res_dir+project_name+"/log/error_rsem.txt")
    
    prefix = pd.Series(os.listdir(out_dir))
    
    count_dir = res_dir+project_name+"/data/rsem/"
   
    for i in range(len(prefix)):
        os.makedirs(count_dir+prefix[i],exist_ok=True)
    
    if genome == 'mouse':
        rsem_index = "/grid/bsr/data/data/utama/genome/mouse_rsem_index_star_gencode_GRCm39_M29_v2.7.10a/mouse_gencode"
        strandBED = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed"
    elif genome == 'human':
        rsem_index = "/grid/bsr/data/data/utama/genome/human_rsem_index_star_gencode_hg38_p13_rel42_v2.7.2b/human_gencode"
        strandBED = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed"
    
    if pairing == "y":
        scriptpath_rsem = '../scripts_DoNotTouch/RSEM/qsub_RSEM_PE.sh'
    else:
        scriptpath_rsem = '../scripts_DoNotTouch/RSEM/qsub_RSEM_SE.sh'

    feature = "gene_id"
    
    count_prefix_list = count_dir+prefix+'/'+prefix
    bam_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord.out.bam'
    bamTranscript_list = out_dir+prefix+'/'+prefix+'Aligned.toTranscriptome.out.bam'

    _save_config_updates(genome=genome, pairing=pairing, out_dir_star=out_dir, out_dir_rsem=count_dir, feature=feature)
    return scriptpath_rsem,rsem_index,bam_list,count_prefix_list,prefix,feature,strandBED,bamTranscript_list

def rsem_PrepDirect():
    
    genome = _saved_or_prompt("genome", "Specify genome:(e.g human, mouse, etc)")
    pairing = _saved_or_prompt("pairing", "Are the reads paired-end:(e.g y/n)")
    out_dir = _saved_or_prompt("out_dir_star",
                               "Specify the path to alignment folder used for quantification:",
                               default=res_dir+project_name+"/data/star/")
    _save_config_updates(genome=genome, pairing=pairing, out_dir_star=out_dir)
    return genome,pairing,_with_slash(out_dir)

def rsem_RunQuantification(scriptpath_rsem,rsem_index,bam_list,count_prefix_list,feature,strandBED,bamTranscript_list):
     
    global project_name
    global res_dir
    
    jobid = []
    for i in range(len(bam_list)):
        stderr = "-e "+res_dir+project_name+"/log/error_rsem.txt"
        stdout = "-o "+res_dir+project_name+"/log/output_rsem.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_rsem+" "+bam_list[i]+" "+rsem_index+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+bamTranscript_list[i]+" "+project_name
        #command = "source "+scriptpath_rsem+" "+bam_list[i]+" "+rsem_index+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+bamTranscript_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])
        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid
    
def rsem_CreateCountMatrix():

    global project_name
    global res_dir
    
    inpath_counts = res_dir+project_name+"/data/rsem/"
    outpath_counts = res_dir+project_name+"/data/counts/"
    print("Count matrix is stored in "+res_dir+project_name+"/data/counts/")

    filelist=sorted([f for f in os.listdir(inpath_counts) if not f.startswith('.')])

    i = 0

    for file in filelist:
        i += 1
        gene_df = pd.read_table(os.path.join(inpath_counts,file,file+'.genes.results'),comment='#',header=[0],index_col=[0])
        isoform_df = pd.read_table(os.path.join(inpath_counts,file,file+'.isoforms.results'),comment='#',header=[0],index_col=[0])
        
        gene_tpm_raw = gene_df.drop(['transcript_id(s)','length','effective_length','expected_count','FPKM'],axis=1)
        gene_fpkm_raw = gene_df.drop(['transcript_id(s)','length','effective_length','expected_count','TPM'],axis=1)
        
        isoform_tpm_raw = isoform_df.drop(['length','effective_length','expected_count','IsoPct','FPKM'],axis=1)
        isoform_fpkm_raw = isoform_df.drop(['length','effective_length','expected_count','IsoPct','TPM'],axis=1)
        
        gene_tpm_raw = gene_tpm_raw.rename({gene_tpm_raw.columns[0]:file},axis='columns')
        gene_fpkm_raw = gene_fpkm_raw.rename({gene_fpkm_raw.columns[0]:file},axis='columns')
        isoform_tpm_raw = isoform_tpm_raw.rename({isoform_tpm_raw.columns[1]:file},axis='columns')
        isoform_fpkm_raw = isoform_fpkm_raw.rename({isoform_fpkm_raw.columns[1]:file},axis='columns')
        
        if i == 1 :
            gene_tpm_matrix = gene_tpm_raw
            gene_fpkm_matrix = gene_fpkm_raw
            isoform_tpm_matrix = isoform_tpm_raw
            isoform_fpkm_matrix = isoform_fpkm_raw
        else :
            gene_tpm_matrix = pd.concat([gene_tpm_matrix,gene_tpm_raw],axis=1)
            gene_fpkm_matrix = pd.concat([gene_fpkm_matrix,gene_fpkm_raw],axis=1)
            
            isoform_tpm_raw = isoform_tpm_raw.drop(['gene_id'],axis=1)
            isoform_fpkm_raw = isoform_fpkm_raw.drop(['gene_id'],axis=1)
            
            isoform_tpm_matrix = pd.concat([isoform_tpm_matrix,isoform_tpm_raw],axis=1)
            isoform_fpkm_matrix = pd.concat([isoform_fpkm_matrix,isoform_fpkm_raw],axis=1)

    gene_tpm_matrix.to_csv(outpath_counts+'gene_tpm_matrix.txt',sep='\t')
    gene_fpkm_matrix.to_csv(outpath_counts+'gene_fpkm_matrix.txt',sep='\t')
    
    isoform_tpm_matrix.to_csv(outpath_counts+'isoform_tpm_matrix.txt',sep='\t')
    isoform_fpkm_matrix.to_csv(outpath_counts+'isoform_fpkm_matrix.txt',sep='\t')

    _save_config_updates(outpath_counts=outpath_counts, feature='gene_id')
    return outpath_counts,gene_tpm_matrix,gene_fpkm_matrix,isoform_tpm_matrix,isoform_fpkm_matrix

def deseq2_Prep(inpath_design):
    
    global project_name
    global res_dir
    
    outpath = res_dir+project_name+"/data/deseq2/"
    os.makedirs(outpath,exist_ok=True)
    print("DESeq2 results are stored in "+res_dir+project_name+"/data/deseq2/")

    design = pd.read_table(inpath_design+'/design_matrix.txt',index_col=0)
    design = design.iloc[:,:len(design.columns)-1]

    print("========================================")
    print("If you have a redundant variable/column in your design matrix,")
    print("please type the name of that variable/column(e.g age)")
    print("Otherwise, leave it blank and press enter/return")
    redundant=input()
    if len(redundant) == 0:
        redundant = "NoRedundant"
    
    print("========================================")
    print("Here's the list of phenotypes/conditions/experiments")
    design_var=[]
    for i in range(len(design.columns)):
        design_var.append(design.columns[i])
        print(design.columns[i]+':')
        print(set(design.iloc[:,i]))
    
    print("========================================")
    print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g control)")
    refcond = input()
    
    print("========================================")
    print("Which phenotype/condition/replicate/batch to compare?(e.g treated)")
    compared = input()
    
    scriptpath_deseq2 = '../scripts_DoNotTouch/DESeq2/qsub_deseq2.sh'
    Rpath_deseq2 = '../scripts_DoNotTouch/DESeq2/DESeq2.R'
    
    _save_config_updates(inpath_design=inpath_design,
                         outpath_deseq2=outpath,
                         refcond=refcond,
                         compared=compared,
                         redundant=redundant)
    return scriptpath_deseq2,Rpath_deseq2,outpath,refcond,compared,design_var,redundant
    
def deseq2_PrepDirect():
    
    outpath_counts = _saved_or_prompt("outpath_counts",
                                      "Specify the path to folder containing count_matrix.txt used for DE:",
                                      default=res_dir+project_name+"/data/counts/")
    inpath_design = _saved_or_prompt("inpath_design",
                                     "Specify the path to folder containing design_matrix.txt used for DE:",
                                     example="../scripts_DoNotTouch/test/manifest/")
    _save_config_updates(outpath_counts=outpath_counts, inpath_design=inpath_design)
    return _with_slash(outpath_counts),_with_slash(inpath_design)

def deseq2_RunDE(scriptpath_deseq2,Rpath_deseq2,inpath_counts,inpath_design,outpath,refcond,compared,redundant):
    
    global project_name
    global res_dir
    
    jobid = []
    stderr = "-e "+res_dir+project_name+"/log/error_deseq2.txt"
    stdout = "-o "+res_dir+project_name+"/log/output_deseq2.txt"
    command = "sbatch "+stderr+" "+stdout+" "+scriptpath_deseq2+" "+Rpath_deseq2+" "+inpath_counts+"/count_matrix.txt"+" "+inpath_design+"/design_matrix.txt"+" "+outpath+" "+refcond+" "+compared+" "+redundant+" "+project_name
    #command = "source "+scriptpath_deseq2+" "+Rpath_deseq2+" "+inpath_counts+"/count_matrix.txt"+" "+inpath_design+"/design_matrix.txt"+" "+outpath+" "+refcond+" "+compared+" "+redundant+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job[0])
    #print(job[1])
    jobid.append(job[0].split(' ')[3])
    #jobid.append(job[1].split(' ')[2])

    return jobid

def shiny_Prep(data_dir=None, inpath_design=None):
    
    global project_name
    global res_dir

    cfg = importlib.reload(config)
    project_name = getattr(cfg, "project_name", project_name)
    res_dir = getattr(cfg, "results_directory", res_dir)

    shiny_dir = "../scripts_DoNotTouch/Shiny/"
    shiny_dir_abs = os.path.abspath(os.path.expanduser(shiny_dir))
    outpath_shiny = res_dir+project_name+"/shiny/"
    os.makedirs(outpath_shiny,exist_ok=True)

    default_data_dir = os.path.abspath(os.path.expanduser(res_dir+project_name+"/data"))
    if data_dir is None or len(str(data_dir).strip()) == 0:
        data_dir = default_data_dir
    else:
        data_dir = os.path.abspath(os.path.expanduser(data_dir))

    design_base = inpath_design if inpath_design is not None else getattr(cfg, "inpath_design", "")
    if len(design_base) > 0:
        design_base = os.path.abspath(os.path.expanduser(design_base))
        if design_base.endswith("design_matrix.txt"):
            design_matrix_path = design_base
        else:
            design_matrix_path = os.path.join(design_base, "design_matrix.txt")
    else:
        design_matrix_path = ""

    if (len(design_matrix_path) == 0) or (not os.path.exists(design_matrix_path)):
        print("========================================")
        print("Paste the full path to design_matrix.txt used for this project:")
        print("If you paste the folder instead, design_matrix.txt will be added automatically.")
        user_design = os.path.expanduser(input()).strip()
        if user_design.endswith("design_matrix.txt"):
            candidate = user_design
        else:
            candidate = os.path.join(user_design, "design_matrix.txt")
        design_matrix_path = os.path.abspath(candidate)
    else:
        design_matrix_path = os.path.abspath(design_matrix_path)

    config_path = outpath_shiny+"shiny_results_config.R"
    with open(config_path, "w") as config_file:
        config_file.write('project_name <- "'+project_name+'"'+'\n')
        config_file.write('results_root <- "'+os.path.abspath(os.path.expanduser(res_dir))+'"'+'\n')
        config_file.write('data_dir <- "'+data_dir+'"'+'\n')
        config_file.write('design_matrix_path <- "'+design_matrix_path+'"'+'\n')
        config_file.write('host <- "0.0.0.0"'+'\n')
        config_file.write('port <- 3838'+'\n')
        config_file.write('logo_search_dirs <- c('+'\n')
        config_file.write('  "'+os.path.abspath(os.path.expanduser("../scripts_DoNotTouch"))+'"'+'\n')
        config_file.write(')'+'\n')

    _save_config_updates(visualizer_data_dir=data_dir, inpath_design=os.path.dirname(design_matrix_path))

    config_path_abs = os.path.abspath(config_path)

    return shiny_dir_abs, outpath_shiny, config_path_abs


def _port_is_available(port, bind_host="127.0.0.1"):

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((bind_host, int(port)))
    except OSError:
        return False
    finally:
        sock.close()
    return True


def _pick_available_port(start_port=3838, end_port=3900):

    for port in range(int(start_port), int(end_port) + 1):
        if _port_is_available(port):
            return port
    raise RuntimeError("No free Shiny port found between {} and {}".format(start_port, end_port))


def shiny_Launch(shiny_dir, shiny_config_path, port=None, proxy_mode="relative", port_min=3838, port_max=3900, host="0.0.0.0"):

    if port is None:
        port = _pick_available_port(start_port=port_min, end_port=port_max)

    port = int(port)
    launch_cmd = [
        "bash",
        os.path.join(shiny_dir, "run_rnaseq_results_explorer.sh"),
        shiny_config_path,
        host,
        str(port)
    ]

    proc = subprocess.Popen(
        launch_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    prefix = os.environ.get("JUPYTERHUB_SERVICE_PREFIX", "/").rstrip("/")
    proxy_relative_suffix = "/proxy/{}/".format(port)
    proxy_absolute_suffix = "/proxy/absolute/{}/".format(port)
    proxy_relative_path = prefix + proxy_relative_suffix if prefix else proxy_relative_suffix
    proxy_absolute_path = prefix + proxy_absolute_suffix if prefix else proxy_absolute_suffix
    proxy_path = proxy_absolute_path if proxy_mode == "absolute" else proxy_relative_path

    return {
        "proc": proc,
        "port": port,
        "proxy_path": proxy_path,
        "proxy_path_relative": proxy_relative_path,
        "proxy_path_absolute": proxy_absolute_path,
        "launch_cmd": launch_cmd
    }


def shiny_LaunchPy(shiny_dir, shiny_config_path, port=None, proxy_mode="relative", port_min=3838, port_max=3900, host="0.0.0.0", python_bin=None):

    if port is None:
        port = _pick_available_port(start_port=port_min, end_port=port_max)

    if python_bin is None or len(str(python_bin).strip()) == 0:
        python_bin = os.environ.get("PYTHON_BIN", "")
    if python_bin is None or len(str(python_bin).strip()) == 0:
        python_bin = "python"

    port = int(port)
    launch_cmd = [
        "bash",
        os.path.join(shiny_dir, "run_rnaseq_results_explorer_py.sh"),
        shiny_config_path,
        host,
        str(port),
        str(python_bin)
    ]

    proc = subprocess.Popen(
        launch_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    prefix = os.environ.get("JUPYTERHUB_SERVICE_PREFIX", "/").rstrip("/")
    proxy_relative_suffix = "/proxy/{}/".format(port)
    proxy_absolute_suffix = "/proxy/absolute/{}/".format(port)
    proxy_relative_path = prefix + proxy_relative_suffix if prefix else proxy_relative_suffix
    proxy_absolute_path = prefix + proxy_absolute_suffix if prefix else proxy_absolute_suffix
    proxy_path = proxy_absolute_path if proxy_mode == "absolute" else proxy_relative_path

    return {
        "proc": proc,
        "port": port,
        "proxy_path": proxy_path,
        "proxy_path_relative": proxy_relative_path,
        "proxy_path_absolute": proxy_absolute_path,
        "launch_cmd": launch_cmd,
        "python_bin": python_bin
    }


def shiny_TerminalCommands(shiny_dir, shiny_config_path, username=None, server_host="bamdev1", port=None, port_min=3838, port_max=3900, server_only=False, minimal=False):

    if (not server_only) and (username is None or len(str(username).strip()) == 0):
        print("========================================")
        print("Enter your CSHL username (for SSH tunnel command):")
        username = input().strip()

    if port is None:
        port = _pick_available_port(start_port=port_min, end_port=port_max)
    port = int(port)

    run_script = os.path.join(shiny_dir, "run_rnaseq_results_explorer.sh")
    pidfile = "~/.rnaseq_shiny_{}.pid".format(port)
    logfile = "~/.rnaseq_shiny_{}.log".format(port)

    server_cmd = (
        "PORT={port}\n"
        "PIDFILE={pidfile}\n"
        "LOGFILE={logfile}\n"
        "nohup bash {script} {cfg} 0.0.0.0 $PORT > $LOGFILE 2>&1 &\n"
        "echo $! > $PIDFILE\n"
        "echo \"Started R Shiny PID $(cat $PIDFILE) on port $PORT\"\n"
        "echo \"Log: $LOGFILE\""
    ).format(
        port=port,
        pidfile=pidfile,
        logfile=logfile,
        script=run_script,
        cfg=shiny_config_path
    )

    local_cmd = None
    if not server_only:
        local_cmd = "ssh -N -L {p}:localhost:{p} {u}@{h}".format(
            p=port,
            u=username,
            h=server_host
        )

    stop_cmd = "kill $(cat {pidfile}) && rm -f {pidfile}".format(pidfile=pidfile)

    if minimal:
        print(server_cmd)
        print(stop_cmd)
    else:
        print("========================================")
        print("Run this on bamdev1 terminal:")
        print(server_cmd)
        if not server_only:
            print("========================================")
            print("Run this on your local Mac terminal:")
            print(local_cmd)
            print("========================================")
            print("Open this in browser:")
            print("http://localhost:{}/".format(port))
        print("========================================")
        print("When done, stop the server on bamdev1 with:")
        print(stop_cmd)

    return {
        "port": port,
        "server_command": server_cmd,
        "local_tunnel_command": local_cmd,
        "browser_url": "http://localhost:{}/".format(port),
        "stop_command": stop_cmd,
        "pidfile": pidfile,
        "logfile": logfile
    }


def shiny_OutsideOneLiner(shiny_dir, shiny_config_path, username=None, server_host="bamdev1", port=None, port_min=3838, port_max=3900, print_stop=True, start_server_here=True):

    if username is None or len(str(username).strip()) == 0:
        username = os.environ.get("USER", "").strip()
    if username is None or len(str(username).strip()) == 0:
        print("Enter your CSHL username:")
        username = input().strip()

    if port is None:
        port = _pick_available_port(start_port=port_min, end_port=port_max)
    port = int(port)

    run_script = os.path.join(shiny_dir, "run_rnaseq_results_explorer.sh")
    pidfile = "~/.rnaseq_shiny_{}.pid".format(port)
    logfile = "~/.rnaseq_shiny_{}.log".format(port)
    pidfile_expanded = os.path.expanduser(pidfile)
    logfile_expanded = os.path.expanduser(logfile)

    if start_server_here:
        log_handle = open(logfile_expanded, "a")
        proc = subprocess.Popen(
            ["bash", run_script, shiny_config_path, "0.0.0.0", str(port)],
            stdin=subprocess.DEVNULL,
            stdout=log_handle,
            stderr=log_handle,
            start_new_session=True
        )
        log_handle.close()
        with open(pidfile_expanded, "w") as f:
            f.write(str(proc.pid))
        one_liner = "ssh -N -L {p}:localhost:{p} {u}@{h}".format(
            p=port,
            u=username,
            h=server_host
        )
    else:
        remote_start_cmd = (
            "nohup bash {script} {cfg} 0.0.0.0 {port} > {log} 2>&1 < /dev/null & "
            "echo $! > {pid}"
        ).format(
            script=shlex.quote(run_script),
            cfg=shlex.quote(shiny_config_path),
            port=port,
            log=logfile,
            pid=pidfile
        )
        one_liner = (
            "PORT={port}; "
            "ssh {user}@{host} {remote}; "
            "ssh -N -L ${{PORT}}:localhost:${{PORT}} {user}@{host}"
        ).format(
            port=port,
            user=username,
            host=server_host,
            remote=shlex.quote(remote_start_cmd)
        )

    if start_server_here:
        stop_cmd = "PORT={}; PIDS=$(lsof -ti :$PORT); if [ -n \"$PIDS\" ]; then kill $PIDS; fi; rm -f {}".format(port, pidfile)
    else:
        stop_cmd = (
            "ssh {user}@{host} \"PORT={port}; "
            "PIDS=\\$(lsof -ti :\\$PORT); "
            "if [ -n \\\"\\$PIDS\\\" ]; then kill \\$PIDS; fi; "
            "rm -f {pid}\""
        ).format(
            user=username,
            host=server_host,
            port=port,
            pid=pidfile
        )

    print(one_liner)
    print("http://localhost:{}/".format(port))
    if print_stop:
        print(stop_cmd)

    return {
        "port": port,
        "one_liner": one_liner,
        "stop_command": stop_cmd,
        "browser_url": "http://localhost:{}/".format(port)
    }


def shiny_ServerFirstCommands(shiny_dir, shiny_config_path, username=None, server_host="bamdev1", port=None, port_min=3838, port_max=3900, print_stop=True, after_login_only=True, minimal_after_login=True):

    if (not after_login_only) and (username is None or len(str(username).strip()) == 0):
        print("Enter your CSHL username:")
        username = input().strip()

    if port is None:
        port = _pick_available_port(start_port=port_min, end_port=port_max)
    port = int(port)

    run_script = os.path.join(shiny_dir, "run_rnaseq_results_explorer.sh")
    pidfile = "~/.rnaseq_shiny_{}.pid".format(port)
    logfile = "~/.rnaseq_shiny_{}.log".format(port)

    login_cmd = None
    if username is not None and len(str(username).strip()) > 0:
        login_cmd = "ssh {user}@{host}".format(user=username, host=server_host)
    server_cmd = (
        "PORT={port}; "
        "nohup bash {script} {cfg} 0.0.0.0 $PORT > {log} 2>&1 < /dev/null & "
        "echo $! > {pid}"
    ).format(
        port=port,
        script=shlex.quote(run_script),
        cfg=shlex.quote(shiny_config_path),
        log=logfile,
        pid=pidfile
    )
    tunnel_cmd = None
    if username is not None and len(str(username).strip()) > 0:
        tunnel_cmd = "ssh -N -L {p}:localhost:{p} {u}@{h}".format(
            p=port,
            u=username,
            h=server_host
        )
    browser_url = "http://localhost:{}/".format(port)
    stop_cmd = "lsof -ti :{} | xargs -r kill".format(port)

    if after_login_only:
        if minimal_after_login:
            print("bash {} {} 0.0.0.0 {} &".format(run_script, shiny_config_path, port))
        else:
            print(server_cmd)
        if print_stop:
            print(stop_cmd)
    else:
        print(login_cmd)
        print(server_cmd)
        print(tunnel_cmd)
        print(browser_url)
        if print_stop:
            print(stop_cmd)

    return {
        "port": port,
        "login_command": login_cmd,
        "server_command": server_cmd,
        "tunnel_command": tunnel_cmd,
        "browser_url": browser_url,
        "stop_command": stop_cmd
    }


def _gseapy_cache_dir(outpath_pathway):
    cache_dir = os.path.join(os.path.dirname(outpath_pathway.rstrip('/')), "_cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir

def _packaged_mouse_human_ortholog_table():
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "reference",
        "mouse_human_orthologs_MGI.tsv"
    )

def _clean_ensembl_ids(values):
    return pd.Series(values, dtype="object").astype(str).str.replace(r"\.\d+$", "", regex=True)

def _default_gtf_paths(genome):
    genome = str(genome).lower()
    paths = []
    if genome == "human":
        env_path = os.environ.get("CSL_HUMAN_GTF", "")
        paths.extend([
            env_path,
            "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
            "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.annotation.gtf",
        ])
    elif genome == "mouse":
        env_path = os.environ.get("CSL_MOUSE_GTF", "")
        paths.extend([
            env_path,
            "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf",
            "/grid/bsr/data/data/utama/genome/mm39_gencode/gencode.vM29.annotation.gtf",
        ])
    return [p for p in paths if p and os.path.exists(os.path.expanduser(p))]

def _open_text_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def _parse_gtf_attributes(attr_text):
    out = {}
    for key, value in re.findall(r'(\S+) "([^"]*)"', attr_text):
        out[key] = value
    return out

def _read_gtf_gene_map(genome):
    for gtf_path in _default_gtf_paths(genome):
        rows = []
        try:
            with _open_text_maybe_gzip(os.path.expanduser(gtf_path)) as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2] != "gene":
                        continue
                    attrs = _parse_gtf_attributes(fields[8])
                    gene_id = attrs.get("gene_id", "")
                    gene_name = attrs.get("gene_name", "")
                    if gene_id and gene_name:
                        rows.append((gene_id.split(".")[0], gene_name))
            if rows:
                df = pd.DataFrame(rows, columns=["ensembl_gene_id", "external_gene_name"]).drop_duplicates()
                print("Using local GTF gene map: "+gtf_path)
                return df
        except Exception as exc:
            print("WARNING: Could not read GTF gene map from "+gtf_path+": "+str(exc))
    return None

def _pick_column(df, candidates):
    normalized = {re.sub(r"[^a-z0-9]+", "", str(col).lower()): col for col in df.columns}
    for candidate in candidates:
        key = re.sub(r"[^a-z0-9]+", "", candidate.lower())
        if key in normalized:
            return normalized[key]
    return None

def _filter_gsea_orthologs(orth):
    mouse_counts = orth["mouse_gene_symbol"].value_counts()
    keep = orth["mouse_gene_symbol"].map(mouse_counts).eq(1)
    filtered = orth.loc[keep, :].copy()
    removed = len(orth) - len(filtered)
    print(
        "GSEA ortholog mapping: retained "
        + str(len(filtered))
        + " one-to-one or many-mouse-to-one-human mappings and excluded "
        + str(removed)
        + " mouse-to-many/many-to-many mappings."
    )
    return filtered

def _read_mouse_human_ortholog_table(cache_dir):
    candidates = [_packaged_mouse_human_ortholog_table()]
    mouse_symbol_names = [
        "mouse_gene_symbol", "mouse_symbol", "mgi_symbol", "marker_symbol",
        "external_gene_name", "mouse_gene_name"
    ]
    mouse_id_names = [
        "mouse_ensembl_gene_id", "mouse_gene_id", "ensembl_gene_id",
        "mmusculus_ensembl_gene_id"
    ]
    human_symbol_names = [
        "human_gene_symbol", "human_symbol", "hgnc_symbol", "human_gene_name",
        "hsapiens_homolog_associated_gene_name"
    ]

    for table_path in candidates:
        if not table_path:
            continue
        table_path = os.path.expanduser(table_path)
        if not os.path.exists(table_path):
            continue
        try:
            if table_path.endswith(".csv"):
                orth = pd.read_csv(table_path)
            else:
                orth = pd.read_table(table_path)
            mouse_symbol_col = _pick_column(orth, mouse_symbol_names)
            mouse_id_col = _pick_column(orth, mouse_id_names)
            human_symbol_col = _pick_column(orth, human_symbol_names)
            if human_symbol_col is None or (mouse_symbol_col is None and mouse_id_col is None):
                print("WARNING: Ortholog table found but required columns were not recognized: "+table_path)
                continue
            out = pd.DataFrame()
            if mouse_symbol_col is not None:
                out["mouse_gene_symbol"] = orth[mouse_symbol_col]
            if mouse_id_col is not None:
                out["mouse_ensembl_gene_id"] = _clean_ensembl_ids(orth[mouse_id_col])
            out["human_gene_symbol"] = orth[human_symbol_col]
            out = out.dropna(subset=["human_gene_symbol"]).drop_duplicates()
            out["human_gene_symbol"] = out["human_gene_symbol"].astype(str).str.strip()
            if "mouse_gene_symbol" in out.columns:
                out = out.dropna(subset=["mouse_gene_symbol"])
                out["mouse_gene_symbol"] = out["mouse_gene_symbol"].astype(str).str.strip()
                out = out[out["mouse_gene_symbol"] != ""]
            out = out[out["human_gene_symbol"] != ""]
            print("Using mouse-human ortholog table: "+table_path)
            return _filter_gsea_orthologs(out)
        except Exception as exc:
            print("WARNING: Could not read mouse-human ortholog table "+table_path+": "+str(exc))
    return None

def _mouse_expression_with_symbols(gene_exp, feature):
    if feature == "gene_name":
        out = gene_exp.copy()
        out["mouse_gene_symbol"] = out.index
        return out
    gtf_map = _read_gtf_gene_map("mouse")
    if gtf_map is None:
        out = gene_exp.copy()
        out["mouse_ensembl_gene_id"] = out.index
        return out
    out = gene_exp.merge(gtf_map, how="left", left_index=True, right_on="ensembl_gene_id")
    out["mouse_ensembl_gene_id"] = out["ensembl_gene_id"]
    out["mouse_gene_symbol"] = out["external_gene_name"]
    mapped = out["mouse_gene_symbol"].notna().sum()
    print("GSEA gene labels: mapped "+str(mapped)+" mouse Ensembl IDs to mouse symbols using local GTF.")
    return out

def _map_mouse_expression_with_ortholog_table(gene_exp, feature, cache_dir):
    orth = _read_mouse_human_ortholog_table(cache_dir)
    if orth is None:
        return None
    mouse_exp = _mouse_expression_with_symbols(gene_exp, feature)
    if "mouse_ensembl_gene_id" in mouse_exp.columns and "mouse_ensembl_gene_id" in orth.columns:
        merged = mouse_exp.merge(orth, how="inner", on="mouse_ensembl_gene_id")
    elif "mouse_gene_symbol" in mouse_exp.columns and "mouse_gene_symbol" in orth.columns:
        merged = mouse_exp.merge(orth, how="inner", on="mouse_gene_symbol")
    else:
        return None
    merged["GENE"] = merged["human_gene_symbol"]
    mapped = merged["GENE"].notna().sum()
    print("GSEA gene labels: mapped "+str(mapped)+" mouse genes to human ortholog symbols using a local ortholog table.")
    return merged

def _collapse_gsea_expression(gene_exp_conv, sample_cols):
    gene_exp_conv = gene_exp_conv.dropna(subset=["GENE"])
    gene_exp_conv["GENE"] = gene_exp_conv["GENE"].astype(str).str.strip()
    gene_exp_conv = gene_exp_conv[gene_exp_conv["GENE"] != ""]
    for col in sample_cols:
        gene_exp_conv[col] = pd.to_numeric(gene_exp_conv[col], errors="coerce")
    gene_exp_conv = gene_exp_conv.dropna(subset=sample_cols, how="all")
    duplicate_gene_rows = int(gene_exp_conv.duplicated("GENE").sum())
    gene_exp_conv = gene_exp_conv.groupby("GENE", as_index=False)[sample_cols].mean()
    if duplicate_gene_rows > 0:
        print("GSEA gene labels: averaged "+str(duplicate_gene_rows)+" duplicate ortholog rows by human gene symbol.")
    gene_exp_conv.insert(1, "NAME", "na")
    return gene_exp_conv[["GENE", "NAME"] + sample_cols]

def _prepare_gseapy_expression(gene_exp, genome, feature, outpath_pathway):
    genome = str(genome).strip().lower()
    feature = str(feature).strip().lower()
    cache_dir = _gseapy_cache_dir(outpath_pathway)
    sample_cols = list(gene_exp.columns)
    gene_exp = gene_exp.copy()
    gene_exp.index = gene_exp.index.astype(str)
    if feature == "gene_id" or gene_exp.index.to_series().str.startswith("ENS").any():
        gene_exp.index = _clean_ensembl_ids(gene_exp.index).values

    if gene_exp.index.to_series().astype(str).str.startswith("ENS").any():
        feature = "gene_id"

    if genome == "human" and feature == "gene_name":
        gene_exp_conv = gene_exp.copy()
        gene_exp_conv["GENE"] = gene_exp_conv.index
        print("GSEA gene labels: using human gene symbols directly.")
        return _collapse_gsea_expression(gene_exp_conv, sample_cols)

    if genome == "human" and feature == "gene_id":
        gtf_map = _read_gtf_gene_map("human")
        if gtf_map is not None:
            gene_exp_conv = gene_exp.merge(gtf_map, how="left", left_index=True, right_on="ensembl_gene_id")
            gene_exp_conv["GENE"] = gene_exp_conv["external_gene_name"]
            mapped = gene_exp_conv["GENE"].notna().sum()
            print("GSEA gene labels: mapped "+str(mapped)+" human Ensembl IDs to symbols using local GTF.")
        else:
            gene_exp_conv = gene_exp.copy()
            gene_exp_conv["GENE"] = gene_exp_conv.index
            print("WARNING: No human GTF map found; using Ensembl IDs for GSEA.")
        return _collapse_gsea_expression(gene_exp_conv, sample_cols)

    if genome == "mouse":
        gene_exp_conv = _map_mouse_expression_with_ortholog_table(gene_exp, feature, cache_dir)
        if gene_exp_conv is not None and gene_exp_conv["GENE"].notna().sum() > 0:
            return _collapse_gsea_expression(gene_exp_conv, sample_cols)
        raise ValueError(
            "Mouse-to-human ortholog mapping requires the bundled MGI table. "
            "GTF can convert mouse Ensembl IDs to mouse symbols, but it cannot define human orthologs."
        )

    gene_exp_conv = gene_exp.copy()
    gene_exp_conv["GENE"] = gene_exp_conv.index
    print("WARNING: Genome '"+str(genome)+"' not recognized for pathway gene mapping; using input gene labels directly.")
    return _collapse_gsea_expression(gene_exp_conv, sample_cols)

def _gseapy_gene_set_prompt():
    saved = _config_value("geneset", "")
    if len(str(saved).strip()) > 0:
        print("Using saved geneset: "+str(saved))
        return os.path.expanduser(str(saved).strip())

    print("========================================")
    print("Specify gene set database or full path to a local .gmt file:")
    print("MSigDB_Hallmark_2020, KEGG_2021_Human, GO_Biological_Process_2025, Reactome_Pathways_2024")
    print("ARCHS4_TFs_Coexp, ENCODE_TF_ChIP-seq_2015, ENCODE_Histone_Modifications_2015")
    print("FANTOM6_lncRNA_KD_DEGs, miRTarBase_2017, TRANSFAC_and_JASPAR_PWMs")
    print("GTEx_Tissues_V8_2023, CellMarker_2024, Cancer_Cell_Line_Encyclopedia")
    print("ClinVar_2019, GTEx_Aging_Signatures_2021, Proteomics_Drug_Atlas_2023")
    geneset = os.path.expanduser(input().strip())
    _save_config_updates(geneset=geneset)
    return geneset
    
def _gseapy_comparison_from_config():
    refcond = _config_value("refcond", "")
    compared = _config_value("compared", "")

    if len(str(refcond).strip()) == 0:
        print("========================================")
        print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g control)")
        refcond = input().strip()
    else:
        print("Using saved reference/baseline: "+str(refcond))

    if len(str(compared).strip()) == 0:
        print("========================================")
        print("Which phenotype/condition/replicate/batch to compare?(e.g treated)")
        compared = input().strip()
    else:
        print("Using saved comparison: "+str(compared))

    _save_config_updates(refcond=refcond, compared=compared)
    return str(refcond), str(compared)

def gseapy_Prep():

    global project_name
    global res_dir

    refcond, compared = _gseapy_comparison_from_config()
    outpath_pathway = res_dir+project_name+"/data/gseapy/"+compared+'_vs_'+refcond+"/"
    os.makedirs(outpath_pathway,exist_ok=True)
    print("GSEApy results are stored in "+outpath_pathway)

    geneset = _gseapy_gene_set_prompt()
    return geneset,outpath_pathway

def gseapy_PrepDirect():

    global project_name
    global res_dir

    refcond, compared = _gseapy_comparison_from_config()
    outpath_pathway = res_dir+project_name+"/data/gseapy/"+compared+'_vs_'+refcond+"/"
    os.makedirs(outpath_pathway,exist_ok=True)
    print("GSEApy results are stored in "+outpath_pathway)

    genome = _saved_or_prompt("genome", "Specify genome:(e.g human, mouse, etc)")
    feature = _config_value("feature", "")
    if len(str(feature).strip()) == 0:
        feature = "gene_name"
        print("No saved feature found; defaulting to gene_name. Ensembl gene IDs will still be detected automatically.")
    else:
        print("Using saved feature: "+str(feature))

    inpath_design = _config_value("inpath_design", "")
    if len(str(inpath_design).strip()) == 0:
        inpath_design = _saved_or_prompt("inpath_design",
                                         "Specify the path to folder containing design_matrix.txt used for DE:",
                                         example="../scripts_DoNotTouch/test/manifest/",
                                         normalize_dir=False)
    inpath_design = _as_config_dir(inpath_design)

    outpath = _config_value("outpath_deseq2", res_dir+project_name+"/data/deseq2/")
    if len(str(outpath).strip()) > 0:
        print("Using saved/default DESeq2 normalized-counts folder: "+str(outpath))
    else:
        outpath = _saved_or_prompt("outpath_deseq2",
                                   "Specify the path to folder containing normalized counts from DE:",
                                   normalize_dir=False)
    outpath = _as_config_dir(outpath)

    geneset = _gseapy_gene_set_prompt()
    _save_config_updates(genome=genome,
                         feature=feature,
                         inpath_design=inpath_design,
                         outpath_deseq2=outpath,
                         refcond=refcond,
                         compared=compared,
                         geneset=geneset)

    return geneset,genome,feature,_with_slash(inpath_design),_with_slash(outpath),outpath_pathway,refcond,compared

def gseapy_RunPathway(geneset,genome,feature,inpath_design,outpath,outpath_pathway,refcond,compared):
    
    global project_name
    
    design = pd.read_table(inpath_design+'/design_matrix.txt',index_col=0)
    vardesign = design.T.loc[( (design.T==refcond) | (design.T==compared) ).any(axis=1),:].T.columns[0]
    #vardesign = design.T[(design.T.iloc[:,0]==refcond) | (design.T.iloc[:,0]==compared)].T.columns[0]
    design = design.loc[design[vardesign].isin([refcond,compared]),:]
    class_vector = list(design[vardesign])
    gene_exp = pd.read_table(outpath+'/normalized_counts_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
    gene_exp = gene_exp.drop(['DESCRIPTION'],axis=1,errors='ignore')
    missing_samples = [sample for sample in design.index if sample not in gene_exp.columns]
    if missing_samples:
        raise ValueError("Normalized counts file is missing samples from design matrix: "+", ".join(missing_samples))
    gene_exp = gene_exp[list(design.index)]
    gene_exp_conv = _prepare_gseapy_expression(gene_exp,genome,feature,outpath_pathway)

    if gene_exp_conv.shape[0] < 10:
        raise ValueError("Too few genes remained after pathway gene mapping ("+str(gene_exp_conv.shape[0])+"). Check genome/feature settings.")

    #############################
    
    gs = GSEA(data=gene_exp_conv,
         gene_sets=geneset,
         classes = class_vector, # cls=class_vector
         # set permutation_type to phenotype if samples >=15
         permutation_type='gene_set',
         permutation_num=1000, # reduce number to speed up test
         outdir=outpath_pathway,
         method='signal_to_noise',
         threads=4, seed= 8)
    
    gs.pheno_pos = compared
    gs.pheno_neg = refcond
    gs_res = gs.run()

    pathways = pd.read_csv(outpath_pathway+'/gseapy.gene_set.gsea.report.csv',index_col=None)
    pathways = pathways.sort_values(by=['FDR q-val'])
    pathways.to_csv(outpath_pathway+'/report.gseapy.'+geneset+'.csv',index=None)
    terms = pathways.Term

    return gs,gs_res,pathways,terms,project_name

def gseapy_DotPlot(outpath_pathway,pathways,geneset):
    
    # to save your figure, make sure that ``ofname`` is not None
    dot = dotplot(pathways,
             column="FDR q-val",
             title=geneset,
             cmap=plt.cm.YlOrRd,
             size=8,
             ofname=outpath_pathway+"DotPlot_Top10."+geneset+'.png',
             figsize=(4,5), cutoff=1)
    dot = dotplot(pathways,
             column="FDR q-val",
             title=geneset,
             cmap=plt.cm.YlOrRd,
             size=8,
             figsize=(4,5), cutoff=1)
    
    return dot

def gseapy_EnrichPlot(pathways,gs):

    terms = pathways.Term

    print("Specify index to visualize Enrichment plot of a selected pathway:(e.g 0)")
    index_files = int(input())

    enrich = gs.plot(terms[index_files]) # If choosing only 1 pathway

    return enrich

def gseapy_heatmap(pathways,gs):

    terms = pathways.Term
    
    print("Specify index to visualize Heatmap plot of a selected pathway:(e.g 0)")
    index_files = int(input())

    genes = pathways.Lead_genes[index_files].split(";")
    # Make sure that ``ofname`` is not None, if you want to save your figure to disk
    heatgsea = heatmap(df = gs.heatmat.loc[genes], z_score=0, title=terms[index_files], figsize=(14,4))

    return heatgsea

def visualization_PrepDirect():
    
    inpath_design = _saved_or_prompt("inpath_design",
                                     "Specify the path to folder containing design_matrix.txt used for DE:",
                                     example="../scripts_DoNotTouch/test/manifest/")
    outpath = _saved_or_prompt("outpath_deseq2",
                               "Specify the path to folder containing DESeq2 results:",
                               default=res_dir+project_name+"/data/deseq2/")
    refcond, compared = _gseapy_comparison_from_config()
    _save_config_updates(inpath_design=inpath_design, outpath_deseq2=outpath, refcond=refcond, compared=compared)
    return _with_slash(inpath_design),_with_slash(outpath),refcond,compared
    
def visualization_heatmap(inpath_design,outpath,refcond,compared):

    print("Provide path to folder containing genelist.txt (max 50 genes). To plot the top 50 differential genes from DESeq2 instead, type 'top50'")
    genepath = input()
    if genepath == "top50":
        genelist = pd.read_table(outpath+'/DEG_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
        genes = genelist.index[:50]
    else:
        genelist = pd.read_table(genepath+'/genelist.txt',index_col=0,header=None)
        genes = genelist.index[:50]
    
    design = pd.read_table(inpath_design+"/design_matrix.txt",index_col=0)
    design = design.iloc[:,:len(design.columns)-1]
    count_norm = pd.read_table(outpath+'/normalized_counts_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
    count_norm = count_norm.drop(['DESCRIPTION'],axis=1)
    count_norm_sig = count_norm[count_norm.index.isin(genes)]
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300 
    for i in range(len(design.columns)):
        lut = dict(zip(design.iloc[:,i].unique(), sns.color_palette("Pastel2")))
        col_colors = design.iloc[:,i].map(lut)
        heat = sns.clustermap(count_norm_sig,z_score=0,cmap='vlag',col_colors=col_colors)
        heat.savefig(outpath+"/heatmap_"+design.columns[i]+"_"+compared+'_vs_'+refcond+'(ref).png')

    return heat

def visualization_pca(design_var):

    global project_name
    global res_dir
    
    outpath = res_dir+project_name+"/data/deseq2/"
    
    print("Here's the list of phenotype/condition/replicate/batch:")
    print("Index")
    design_var = pd.DataFrame(design_var)
    print(design_var)
    print("==================================")
    print("Which index of phenotype/condition/replicate/batch to view PC plot?")
    inpca = int(input())
    
    return outpath+"/",inpca,design_var
