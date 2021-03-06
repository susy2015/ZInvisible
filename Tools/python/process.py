# process.py

import shutil
import subprocess
import argparse
import json
import os
import re

def removeContent(folder):
    for f in os.listdir(folder):
        path = os.path.join(folder, f)
        try:
            # remove files in folder
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
            # remove subfolders
            elif os.path.isdir(path):
                shutil.rmtree(path)
        except Exception as e:
            print "Failed to delete {0}. Exception: {1}".format(path, e)

def process():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--no_plots",     "-n", default = False, action = "store_true", help="do not make plots; only move and hadd root files, etc.")
    parser.add_argument("--plot_only",    "-p", default = False, action = "store_true", help="only make plots; do not move and hadd root files, etc.")
    parser.add_argument("--run_two_only", "-r", default = False, action = "store_true", help="combine with plot_only to only plot Run 2")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")

    options         = parser.parse_args()
    json_file       = options.json_file
    no_plots        = options.no_plots
    plot_only       = options.plot_only
    run_two_only    = options.run_two_only
    verbose         = options.verbose

    # make combined Run2
    do_run2 = True
    
    # flag for not making plots
    noPlotFlag  = ""
    if no_plots:
        noPlotFlag  = "-n"

    resultFileMap = {}
    runMap = {}
    
    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return

    # remove plots if they exist
    folder = "plots"
    removeContent(folder)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # process all eras
        for era in runMap:
            # skip all eras except Run 2 if only plotting Run 2
            if plot_only and run_two_only and era != "Run2":
                continue
            directory = runMap[era]
            print "Processing Era {0}, directory {1}".format(era, directory)
            resultFile = "condor/{0}/result.root".format(directory)
            resultFileMap[era] = resultFile
            file_exists = os.path.exists(resultFile)
            command = ""
            if plot_only:
                # Only run if file exists
                if file_exists:
                    command = "./makePlots -q -f -I {0} -Y {1} -R Data_MET_{1} | grep -v LHAPDF".format(resultFile, era) 
                else:
                    print "The file {0} does not exist".format(resultFile)
                    exit(1)
            else:
                # WARNING: running processResults.sh will move result.root
                # DO NOT run processResults.sh if result.root exists
                if file_exists:
                    #print "WARNING: The file {0} exists. Skipping {1}.".format(resultFile, era)
                    #continue
                    print "WARNING: The file {0} exists.".format(resultFile)
                    print "- DO NOT run processResults.sh if result.root exists because processResults.sh will move result.root!"
                    print "- Please use the -p option to only run makePlots instead of processResults.sh."
                    #exit(1)
                    print "- Skipping {0}.".format(resultFile)
                    continue
                else:
                    command = "./condor/processResults.sh {0} {1} {2}".format(directory, era, noPlotFlag) 
            returncode = subprocess.check_call(command, shell=True)

            if not no_plots:
                # Move plots to directory for each era
                returncode = subprocess.check_call("mkdir -p {0}/{1}".format(folder, era),     shell=True)
                returncode = subprocess.check_call("mv {0}/*.png {0}/{1}".format(folder, era), shell=True)
                returncode = subprocess.check_call("mv {0}/*.pdf {0}/{1}".format(folder, era), shell=True)
    
    # --- not needed if only making plots --- #
    if do_run2 and not plot_only: 
        
        # ----------------------- # 
        # --- 2016 + 2017 and --- #
        # --- Combined Run 2  --- #
        # ----------------------- # 

        filesToCombine = {}
        filesToCombine["2016and2017"] = " ".join([resultFileMap["2016"], resultFileMap["2017"]])
        filesToCombine["Run2"]        = " ".join([resultFileMap["2016"], resultFileMap["2017"], resultFileMap["2018"]])
        
        for era in ["2016and2017", "Run2"]:
            # make directory if it does not exist
            date = re.match("runs/submission_(.*).json", json_file).group(1)
            directory = "DataMC_{0}_submission_{1}".format(era, date)
            command = "mkdir -p condor/{0}".format(directory)
            returncode = subprocess.check_call(command, shell=True)
            
            # hadd if combined file does not exist
            resultFile = "condor/{0}/result.root".format(directory)
            file_exists = os.path.exists(resultFile)
            if not file_exists:
                files = filesToCombine[era]
                command = "time ahadd.py {0} {1}".format(resultFile, files)
                returncode = subprocess.check_call(command, shell=True)

            if not no_plots:
                # run make plots 
                command = "./makePlots -q -f -I {0} -Y {1} -R Data_MET_{1} | grep -v LHAPDF".format(resultFile, era) 
                returncode = subprocess.check_call(command, shell=True)
                
                # Move plots to directory for each era
                returncode = subprocess.check_call("mkdir -p {0}/{1}".format(folder, era),     shell=True)
                returncode = subprocess.check_call("mv {0}/*.png {0}/{1}".format(folder, era), shell=True)
                returncode = subprocess.check_call("mv {0}/*.pdf {0}/{1}".format(folder, era), shell=True)

            # add this era directory to map
            runMap[era]   = directory
        
        # make Run 2 json file
        new_json_file = "runs/submission_{0}_{1}.json".format("Run2", date) 
        file_exists   = os.path.exists(new_json_file)
        if not file_exists:
            print "Creating {0}".format(new_json_file)
            with open (new_json_file, "w") as j:
                json.dump(runMap, j, sort_keys=True, indent=4, separators=(',', ' : '))


if __name__ == "__main__":
    process()


