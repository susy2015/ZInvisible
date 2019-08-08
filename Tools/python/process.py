# process.py

import subprocess
import argparse
import json
import os

def process():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--plot_only",    "-p", default = False, action = "store_true", help="only make plots; do not move and hadd root files, etc.")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")

    options     = parser.parse_args()
    json_file   = options.json_file
    plot_only   = options.plot_only
    verbose     = options.verbose
    
    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return

    # remove plots if they exist
    print "Removing plots"
    folder = "plots"
    for f in os.listdir(folder):
        path = os.path.join(folder, f)
        try:
            if os.path.isfile(path):
                os.unlink(path)
        except Exception as e:
            print e
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # process all eras
        for era in runMap:
            directory = runMap[era]
            print "Processing Era {0}, directory {1}".format(era, directory)
            resultFile = "condor/{0}/result.root".format(directory)
            file_exists = os.path.exists(resultFile)
            command = ""
            if plot_only:
                # Only run if file exists
                if file_exists:
                    command = "./makePlots -f -I {0} -Y {1} -R Data_MET_{1} | grep -v LHAPDF".format(resultFile, era) 
                else:
                    print "The file {0} does not exist".format(resultFile)
                    exit(1)
            else:
                # WARNING: running processResults.sh will move result.root
                # DO NOT run processResults.sh if result.root exists
                if file_exists:
                    print "WARNING: The file {0} exists. DO NOT run processResults.sh if result.root exists because processResults.sh will move result.root!".format(resultFile)
                    exit(1)
                else:
                    command = "./condor/processResults.sh {0} {1}".format(directory, era) 
            returncode = subprocess.check_call(command, shell=True)

if __name__ == "__main__":
    process()



