# process.py

import subprocess
import argparse
import json
import os

def process():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")

    options     = parser.parse_args()
    json_file   = options.json_file
    verbose     = options.verbose
    
    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return

    # remove plots if they exist
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
            command = "./condor/processResults.sh {0} {1}".format(directory, era) 
            returncode = subprocess.check_call(command, shell=True)

if __name__ == "__main__":
    process()



