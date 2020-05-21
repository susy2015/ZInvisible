# quick_plot.py

import ROOT
import json
import os
import argparse
import numpy as np
from tools import plot

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def run(era, result_file, verbose):
    if verbose:
        print "{0}: {1}".format(era, result_file)


def main():
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
    
    #eras = ["2016", "2017", "2018", "Run2"]
    eras = ["Run2"]
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            run(era, result_file, verbose)

if __name__ == "__main__":
    main()

