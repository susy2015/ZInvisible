# run_modules.py

import argparse
import json
import os
from cutflow_plot import makeCutflows 
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from search_bins import SearchBins
from search_bins import ValidationBins

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options     = parser.parse_args()
    json_file   = options.json_file
    verbose     = options.verbose
    
    doCutflows = True
    doPhotons = True
    useNbNsvSelection = True
    draw = True

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    eras = ["2016", "2017", "2018_AB", "2018_CD"]
    plot_dir  = "more_plots"
    # add "/" to directory if not present
    if plot_dir[-1] != "/":
        plot_dir += "/"
    # make directory if it does not exist
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    latex_dir = "latex_files"
    # add "/" to directory if not present
    if latex_dir[-1] != "/":
        latex_dir += "/"
    # make directory if it does not exist
    if not os.path.exists(latex_dir):
        os.makedirs(latex_dir)
    

    N = Normalization(useNbNsvSelection, verbose)
    S = Shape(plot_dir, draw, verbose)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            if doCutflows:
                makeCutflows(result_file, era, plot_dir, doPhotons)
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)

        N.makeTexFile(latex_dir + "normalization_Zmass.tex")
        
        # search bins
        SB = SearchBins(N, S, eras, plot_dir, verbose)
        # validation bins
        VB = ValidationBins(N, S, eras, plot_dir, verbose)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            VB.getValues(result_file, era)
        
        VB.makeTexFile(latex_dir + "zinv_prediction.tex")


if __name__ == "__main__":
    main()



