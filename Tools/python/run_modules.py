# run_modules.py

import json
import os
from cutflow_plot import makeCutflows 
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from search_bins import SearchBins
from search_bins import ValidationBins

def main():
    doCutflows = True
    doPhotons = True
    useAllMC = True
    draw = True
    verbose = False
    
    eras = ["2016", "2017", "2018_AB", "2018_CD"]
    json_file = "run_2019-06-23.json" 
    
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
    

    N = Normalization(useAllMC, verbose)
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



