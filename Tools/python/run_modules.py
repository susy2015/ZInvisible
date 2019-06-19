# runModules.py

import json
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
    verbose = True
    
    plot_dir = "more_plots"
    #json_file = "run_1.json" 
    #eras = ["2016", "2017", "2018_AB", "2018_CD"]
    json_file = "run_5.json" 
    eras = ["2016"]

    #N = Normalization(useAllMC, verbose)
    #S = Shape(plot_dir, draw, verbose)
    N = Normalization(useAllMC, False)
    S = Shape(plot_dir, draw, False)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            if doCutflows:
                makeCutflows(result_file, era, plot_dir, doPhotons)
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)

        N.makeTexFile("normalization_Zmass.tex")
        
        # search bins
        SB = SearchBins(N, S, eras, plot_dir, verbose)
        # validation bins
        VB = ValidationBins(N, S, eras, plot_dir, verbose)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            VB.getValues(result_file, era)
        
        VB.makeTexFile("zinv_prediction_v5.tex")


if __name__ == "__main__":
    main()



