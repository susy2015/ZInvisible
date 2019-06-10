# runModules.py

import json
from cutflowPlot import makeCutflows 
from calculateNormalizationFromLL import Normalization
from calculateShapeFromPhoton import Shape
from searchBins import SearchBins
from searchBins import ValidationBins

def main():
    doCutflows = False
    doPhotons = True
    useHEMVeto = True
    useAllMC = True
    draw = False
    verbose = True
    
    plot_dir = "more_plots"
    #json_file = "run_1.json" 
    #eras = ["2016", "2017", "2018_AB", "2018_CD"]
    json_file = "run_3.json" 
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
                makeCutflows(result_file, era, plot_dir, doPhotons, useHEMVeto)
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)

        N.makeTexFile("normalization_Zmass.tex")
        
        # search bins
        SB = SearchBins(N, S, eras, verbose)
        # validation bins
        VB = ValidationBins(N, S, eras, verbose)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            VB.getValues(result_file, era)
        
        VB.makeTexFile("zinv_prediction_v3.tex")


if __name__ == "__main__":
    main()



