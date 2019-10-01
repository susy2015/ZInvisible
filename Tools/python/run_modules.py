# run_modules.py

import argparse
import json
import os
from cutflow_plot import makeCutflows 
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from search_bins import SearchBins
from search_bins import ValidationBins
from data_card import makeDataCard
from make_table import Table

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options     = parser.parse_args()
    json_file   = options.json_file
    verbose     = options.verbose

    doCutflows = False
    doPhotons = True
    draw = True

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    #eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM"]
    eras = ["2016"]
    dirList = []
    plot_dir                = "more_plots"
    latex_dir               = "latex_files"
    dataCard_dir            = "data_cards"
    dataCardValidation_dir  = dataCard_dir + "/validation"
    dataCardSearch_dir      = dataCard_dir + "/search"
    dirList.append(plot_dir)
    dirList.append(latex_dir)
    dirList.append(dataCard_dir)
    dirList.append(dataCardValidation_dir)
    dirList.append(dataCardSearch_dir)
    # add "/" to directory if not present
    if plot_dir[-1]  != "/":                plot_dir     += "/"
    if latex_dir[-1] != "/":                latex_dir    += "/"
    if dataCard_dir[-1] != "/":             dataCard_dir += "/"
    if dataCardValidation_dir[-1] != "/":   dataCardValidation_dir += "/"
    if dataCardSearch_dir[-1] != "/":       dataCardSearch_dir += "/"
    
    for d in dirList:
        # make directory if it does not exist
        if not os.path.exists(d):
            os.makedirs(d)

    N = Normalization(verbose)
    #S = Shape(plot_dir, draw, verbose)
    S = Shape(plot_dir, draw, True)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # validation bins
        VB = ValidationBins(N, S, eras, plot_dir, verbose)
        # search bins
        SB = SearchBins(N, S, eras, plot_dir, verbose)
        # loop over eras
        for era in eras:
            print "|---------- Era: {0} ----------|".format(era)
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            if doCutflows:
                makeCutflows(result_file, era, plot_dir, doPhotons)
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)
            VB.getValues(result_file, era)
            SB.getValues(result_file, era)
            makeDataCard(VB, dataCardValidation_dir, era)
            makeDataCard(SB, dataCardSearch_dir,     era)

    N.makeTexFile("validation", latex_dir + "validationBins_normalization_Zmass.tex")
    N.makeTexFile("search",     latex_dir + "searchBins_normalization_Zmass.tex")
    VB.makeTexFile("Z Invisible Prediction for Validation Bins", latex_dir + "zinv_prediction_validation_bins.tex")
    SB.makeTexFile("Z Invisible Prediction for Search Bins",     latex_dir + "zinv_prediction_search_bins.tex")
    
    # total Run2 prediction
    # root files to save histograms
    total_era = "Run2"
    validation_file = "validationBinsZinv_" + total_era + ".root"
    search_file     = "searchBinsZinv_"     + total_era + ".root"
    VB.makeTotalPred( validation_file,  "Validation Bin",   "validation", total_era   )
    SB.makeTotalPred( search_file,      "Search Bin",       "search",     total_era   )

    # TODO: making data card for Run 2 does not work because we have not run calcPrediction() for Run 2
    #       calcPrediction() depends on norm and shape (which we calculate per era, not for all of Run 2)
    #       make a way (possibly another function) to calculate values for Run 2 data card  
    #makeDataCard(VB, dataCardValidation_dir, total_era)
    #makeDataCard(SB, dataCardSearch_dir,     total_era)


if __name__ == "__main__":
    main()


