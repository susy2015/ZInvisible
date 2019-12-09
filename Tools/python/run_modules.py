# run_modules.py

import argparse
import json
import os
from cutflow_plot import makeCutflows 
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from systematics import Systematic
from search_bins import SearchBins, ValidationBins, CRUnitBins, SRUnitBins
from data_card import makeDataCard
from units import saveResults
from make_table import Table

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options     = parser.parse_args()
    json_file   = options.json_file
    verbose     = options.verbose

    doUnits     = True
    doCutflows  = False
    doPhotons   = False
    draw        = False
    

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM", "Run2"]
    #eras = ["2016"]
    dirList = []
    plot_dir                = "more_plots"
    latex_dir               = "latex_files"
    results_dir             = "results"
    dataCard_dir            = "data_cards"
    dataCardValidation_dir  = dataCard_dir + "/validation"
    dataCardSearch_dir      = dataCard_dir + "/search"
    dirList.append(plot_dir)
    dirList.append(latex_dir)
    dirList.append(results_dir)
    dirList.append(dataCard_dir)
    dirList.append(dataCardValidation_dir)
    dirList.append(dataCardSearch_dir)
    # add "/" to directory if not present
    if plot_dir[-1]  != "/":                plot_dir                += "/"
    if latex_dir[-1] != "/":                latex_dir               += "/"
    if results_dir[-1] != "/":              results_dir             += "/"
    if dataCard_dir[-1] != "/":             dataCard_dir            += "/"
    if dataCardValidation_dir[-1] != "/":   dataCardValidation_dir  += "/"
    if dataCardSearch_dir[-1] != "/":       dataCardSearch_dir      += "/"
    
    for d in dirList:
        # make directory if it does not exist
        if not os.path.exists(d):
            os.makedirs(d)
    # normalization
    N = Normalization(plot_dir, verbose)
    # shape
    S = Shape(plot_dir, draw, doUnits, verbose)
    # validation bins
    VB = ValidationBins(    N, S, eras, plot_dir, verbose, draw=True, saveRootFile=True )
    # search bins
    SB = SearchBins(        N, S, eras, plot_dir, verbose, draw=True, saveRootFile=True )
    if doUnits:
        # control region unit bins  
        CRunits = CRUnitBins(N, S, eras, plot_dir, verbose) 
        # search region unit bins  
        SRunits = SRUnitBins(N, S, eras, plot_dir, verbose) 
    # systematics
    Syst = Systematic(plot_dir, N, S)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
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
            if doUnits:
                CRunits.getValues(result_file, era)
                SRunits.getValues(result_file, era)
            Syst.makeZvsPhoton(result_file, era, False)
            Syst.makeZvsPhoton(result_file, era, True)
            makeDataCard(VB, dataCardValidation_dir, era)
            makeDataCard(SB, dataCardSearch_dir,     era)

    N.makeTexFile("validation", latex_dir + "validationBins_normalization_Zmass.tex")
    N.makeTexFile("search",     latex_dir + "searchBins_normalization_Zmass.tex")
    N.makeComparison("validation")
    N.makeComparison("search")
    S.makeComparison("validation")
    S.makeComparison("search")
    VB.makeTexFile("Z Invisible Per Era Prediction for Validation Bins", latex_dir + "zinv_per_era_prediction_validation_bins.tex")
    SB.makeTexFile("Z Invisible Per Era Prediction for Search Bins",     latex_dir + "zinv_per_era_prediction_search_bins.tex")
    
    # total Run2 prediction
    # root files to save histograms
    total_era = "Run2"
    validation_file = "validationBinsZinv_" + total_era + ".root"
    search_file     = "searchBinsZinv_"     + total_era + ".root"
    # WARNING: only run makeTotalPred() if you do not already have Run2 combined histograms; otherwise you will double count!
    #VB.makeTotalPred( validation_file,  "Validation Bin",   "validation", total_era   )
    #SB.makeTotalPred( search_file,      "Search Bin",       "search",     total_era   )
    VB.makeTexFile("Z Invisible Total Prediction for Validation Bins", latex_dir + "zinv_total_prediction_validation_bins.tex", total_era)
    SB.makeTexFile("Z Invisible Total Prediction for Search Bins",     latex_dir + "zinv_total_prediction_search_bins.tex",     total_era)

    # make json files
    VB.makeJson(VB.binValues,           results_dir + "ValidationBinResults.json")
    SB.makeJson(SB.binValues,           results_dir + "SearchBinResults.json")
    if doUnits:
        CRunits.makeJson(CRunits.binValues, results_dir + "CRUnitsResults.json")
        SRunits.makeJson(SRunits.binValues, results_dir + "SRUnitsResults.json")
        # saveResults(inFile, outFile, CRunits, SRunits, era)
        saveResults("dc_BkgPred_BinMaps_master.json", results_dir + "zinv_yields_" + total_era + ".json", CRunits, SRunits, total_era)

    # TODO: making data card for Run 2 does not work because we have not run calcPrediction() for Run 2
    #       calcPrediction() depends on norm and shape (which we calculate per era, not for all of Run 2)
    #       make a way (possibly another function) to calculate values for Run 2 data card  
    #makeDataCard(VB, dataCardValidation_dir, total_era)
    #makeDataCard(SB, dataCardSearch_dir,     total_era)


if __name__ == "__main__":
    main()


