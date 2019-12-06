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
    parser.add_argument("--runs_json",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--syst_json",    "-s", default="",                             help="json file containing systematics")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options         = parser.parse_args()
    runs_json       = options.runs_json
    syst_json       = options.syst_json
    verbose         = options.verbose
    doUnits         = False
    doPhotons       = False
    draw            = False
    saveRootFile    = False
    runMap          = {}
    systMap         = {}
    histMap         = {}

    if not os.path.exists(runs_json):
        print "The json file \"{0}\" containing runs does not exist.".format(runs_json)
        return
    if not os.path.exists(syst_json):
        print "The json file \"{0}\" containing systematics does not exist.".format(syst_json)
        return
    with open(runs_json, "r") as input_file:
        runMap  = json.load(input_file)
    with open(syst_json, "r") as input_file:
        systMap = json.load(input_file)
    
    #eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM", "Run2"]
    eras        = ["Run2"]
    dirList     = []
    plot_dir    = "syst_plots"
    dirList.append(plot_dir)
    # add "/" to directory if not present
    if plot_dir[-1]  != "/":  plot_dir += "/"
    
    for d in dirList:
        # make directory if it does not exist
        if not os.path.exists(d):
            os.makedirs(d)

    N  = Normalization(plot_dir, verbose)
    S  = Shape(plot_dir, draw, doUnits, verbose)
    VB = ValidationBins( N, S, eras, plot_dir, verbose, draw, saveRootFile )
    SB = SearchBins(     N, S, eras, plot_dir, verbose, draw, saveRootFile )
    
    with open(runs_json, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            histMap[era] = {}
            print "|---------- Era: {0} ----------|".format(era)
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)
            VB.getValues(result_file, era)
            SB.getValues(result_file, era)
            # get histograms
            # central prediction
            for region in VB.histograms[era]:
                histMap[era][region] = {}
                h = VB.histograms[era][region]["pred"]
                histMap[era][region]["pred"] = h
                # print for testing
                nBins = h.GetNbinsX()
                for i in xrange(1, nBins + 1):
                    print "{0}, {1}, bin {2}: {3}".format(era, region, i, h.GetBinContent(i))
            # syst up/down predictions


if __name__ == "__main__":
    main()




