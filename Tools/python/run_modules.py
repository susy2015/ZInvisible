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

    validation = True
    doCutflows = False
    doPhotons = True
    useNbNsvSelection = True
    draw = False

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM"]
    dirList = []
    plot_dir        = "more_plots"
    latex_dir       = "latex_files"
    dataCard_dir    = "data_cards"
    dirList.append(plot_dir)
    dirList.append(latex_dir)
    dirList.append(dataCard_dir)
    # add "/" to directory if not present
    if plot_dir[-1]  != "/":    plot_dir     += "/"
    if latex_dir[-1] != "/":    latex_dir    += "/"
    if dataCard_dir[-1] != "/": dataCard_dir += "/"
    
    for d in dirList:
        # make directory if it does not exist
        if not os.path.exists(d):
            os.makedirs(d)

    N = Normalization(validation, useNbNsvSelection, verbose)
    S = Shape(plot_dir, draw, verbose)

    table_file = open("njets_table.txt", "w+") 
    T = Table(table_file)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # search bins
        SB = SearchBins(N, S, eras, plot_dir, verbose)
        # validation bins
        VB = ValidationBins(N, S, eras, plot_dir, verbose)
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
            makeDataCard(VB, dataCard_dir, era)
            T.makeTable(result_file, era)

    N.makeTexFile(latex_dir + "normalization_Zmass.tex")
    VB.makeTexFile(latex_dir + "zinv_prediction.tex")

    table_file.close()

if __name__ == "__main__":
    main()



