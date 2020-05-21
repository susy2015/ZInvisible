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
    f = ROOT.TFile(result_file, "read")

    # FatJet_nGenPart histograms
    h_TTbar_MergedTop_nGenPart_nmteq1 = f.Get("FatJet_nGenPart/TTbar_MergedTop_nGenPart_recalc_nmteq1FatJet_nGenPartFatJet_nGenPart{FatJet_Stop0l=1}TTbar recalcsingle")
    h_ZNuNu_MergedTop_nGenPart_nmteq1 = f.Get("FatJet_nGenPart/ZNuNu_MergedTop_nGenPart_recalc_nmteq1FatJet_nGenPartFatJet_nGenPart{FatJet_Stop0l=1}ZJetsToNuNu recalcsingle")
    h_GJets_MergedTop_nGenPart_nmteq1 = f.Get("FatJet_nGenPart/GJets_MergedTop_nGenPart_recalc_nmteq1FatJet_nGenPartFatJet_nGenPart{FatJet_Stop0l=1}GJets recalcsingle")
    
    histograms = [h_TTbar_MergedTop_nGenPart_nmteq1, h_ZNuNu_MergedTop_nGenPart_nmteq1, h_GJets_MergedTop_nGenPart_nmteq1]
    labels  = ["TTbar", "ZNuNu", "GJets"]
    #histograms = [h_ZNuNu_MergedTop_nGenPart_nmteq1, h_GJets_MergedTop_nGenPart_nmteq1]
    #labels  = ["ZNuNu", "GJets"]
    name_1  = "MergedTop_nGenPart_nmteq1" 
    name_2  = "MergedTop_nGenPart_nmteq1_norm" 
    title_1 = "MergedTop_nGenPart, nmt=1, {0}".format(era)
    title_2 = "MergedTop_nGenPart, nmt=1, {0}, norm.".format(era)
    x_title = "MergedTop_nGenPart"
    x_min   = 0
    x_max   = 11
    y_min_1 = 0
    y_max_1 = 300
    y_min_2 = 0
    y_max_2 = 0.3
    # plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, era, showStats=False, normalize=False)
    plot(histograms, labels, name_1, title_1, x_title, x_min, x_max, y_min_1, y_max_1, era, showStats=True, normalize=False)
    plot(histograms, labels, name_2, title_2, x_title, x_min, x_max, y_min_2, y_max_2, era, showStats=True, normalize=True)


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
    
    eras = ["2016", "2017", "2018", "Run2"]
    #eras = ["Run2"]
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            run(era, result_file, verbose)

if __name__ == "__main__":
    main()

