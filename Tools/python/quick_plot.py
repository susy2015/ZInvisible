# quick_plot.py

import ROOT
import copy
import json
import os
import argparse
import numpy as np
from tools import plot

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# plot merged top partons
def plotMergedTopPartons(era, result_file, verbose):
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
    # y limits for TTbar, ZNuNu, GJets
    y_min_1 = 10**-1
    y_max_1 = 10**5
    y_min_2 = 0
    y_max_2 = 0.8
    # y limits for ZNuNu, GJets
    #y_min_1 = 0
    #y_max_1 = 300
    #y_min_2 = 0
    #y_max_2 = 0.3
    
    # WARNING: if using setLog=True, do not use y_min = 0
    # WARNING: currently stats do not show properly on log scale... this would need work to fix
    # plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, era, showStats=False, normalize=False, setLog=False)
    plot(histograms, labels, name_1, title_1, x_title, x_min, x_max, y_min_1, y_max_1, era, showStats=False, normalize=False, setLog=True)
    plot(histograms, labels, name_2, title_2, x_title, x_min, x_max, y_min_2, y_max_2, era, showStats=True, normalize=True)

# plot Nb, Nsv, Nj for different eras
def plotVars(eras, runMap, verbose):
    var = "Nb"
    histograms = []
    labels     = []
    for era in eras:
        runDir = runMap[era]
        result_file = "condor/" + runDir + "/result.root"
        if verbose:
            print "{0}: {1}".format(era, result_file)
        # use deepcopy so that histogram exists after file is closed / reassigned
        f = ROOT.TFile(result_file, "read")
        h_name = "nBottoms_drLeptonCleaned_jetpt30/DataMC_Electron_Baseline_nb_jetpt30nBottoms_drLeptonCleaned_jetpt30nBottoms_drLeptonCleaned_jetpt30Datadata"
        histograms.append(copy.deepcopy(f.Get(h_name)))
        labels.append(era + " data")
        if not histograms[-1]:
            print "ERROR: Unable to load histogram {0}".format(h_name) 
    name = var + "_data"
    title = "Data comparison for " + var
    x_title = var
    x_min = 0
    x_max = 6
    #y_min = 10**-1
    #y_max = 10**5
    y_min = 10**-2
    y_max = 10**1
    # plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, era, showStats=False, normalize=False, setLog=False)
    plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, "Run2", showStats=False, normalize=True, setLog=True)

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
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            plotMergedTopPartons(era, result_file, verbose)
        plotVars(eras, runMap, verbose)

if __name__ == "__main__":
    main()

