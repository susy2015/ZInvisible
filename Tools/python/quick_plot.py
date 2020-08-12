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
    # plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, era, plot_dir, showStats=False, normalize=False, setLog=False)
    plot(histograms, labels, name_1, title_1, x_title, x_min, x_max, y_min_1, y_max_1, era, "more_plots", showStats=False, normalize=False, setLog=True)
    plot(histograms, labels, name_2, title_2, x_title, x_min, x_max, y_min_2, y_max_2, era, "more_plots", showStats=True, normalize=True)

# plot Nb, Nsv, Nj for different eras
def plotVars(var, particle, eras, runMap, varMap, verbose):
    var_label = varMap[var]["label"] 
    h_var     = varMap[var]["h_var"]  
    histogramsData  = []
    histogramsMC    = []
    histogramsRatio = []
    labelsData      = []
    labelsMC        = []
    labelsRatio     = []
    for era in eras:
        runDir = runMap[era]
        result_file = "condor/" + runDir + "/result.root"
        if verbose:
            print "{0}: {1}".format(era, result_file)
        # read input file
        f = ROOT.TFile(result_file, "read")
        # --- Data --- #
        h_name = "{0}/DataMC_{1}_Baseline_{2}_jetpt30{0}{0}Datadata".format(h_var, particle, var)
        h_data = f.Get(h_name)
        # use deepcopy so that histogram exists after file is closed / reassigned
        histogramsData.append(copy.deepcopy(h_data))
        labelsData.append(era + " Data")
        if not histogramsData[-1]:
            print "ERROR: Unable to load histogram {0}".format(h_name) 
        # --- MC --- #
        h1_name = "{0}/DataMC_{1}_Baseline_{2}_jetpt30{0}{0}DYstack".format(h_var, particle, var)
        h2_name = "{0}/DataMC_{1}_Baseline_{2}_jetpt30{0}{0}t#bar{{t}}stack".format(h_var, particle, var)
        h3_name = "{0}/DataMC_{1}_Baseline_{2}_jetpt30{0}{0}Single tstack".format(h_var, particle, var)
        h4_name = "{0}/DataMC_{1}_Baseline_{2}_jetpt30{0}{0}Rarestack".format(h_var, particle, var)
        h1 = f.Get(h1_name)
        h2 = f.Get(h2_name)
        h3 = f.Get(h3_name)
        h4 = f.Get(h4_name)
        h_mc = h1.Clone("h_mc")
        h_mc.Add(h2)
        h_mc.Add(h3)
        h_mc.Add(h4)
        # use deepcopy so that histogram exists after file is closed / reassigned
        histogramsMC.append(copy.deepcopy(h_mc))
        labelsMC.append(era + " MC")
        if not histogramsMC[-1]:
            print "ERROR: Unable to load histogram {0}".format(h_name) 
        # --- Ratio --- #
        h_ratio = h_data.Clone("h_ratio")
        h_ratio.Divide(h_mc)
        # use deepcopy so that histogram exists after file is closed / reassigned
        histogramsRatio.append(copy.deepcopy(h_ratio))
        labelsRatio.append(era + " Data/MC")
        if not histogramsRatio[-1]:
            print "ERROR: Unable to load histogram {0}".format(h_name) 
    nameData    = "{0}_data_{1}".format(particle, var_label)
    titleData   = "{0} data comparison for {1}".format(particle, var_label)
    nameMC      = "{0}_mc_{1}".format(particle, var_label)
    titleMC     = "{0} MC comparison for {1}".format(particle, var_label)
    nameRatio   = "{0}_ratio_{1}".format(particle, var_label)
    titleRatio  = "{0} Data/MC comparison for {1}".format(particle, var_label)
    x_title     = var_label
    x_min = 0
    x_max = 11
    # standard y-axis limits
    #y_min_1 = 10**-2
    #y_max_1 = 10**6
    # normalized y-axis limits
    y_min_1 = 10**-7
    y_max_1 = 10**1
    y_min_2 = 0.0
    y_max_2 = 3.0
    # plot(histograms, labels, name, title, x_title, x_min, x_max, y_min, y_max, era, showStats=False, normalize=False, setLog=False)
    plot(histogramsData, labelsData, nameData, titleData, x_title, x_min, x_max, y_min_1, y_max_1, "Run2", showStats=False, normalize=True, setLog=True)
    plot(histogramsMC, labelsMC, nameMC, titleMC, x_title, x_min, x_max, y_min_1, y_max_1, "Run2", showStats=False, normalize=True, setLog=True)
    plot(histogramsRatio, labelsRatio, nameRatio, titleRatio, x_title, x_min, x_max, y_min_2, y_max_2, "Run2", showStats=False, normalize=False, setLog=False)

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
    varMap = {}
    varMap["nb"] = {}
    varMap["nb"]["label"] = "Nb"
    varMap["nb"]["h_var"] = "nBottoms_drLeptonCleaned_jetpt30"
    varMap["nsv"] = {}
    varMap["nsv"]["label"] = "Nsv"
    varMap["nsv"]["h_var"] = "nSoftBottoms_drLeptonCleaned_jetpt30"
    varMap["nj"] = {}
    varMap["nj"]["label"] = "Nj"
    varMap["nj"]["h_var"] = "nJets_drLeptonCleaned_jetpt30"
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            plotMergedTopPartons(era, result_file, verbose)
        plotVars("nb",  "Electron", eras, runMap, varMap, verbose)
        plotVars("nb",  "Muon",     eras, runMap, varMap, verbose)
        plotVars("nsv", "Electron", eras, runMap, varMap, verbose)
        plotVars("nsv", "Muon",     eras, runMap, varMap, verbose)
        plotVars("nj",  "Electron", eras, runMap, varMap, verbose)
        plotVars("nj",  "Muon",     eras, runMap, varMap, verbose)

if __name__ == "__main__":
    main()

