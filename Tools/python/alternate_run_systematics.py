# run_systematics.py

import copy
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
from tools import setupHist
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir):
    eraTag = "_" + era
    draw_option = "hist"
            
    # colors
    color_red    = "vermillion"
    color_blue   = "electric blue"
    color_green  = "irish green" 
    color_purple = "violet"
    color_black  = "black"
    
    # legend: TLegend(x1,y1,x2,y2)
    legend_x1 = 0.6
    legend_x2 = 0.9 
    legend_y1 = 0.7
    legend_y2 = 0.9 
    
    c = ROOT.TCanvas("c", "c", 800, 800)
    c.Divide(1, 2)
    
    name = "{0}_{1}_syst".format(bintype, mySyst)
    
    h_ratio_up      = h_up.Clone("h_ratio_up") 
    h_ratio_down    = h_down.Clone("h_ratio_down") 
    h_ratio_up.Divide(h)
    h_ratio_down.Divide(h)
    
    title = "Z to Invisible: " + name + " in " + region + " for " + era
    x_title = bintype + " bins"
    setupHist(h,                title, x_title, "Events",               color_black,  10.0 ** -2, 10.0 ** 5)
    setupHist(h_up,             title, x_title, "Events",               color_red,    10.0 ** -2, 10.0 ** 5)
    setupHist(h_down,           title, x_title, "Events",               color_blue,   10.0 ** -2, 10.0 ** 5)
    setupHist(h_ratio_up,       title, x_title, "variation / nominal",  color_red,    0.5, 1.5)
    setupHist(h_ratio_down,     title, x_title, "variation / nominal",  color_blue,   0.5, 1.5)
    
    # histograms
    c.cd(1)
    ROOT.gPad.SetLogy(1) # set log y
    # draw
    h.Draw(draw_option)
    h_up.Draw(draw_option + " same")
    h_down.Draw(draw_option + " same")
    
    # legend: TLegend(x1,y1,x2,y2)
    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
    legend.AddEntry(h,              "Z#rightarrow#nu#nu MC",  "l")
    legend.AddEntry(h_up,           "syst up",                "l")
    legend.AddEntry(h_down,         "syst down",              "l")
    legend.Draw()
    
    # ratios
    pad = c.cd(2)
    pad.SetGrid()
    # draw
    h_ratio_up.Draw(draw_option)
    h_ratio_down.Draw(draw_option + " same")
    
    # save histograms
    plot_name = plot_dir + name + "_" + region + eraTag
    c.Update()
    c.SaveAs(plot_name + ".png")



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
    validationHistMap         = {}

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
    #eras        = ["2018_PreHEM"]
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
            validationHistMap[era] = {}
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
                validationHistMap[era][region] = {}
                validationHistMap[era][region]["pred"] = VB.histograms[era][region]["pred"].Clone()
                
                # print for testing
                # nBins = h.GetNbinsX()
                # for i in xrange(1, nBins + 1):
                #     print "{0}, {1}, bin {2}: pred = {3}".format(era, region, i, h.GetBinContent(i))
                
                # syst up/down predictions
                #mySyst = "jes"
                #mySyst = "btag"
                #for mySyst in systMap: 
                # TODO: fix eff_restoptag_syst, isr_syst, pdf_syst
                systematics = ["jes", "btag", "pileup", "prefire"]
                for mySyst in systematics: 
                    validationHistMap[era][region]["syst_" + mySyst] = {}
                    # TODO: fix bug; final histogram plotted is wrong, e.g. last systematic, high dm, down
                    #for direction in ["down", "up"]:
                    for direction in ["up", "down"]:
                        # hist name format for reference
                        # std::string histSuffixSyst = "_" + syst + "_syst_" + w.first + JetPtCut;
                        # _name_syst_direction
                        systTag = "_{0}_syst_{1}".format(mySyst, direction)
                        print "region: {0}, systTag: {1}".format(region, systTag)
                        N.getNormAndError( result_file, era, systTag )
                        S.getShape(        result_file, era, systTag )
                        VB.getValues(      result_file, era, systTag )
                        SB.getValues(      result_file, era, systTag )
                        validationHistMap[era][region]["syst_" + mySyst][direction] = VB.histograms[era][region]["pred"].Clone()
                        
                        # print for testing
                        # nBins = h.GetNbinsX()
                        # for i in xrange(1, nBins + 1):
                        #     print "{0}, {1}, bin {2}: {3} = {4}".format(era, region, i, mySyst + "_" + direction, h.GetBinContent(i))
                    
                    # ---------------------- #
                    # --- Draw Histogram --- #
                    # ---------------------- #
                    h       = validationHistMap[era][region]["pred"]
                    h_up    = validationHistMap[era][region]["syst_" + mySyst]["up"]
                    h_down  = validationHistMap[era][region]["syst_" + mySyst]["down"]
                    plot(h, h_up, h_down, mySyst, "validation", region, era, plot_dir)


if __name__ == "__main__":
    main()


