# run_modules.py

import numpy as np
import json
import argparse
import os
import ROOT
from Norm_lepton_zmass import Normalization
from Shape_photon_met import Shape
from Search_bins import  ValidationBins, SearchBins
from tools import invert

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


#TODO: combined all systematics to get total syst up/down
# first combined syst in validation bins and give to Hui
# loop over validation bins
# loop over systematics
# syst_up_i       = (p_up - p)   / p
# syst_down_i     = (p - p_down) / p
# syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
# syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 

def calcTotalSyst(h_pred, systHistoMap):
    pass

def writeToConf(outFile, searchBinMap, syst, h, h_up, h_down, offset):
    zinv = "znunu"
    # sb_i = bin_i - 1 + offset
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = searchBinMap[str(i)]
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        #s_up   = (p_up - p)   / p
        #s_down = (p - p_down) / p
        # If there is a zero prediction for nominal, up, and down can you set the systematic to 1.
        # The error is a deviation from 1.
        r_up   = 1
        r_down = 1
        if p != 0:
            r_up   = p_up   / p
            r_down = p_down / p
        else:
            print "WARNING: pred = 0 for search bin {0}".format(sb_i)
        outFile.write("{0}  {1}_Up  {2}  {3}\n".format(   sb_name, syst, zinv, r_up))
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, syst, zinv, r_down))


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
    units_json      = "dc_BkgPred_BinMaps_master.json"
   
    doUnits      = False
    draw         = False
    runMap       = {}
    systMap      = {}
    unitMap      = {}

    # check if files exist
    if not os.path.exists(runs_json):
        print "The json file \"{0}\" containing runs does not exist.".format(runs_json)
        return
    if not os.path.exists(syst_json):
        print "The json file \"{0}\" containing systematics does not exist.".format(syst_json)
        return
    if not os.path.exists(units_json):
        print "The json file \"{0}\" containing systematics does not exist.".format(units_json)
        return
    
    # load json files
    with open(runs_json, "r") as input_file:
        runMap  = json.load(input_file)
    with open(syst_json, "r") as input_file:
        systMap = json.load(input_file) 
    with open(units_json, "r") as input_file:
        unitMap = json.load(input_file) 
    
    # map search bin numbers to string names
    searchBinMap = invert(unitMap["binNum"])

    era          = "Run2"
    runDir = runMap[era]
    result_file  = "condor/" + runDir + "/result.root"
    conf_file    = "zinv_syst_" + era + ".conf"
    out_dir      = "angel_syst_plots/" 
    syst_dir     = "angel_syst_plots/"
    variable     = "pred"
    regions      = ["lowdm", "highdm"]
    directions   = ["up", "", "down"]
    bintypes     = ["validation", "search"]
    systematics  = ["btag","eff_restoptag","eff_sb","eff_toptag","eff_wtag","met_trig","pileup"] # all systematics available from Znunu_nValidationBin
    
    histo_tmp  = {region:dict.fromkeys(directions) for region in regions} # histo_tmp[region][direction]
    histo      = {bintype:histo_tmp for bintype in bintypes}              # histo[bintype][region][direction]
    syst_histo = {syst:histo for syst in systematics}                     # syst_histo[systemaitc][bintype][region][direction]
    
    #-------------------------------------------------------
    # Class instanceses summoning 
    #-------------------------------------------------------
    
    N = Normalization(out_dir, verbose)
    S = Shape(out_dir, draw, doUnits, verbose)
    VB = ValidationBins(N, S, era, out_dir, verbose)
    SB = SearchBins(N, S, era, out_dir, verbose)
    
    #-------------------------------------------------------
    # Normal predictions (no systematics) 
    #-------------------------------------------------------
    
    N.getNormAndError(result_file, "", era)
    S.getShape(result_file, "", era)
    VB.getValues(result_file, "", era)
    SB.getValues(result_file, "", era)
    
    histo["validation"]["lowdm"][""]   =  VB.histograms[era]["lowdm"][variable].Clone()
    histo["validation"]["highdm"][""]  =  VB.histograms[era]["highdm"][variable].Clone()
    histo["search"]["lowdm"][""]       =  SB.histograms[era]["lowdm"][variable].Clone()
    histo["search"]["highdm"][""]      =  SB.histograms[era]["highdm"][variable].Clone()
    
    #-------------------------------------------------------
    # Calculate normalization and shape factors
    #-------------------------------------------------------
    with open(conf_file, "w") as outFile:
           
        for syst in systematics:
        
            # --- syst up --- #
            N.getNormAndError(result_file, syst + "_syst_up", era)
            S.getShape(result_file, syst + "_syst_up", era)
            VB.getValues(result_file, syst + "_syst_up", era)
            SB.getValues(result_file, syst + "_syst_up", era)
            
            histo["validation"]["highdm"]["up"]  = VB.histograms[era]["highdm"][variable].Clone()
            histo["validation"]["lowdm"]["up"]   = VB.histograms[era]["lowdm"][variable].Clone()
            histo["search"]["highdm"]["up"]      = SB.histograms[era]["highdm"][variable].Clone()
            histo["search"]["lowdm"]["up"]       = SB.histograms[era]["lowdm"][variable].Clone()
            syst_histo[syst]["validation"]["lowdm"]["up"]  = histo["validation"]["lowdm"]["up"]
            syst_histo[syst]["validation"]["highdm"]["up"] = histo["validation"]["highdm"]["up"]
            syst_histo[syst]["search"]["lowdm"]["up"]      = histo["search"]["lowdm"]["up"]
            syst_histo[syst]["search"]["highdm"]["up"]     = histo["search"]["highdm"]["up"]
            
            # --- syst down --- #
            N.getNormAndError(result_file, syst + "_syst_down", era)
            S.getShape(result_file, syst + "_syst_down", era)
            VB.getValues(result_file, syst + "_syst_down", era)
            SB.getValues(result_file, syst + "_syst_down", era)
            
            histo["validation"]["highdm"]["down"]   = VB.histograms[era]["highdm"][variable].Clone()
            histo["validation"]["lowdm"]["down"]    = VB.histograms[era]["lowdm"][variable].Clone()
            histo["search"]["highdm"]["down"]       = SB.histograms[era]["highdm"][variable].Clone()
            histo["search"]["lowdm"]["down"]        = SB.histograms[era]["lowdm"][variable].Clone()
            syst_histo[syst]["validation"]["lowdm"]["down"]  = histo["validation"]["lowdm"]["down"]
            syst_histo[syst]["validation"]["highdm"]["down"] = histo["validation"]["highdm"]["down"]
            syst_histo[syst]["search"]["lowdm"]["down"]      = histo["search"]["lowdm"]["down"]
            syst_histo[syst]["search"]["highdm"]["down"]     = histo["search"]["highdm"]["down"]
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------
            #writeToConf(outFile, searchBinMap, syst, h, h_up, h_down, offset)
            systForConf = systMap[syst]["name"]  
            writeToConf(outFile, searchBinMap, systForConf, histo["search"]["lowdm"][""],  histo["search"]["lowdm"]["up"],  histo["search"]["lowdm"]["down"],  0)
            writeToConf(outFile, searchBinMap, systForConf, histo["search"]["highdm"][""], histo["search"]["highdm"]["up"], histo["search"]["highdm"]["down"], SB.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------
            
            for bintype in bintypes:
                for region in regions:
                
                    # legend: TLegend(x1,y1,x2,y2)
                    legend_x1 = 0.7 
                    legend_x2 = 0.9 
                    legend_y1 = 0.7 
                    legend_y2 = 0.9 
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                
                    c = ROOT.TCanvas("c", "c", 800, 800)
                    c.Divide(1, 2)
                
                    c.cd(1)
                    
                    histo[bintype][region][""].SetLineColor(ROOT.kBlue)
                    histo[bintype][region][""].SetTitle(bintype + " bins, " + syst + " systematic, " + era)
                    histo[bintype][region][""].Draw("hist")
                    histo[bintype][region]["up"].SetLineColor(ROOT.kRed)
                    histo[bintype][region]["up"].Draw("same")
                    histo[bintype][region]["down"].SetLineColor(ROOT.kViolet)
                    histo[bintype][region]["down"].Draw("same")
                    
                    ROOT.gPad.SetLogy(1) # set log y
                    
                    legend.AddEntry(histo[bintype][region]["up"]   , "up"      , "l")
                    legend.AddEntry(histo[bintype][region][""]     , "nominal" , "l")
                    legend.AddEntry(histo[bintype][region]["down"] , "down"    , "l")
                    legend.Draw()
                    
                    c.cd(2)
                
                    ratio_up = histo[bintype][region]["up"].Clone()
                    ratio_up.Divide(histo[bintype][region][""])
                    ratio_up.SetLineColor(ROOT.kRed)
                
                    ratio_down = histo[bintype][region]["down"].Clone()
                    ratio_down.Divide(histo[bintype][region][""])
                    ratio_down.SetLineColor(ROOT.kViolet)
                
                    ratio_up.GetYaxis().SetRangeUser(0,1.5)
                
                    ratio_up.Draw("hist")
                    ratio_down.Draw("hist same")
                    
                    c.Update()
                    
                    file_name = syst_dir + bintype + "_syst_" + syst + "_" + region
                    c.SaveAs(file_name + ".png")

                    del c
    
    #-------------------------------------------------------
    # Calculate total systematic up/down
    #-------------------------------------------------------

    #TODO: combined all systematics to get total syst up/down
    # first combined syst in validation bins and give to Hui
    # loop over validation bins
    # loop over systematics
    # syst_up_i       = (p_up - p)   / p
    # syst_down_i     = (p - p_down) / p
    # syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
    # syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 
    
    validationBinMap = {}
    validationBinMap["lowdm"]   = VB.low_dm_bins
    validationBinMap["highdm"]  = VB.high_dm_bins
    
    h_pred_lowdm    = histo["validation"]["lowdm"][""]
    h_pred_highdm   = histo["validation"]["highdm"][""]

    # validation bins
    # bins are list of strings starting at 0
    # loop over regions (lowdm and highdm)
    print "# validation bin systematics"
    for region in regions:
        print "--- {0} ---".format(region)
        for b in validationBinMap[region]:
            b_i = int(b) + 1
            p = histo["validation"][region][""].GetBinContent(b_i)
            syst_up_sum   = 0.0
            syst_down_sum = 0.0
            if p != 0:
                for syst in systematics:
                    # syst_histo[systemaitc][bintype][region][direction]
                    h_up    = syst_histo[syst]["validation"][region]["up"]
                    h_down  = syst_histo[syst]["validation"][region]["down"]
                    p_up    = h_up.GetBinContent(b_i)
                    p_down  = h_down.GetBinContent(b_i)
                    syst_up   = (p_up - p  ) / p
                    syst_down = (p - p_down) / p
                    syst_up_sum     += syst_up**2
                    syst_down_sum   += syst_down**2
            syst_up_total   = np.sqrt(syst_up_sum)
            syst_down_total = np.sqrt(syst_down_sum)
            final_up   = 1.0 + syst_up_total
            final_down = 1.0 - syst_down_total
            print "bin {0}, pred={1}, syst_up={2}, syst_down={3}".format(b_i, p, final_up, final_down)



if __name__ == "__main__":
    main()



