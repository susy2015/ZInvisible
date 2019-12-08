# run_modules.py

import json
import argparse
import os
from Norm_lepton_zmass import Normalization
from Shape_photon_met import Shape
from Search_bins import  ValidationBins
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

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
   
    doUnits      = False
    draw         = False
    systMap      = {}
    histMap      = {}

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
    
    era          = "Run2"
    runDir = runMap[era]
    result_file = "condor/" + runDir + "/result.root"
    out_dir      = "validation/" 
    syst_dir     = "angel_syst_plots/"
    variable     = "pred"
    regions      = ["lowdm", "highdm"]
    directions   = ['up', '', 'down']
    systematics  = ['btag','eff_restoptag','eff_sb','eff_toptag','eff_wtag','met_trig','pileup'] # all systematics available from Znunu_nValidationBin
    
    histo = {region:dict.fromkeys(directions) for region in regions} # histo[regioin][direction]
    
    #-------------------------------------------------------
    # Class instanceses summoning 
    #-------------------------------------------------------
    
    N = Normalization(out_dir, verbose)
    S = Shape(out_dir, draw, doUnits, verbose)
    VB = ValidationBins(N, S, era, out_dir, verbose)
    
    #-------------------------------------------------------
    # Normal predictions (no systematics) 
    #-------------------------------------------------------
    
    N.getNormAndError(result_file, "", era)
    S.getShape(result_file, "", era)
    VB.getValues(result_file, "", era)
    
    histo["lowdm"][""]   =  VB.histograms[era]["lowdm"][variable].Clone()
    histo["highdm"][""]  =  VB.histograms[era]["highdm"][variable].Clone()
    
    #-------------------------------------------------------
    # Calculate normalization and shape factors
    #-------------------------------------------------------
    
    for syst in systematics:
    
        N.getNormAndError(result_file, syst + "_syst_up", era)
        S.getShape(result_file, syst + "_syst_up", era)
        print("Everything is fine so far")
        VB.getValues(result_file, syst + "_syst_up", era)
        
        histo["highdm"]["up"]  =  VB.histograms[era]["highdm"][variable].Clone()
        histo["lowdm"]["up"]   =  VB.histograms[era]["lowdm"][variable].Clone()
        
        N.getNormAndError(result_file, syst + "_syst_down", era)
        S.getShape(result_file, syst + "_syst_down", era)
        VB.getValues(result_file, syst + "_syst_down", era)
        
        histo["highdm"]["down"]  =  VB.histograms[era]["highdm"][variable].Clone()
        histo["lowdm"]["down"]   =  VB.histograms[era]["lowdm"][variable].Clone()
        
        #-------------------------------------------------------
        # Canvas 
        #-------------------------------------------------------
        
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
            
            histo[region][""].Draw("hist")
            histo[region][""].SetLineColor(ROOT.kBlue)
            
            histo[region]["up"].Draw("same")
            histo[region]["up"].SetLineColor(ROOT.kRed)
            
            histo[region]["down"].Draw("same")
            histo[region]["down"].SetLineColor(ROOT.kViolet)
            
            ROOT.gPad.SetLogy(1) # set log y
            
            legend.AddEntry(histo[region]["up"]   , "up"      , "l")
            legend.AddEntry(histo[region][""]     , "nominal" , "l")
            legend.AddEntry(histo[region]["down"] , "down"    , "l")
            legend.Draw()
            
            c.cd(2)
        
            ratio_up = histo[region]["up"].Clone()
            ratio_up.Divide(histo[region][""])
            ratio_up.SetLineColor(ROOT.kRed)
        
            ratio_down = histo[region]["down"].Clone()
            ratio_down.Divide(histo[region][""])
            ratio_down.SetLineColor(ROOT.kViolet)
        
            ratio_up.GetYaxis().SetRangeUser(0,1.5)
        
            ratio_up.Draw("hist")
            ratio_down.Draw("hist same")
            
            c.Update()
            
            file_name = syst_dir + "Syst_" + syst + "_" + region
            c.SaveAs(file_name + ".png")


if __name__ == "__main__":
    main()



