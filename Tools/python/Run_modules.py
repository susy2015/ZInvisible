# run_modules.py

import os
from Norm_lepton_zmass import Normalization
from Shape_photon_met import Shape
from Search_bins import  ValidationBins
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
    
era     = "2018_PreHEM"
verbose  = False 
doUnits  = False
draw     = False

out_dir     = "validation/" 
result_file = "result.root"
regions = ["lowdm", "highdm"]
directions = ['up', '', 'down']

histo = {region:dict.fromkeys(directions) for region in regions} # histo[regioin][direction]

#-------------------------------------------------------
# Class instanceses construction
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

histo["lowdm"][""]   =  VB.histograms[era]["lowdm"]["mc"].Clone()
histo["highdm"][""]  =  VB.histograms[era]["highdm"]["mc"].Clone()

#-------------------------------------------------------
# Calculate normalization and shape factors
#-------------------------------------------------------

#syst = "btag"
syst = "eff_sb_photon"

N.getNormAndError(result_file, syst + "_syst_up", era)
S.getShape(result_file, syst + "_syst_up", era)
VB.getValues(result_file, syst + "_syst_up", era)

histo["highdm"]["up"]  =  VB.histograms[era]["highdm"]["mc"].Clone()
histo["lowdm"]["up"]   =  VB.histograms[era]["lowdm"]["mc"].Clone()

N.getNormAndError(result_file, syst + "_syst_down", era)
S.getShape(result_file, syst + "_syst_down", era)
VB.getValues(result_file, syst + "_syst_down", era)

histo["highdm"]["down"]  =  VB.histograms[era]["highdm"]["mc"].Clone()
histo["lowdm"]["down"]   =  VB.histograms[era]["lowdm"]["mc"].Clone()

#-------------------------------------------------------
# Canvas 
#-------------------------------------------------------

for region in regions:

    c = ROOT.TCanvas("c", "c", 800, 800)
    
    histo[region][""].Draw("hist")
    histo[region][""].SetLineColor(ROOT.kBlue)
    
    histo[region]["up"].Draw("same")
    histo[region]["up"].SetLineColor(ROOT.kRed)
    
    histo[region]["down"].Draw("same")
    histo[region]["down"].SetLineColor(ROOT.kViolet)
    
    ROOT.gPad.SetLogy(1) # set log y
    
    # legend: TLegend(x1,y1,x2,y2)
    legend_x1 = 0.7 
    legend_x2 = 0.9 
    legend_y1 = 0.7 
    legend_y2 = 0.9 
    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
    legend.AddEntry(histo[region]["up"]   , "up"      , "l")
    legend.AddEntry(histo[region][""]     , "nominal" , "l")
    legend.AddEntry(histo[region]["down"] , "down"    , "l")
    
    legend.Draw()
    
    c.Update()
    
    file_name = "Syst_" + syst + "_" + region
    c.SaveAs(file_name + ".png")

