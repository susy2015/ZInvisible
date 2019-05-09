# cutflowPlot.py

import os
import ROOT
from colors import getColorIndex

def setupHist(hist, labels, title, color, y_min, y_max):
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    # label bins
    for i, label in enumerate(labels, 1):
        #print i, label
        x_axis.SetBinLabel(i, label)
    
    y_axis.SetRangeUser(y_min, y_max)
    hist.SetTitle(title)
    hist.SetStats(ROOT.kFALSE)
    hist.SetLineColor(getColorIndex(color))
    hist.SetLineWidth(3)

def makePlots(f_name, year):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    f_name_Electron = f_name
    f_name_Muon     = f_name
    plot_dir = "cutflows_"+year+"/"
    h_dir = "CutFlows/"
    # check that files exist
    f_names = [f_name_Electron, f_name_Muon]
    for f_name in f_names:
        exists = os.path.isfile(f_name)
        if not exists: 
            print "The file {0} does not exist".format(f_name)
            return
    # make directory for plots if it does not exist
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    # colors
    red_color    = "vermillion"
    blue_color   = "electric blue"
    green_color  = "irish green" 
    purple_color = "violet"

    # y axis limits
    y_min = 10.0**1
    y_max = 10.0**10

    c = ROOT.TCanvas("c", "c", 800, 800)
    f_Electron = ROOT.TFile(f_name_Electron)
    f_Muon     = ROOT.TFile(f_name_Muon)
    # Data: apply trigger; MC: no trigger
    cuts = [
        "No cuts",
        "Trigger",
        "Filter",
        "JetID",
        "MET",
        "HT",
        "nJets",
        "Baseline",
    ]

    # Electron
    cutsElectronLowDM = cuts[:]
    cutsElectronLowDM.append("Low #Deltam")
    cutsElectronLowDM.append("Muon veto")
    cutsElectronLowDM.append("Z to ee")
    cutsElectronHighDM = cuts[:]
    cutsElectronHighDM.append("High #Deltam")
    cutsElectronHighDM.append("Muon veto")
    cutsElectronHighDM.append("Z to ee")
    # Muon 
    cutsMuonLowDM = cuts[:]
    cutsMuonLowDM.append("Low #Deltam")
    cutsMuonLowDM.append("Electron veto")
    cutsMuonLowDM.append("Z to #mu#mu")
    cutsMuonHighDM = cuts[:]
    cutsMuonHighDM.append("High #Deltam")
    cutsMuonHighDM.append("Electron veto")
    cutsMuonHighDM.append("Z to #mu#mu")
    # Lepton
    cutsLeptonLowDM = cuts[:]
    cutsLeptonLowDM.append("Low #Deltam")
    cutsLeptonLowDM.append("Lepton veto")
    cutsLeptonLowDM.append("Z to LL")
    cutsLeptonHighDM = cuts[:]
    cutsLeptonHighDM.append("High #Deltam")
    cutsLeptonHighDM.append("Lepton veto")
    cutsLeptonHighDM.append("Z to LL")
    
    h_map = {
        "CutFlow_Electron_LowDM"  : {"data" : "CutFlow_Data_Electron_LowDM_met",  "mc" : "CutFlow_MC_Electron_LowDM_met",  "cuts" : cutsElectronLowDM,  "file" : f_Electron},
        "CutFlow_Electron_HighDM" : {"data" : "CutFlow_Data_Electron_HighDM_met", "mc" : "CutFlow_MC_Electron_HighDM_met", "cuts" : cutsElectronHighDM, "file" : f_Electron},
        "CutFlow_Muon_LowDM"  : {"data" : "CutFlow_Data_Muon_LowDM_met",  "mc" : "CutFlow_MC_Muon_LowDM_met",  "cuts" : cutsMuonLowDM,  "file" : f_Muon},
        "CutFlow_Muon_HighDM" : {"data" : "CutFlow_Data_Muon_HighDM_met", "mc" : "CutFlow_MC_Muon_HighDM_met", "cuts" : cutsMuonHighDM, "file" : f_Muon}
    }
    plot_map = {
        "CutFlow_LowDM"  : {"Electron" : "CutFlow_Electron_LowDM",  "Muon" : "CutFlow_Muon_LowDM",  "cuts" :  cutsLeptonLowDM},
        "CutFlow_HighDM" : {"Electron" : "CutFlow_Electron_HighDM", "Muon" : "CutFlow_Muon_HighDM", "cuts" :  cutsLeptonHighDM}
    }
    
    for key in h_map:
        cutList = h_map[key]["cuts"]    
        f       = h_map[key]["file"]
        h_mc    = f.Get(h_dir + h_map[key]["mc"])
        h_data  = f.Get(h_dir + h_map[key]["data"])
        
        # setup histograms
        #setupHist(hist, labels, title, color, y_min, y_max):
        setupHist(h_mc,   cutList, key, blue_color, y_min, y_max)
        setupHist(h_data, cutList, key, red_color,  y_min, y_max)
        
        # draw histograms
        h_mc.Draw("hist")
        h_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        #legend = ROOT.TLegend()
        legend.AddEntry(h_data, "Data", "l")
        legend.AddEntry(h_mc,   "MC",   "l")
        legend.Draw()
        
        c.SetLogy()
        c.Update()
        c.SaveAs(plot_dir + key + ".pdf")

    for key in plot_map:
        cutList         = plot_map[key]["cuts"]
        key_Electron    = plot_map[key]["Electron"]
        key_Muon        = plot_map[key]["Muon"]
        f_Electron      = h_map[key_Electron]["file"]
        f_Muon          = h_map[key_Muon]["file"]
        # histograms
        h_Electron_mc   = f_Electron.Get(h_dir + h_map[key_Electron]["mc"])
        h_Electron_data = f_Electron.Get(h_dir + h_map[key_Electron]["data"])
        h_Muon_mc       = f_Muon.Get(h_dir + h_map[key_Muon]["mc"])
        h_Muon_data     = f_Muon.Get(h_dir + h_map[key_Muon]["data"])
        # ratios
        h_ratio_mc          = h_Electron_mc.Clone("h_ratio_mc")
        h_ratio_data        = h_Electron_data.Clone("h_ratio_data")
        h_ratio_Electron    = h_Electron_data.Clone("h_ratio_Electron")
        h_ratio_Muon        = h_Muon_data.Clone("h_ratio_Muon")
        h_ratio_mc.Divide(h_Muon_mc)
        h_ratio_data.Divide(h_Muon_data)
        h_ratio_Electron.Divide(h_Electron_mc)
        h_ratio_Muon.Divide(h_Muon_mc)
        
        # setup histograms
        #setupHist(hist, labels, title, color, y_min, y_max):
        setupHist(h_Electron_mc,    cutList, key, blue_color,   y_min, y_max)
        setupHist(h_Electron_data,  cutList, key, red_color,    y_min, y_max)
        setupHist(h_Muon_mc,        cutList, key, green_color,  y_min, y_max)
        setupHist(h_Muon_data,      cutList, key, purple_color, y_min, y_max)
        setupHist(h_ratio_mc,       cutList, key, blue_color,   0.0, 2.0)
        setupHist(h_ratio_data,     cutList, key, red_color,    0.0, 2.0)
        setupHist(h_ratio_Electron, cutList, key, red_color,    0.0, 2.0)
        setupHist(h_ratio_Muon,     cutList, key, purple_color, 0.0, 2.0)
        
        # --- draw histograms --- #
        h_Electron_mc.Draw("hist")
        h_Electron_data.Draw("hist same")
        h_Muon_mc.Draw("hist same")
        h_Muon_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        legend.AddEntry(h_Electron_data, "Electron Data", "l")
        legend.AddEntry(h_Electron_mc,   "Electron MC",   "l")
        legend.AddEntry(h_Muon_data,     "Muon Data",     "l")
        legend.AddEntry(h_Muon_mc,       "Muon MC",       "l")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_dir + key + ".pdf")

        # --- draw histograms --- #
        h_ratio_mc.Draw("hist")
        h_ratio_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        legend.AddEntry(h_ratio_data, "Electron/Muon Data", "l")
        legend.AddEntry(h_ratio_mc,   "Electron/Muon MC",   "l")
        legend.Draw()
        
        c.SetLogy(0) # unset log y
        c.Update()
        c.SaveAs(plot_dir + key + "_ElectronMuonRatios.pdf")
        
        # --- draw histograms --- #
        h_ratio_Electron.Draw("hist")
        h_ratio_Muon.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        legend.AddEntry(h_ratio_Electron, "Electron Data/MC", "l")
        legend.AddEntry(h_ratio_Muon,     "Muon Data/MC",     "l")
        legend.Draw()
        
        c.SetLogy(0) # unset log y
        c.Update()
        c.SaveAs(plot_dir + key + "_DataMCRatios.pdf")

        # log scale of same ratio plot
        # set y axis for log y
        setupHist(h_ratio_Electron, cutList, key, red_color,    0.1, 100.0)
        setupHist(h_ratio_Muon,     cutList, key, purple_color, 0.1, 100.0)
        
        # --- draw histograms --- #
        h_ratio_Electron.Draw("hist")
        h_ratio_Muon.Draw("hist same")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_dir + key + "_DataMCRatios_LogScale.pdf")

def main():
    #makePlots("condor/DataMC_2016_submission_2019-05-05_21-57-41/result.root", "2016")
    #makePlots("condor/DataMC_2017_submission_2019-05-05_22-28-09/result.root", "2017")
    #makePlots("condor/DataMC_2018_submission_2019-05-05_22-44-26/result.root", "2018")
    makePlots("condor/DataMC_2018_submission_2019-05-09_12-18-12/result.root", "2018")

if __name__ == "__main__":
    main()



