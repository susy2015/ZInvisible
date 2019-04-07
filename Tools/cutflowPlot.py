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

def main():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    f_name_Electron = "quickResult_Electron_v2.root"
    f_name_Muon     = "quickResult_Muon_v1.root"
    #f_name_Electron = "ElectronCutFlow_v5.root"
    #f_name_Muon     = "MuonCutFlow_v1.root"
    plot_dir = "cutflows/"
    h_dir = "CutFlows/"
    f_names = [f_name_Electron, f_name_Muon]
    for f_name in f_names:
        exists = os.path.isfile(f_name)
        if not exists: 
            print "The file {0} does not exist".format(f_name)
            return
    
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
        setupHist(h_mc,   cutList, key, "blue", 0.1, 10.0**9)
        setupHist(h_data, cutList, key, "red",  0.1, 10.0**9)
        
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
        h_ratio_mc   = h_Electron_mc.Clone("h_ratio_mc")
        h_ratio_data = h_Electron_data.Clone("h_ratio_data")
        h_ratio_mc.Divide(h_Muon_mc)
        h_ratio_data.Divide(h_Muon_data)
        
        # setup histograms
        #setupHist(hist, labels, title, color, y_min, y_max):
        setupHist(h_Electron_mc,   cutList, key, "blue",   0.1, 10.0**9)
        setupHist(h_Electron_data, cutList, key, "red",    0.1, 10.0**9)
        setupHist(h_Muon_mc,       cutList, key, "green",  0.1, 10.0**9)
        setupHist(h_Muon_data,     cutList, key, "purple", 0.1, 10.0**9)
        setupHist(h_ratio_mc,      cutList, key, "blue",   0.0, 2.0)
        setupHist(h_ratio_data,    cutList, key, "red",    0.0, 2.0)
        
        # draw histograms
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

        # draw ratios
        h_ratio_mc.Draw("hist")
        h_ratio_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        legend.AddEntry(h_ratio_data, "Electron/Muon Data", "l")
        legend.AddEntry(h_ratio_mc,   "Electron/Muon MC",   "l")
        legend.Draw()
        
        c.SetLogy(0) # unset log y
        c.Update()
        c.SaveAs(plot_dir + key + "_ratios.pdf")



if __name__ == "__main__":
    main()







