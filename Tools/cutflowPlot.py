# cutflowPlot.py

import os
import ROOT
from colors import getColorIndex

def setupHist(hist, labels, title, color):
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    # label bins
    for i, label in enumerate(labels, 1):
        #print i, label
        x_axis.SetBinLabel(i, label)
    
    y_axis.SetRangeUser(0.1, 10.0**9)
    hist.SetTitle(title)
    hist.SetStats(ROOT.kFALSE)
    hist.SetLineColor(getColorIndex(color))
    hist.SetLineWidth(3)

def main():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    f_name_Electron = "quickResult_Electron_v1.root"
    f_name_Muon     = "quickResult_Muon_v1.root"
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
    
    for key in h_map:
        cutList = h_map[key]["cuts"]    
        f       = h_map[key]["file"]
        h_mc    = f.Get(h_dir + h_map[key]["mc"])
        h_data  = f.Get(h_dir + h_map[key]["data"])
        
        # setup histograms
        #setupHist(hist, labels, title, color):
        setupHist(h_mc,   cutList, key, "blue")
        setupHist(h_data, cutList, key, "red")
        
        # draw histograms
        h_mc.Draw("hist")
        h_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        #legend = ROOT.TLegend()
        legend.AddEntry(h_data, "Data","l")
        legend.AddEntry(h_mc,   "MC",  "l")
        legend.Draw()
        
        c.SetLogy()
        c.Update()
        c.SaveAs(plot_dir + key + ".pdf")


if __name__ == "__main__":
    main()

