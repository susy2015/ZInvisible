# cutflow_plot.py

import os
import ROOT
from colors import getColorIndex
    
# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

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

def makeCutflows(file_name, era, plot_dir, doPhotons):
    eraTag = "_" + era
    h_dir = "CutFlows/"
    if plot_dir[-1] != "/":
        plot_dir += "/"
    if h_dir[-1] != "/":
        h_dir += "/"
    # check that the file exists
    if not os.path.isfile(file_name): 
        print "The file {0} does not exist".format(file_name)
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

    rootFile = ROOT.TFile(file_name)
    c = ROOT.TCanvas("c", "c", 800, 800)
        
    # legend: TLegend(x1,y1,x2,y2)
    legend_x1 = 0.5
    legend_x2 = 0.9 
    legend_y1 = 0.7 
    legend_y2 = 0.9 

    cuts = []
    # Z NuNu in validation region
    cutsZNuNu = [
        "No cuts",
        "Lepton veto",
        "HEM veto",
        "JetID",
        "Event filter",
        "MET",
        "HT",
        "nJets",
    ]
    cutsLepton = [
        "No cuts",
        "Trigger",
        "Lepton veto",
        "LL pt",
        "LL charge",
        "m_LL > 50",
        "HEM veto",
        "JetID",
        "Event filter",
        "MET",
        "HT",
        "nJets",
        "dPhi",
    ]
    cutsPhoton = [
        "No cuts",
        "Trigger",
        "Lepton veto",
        "Photon p_{T}",
        "HEM veto",
        "JetID",
        "Event filter",
        "original MET < 250",
        "modified MET > 250",
        "HT",
        "nJets",
        "dPhi",
    ]

    # copy list instead of pointer to list
    cutsZNuNuLowDM = cutsZNuNu[:]
    cutsZNuNuLowDM.append("dPhi")
    cutsZNuNuLowDM.append("Low #Deltam")
    cutsZNuNuLowDMHighMET = cutsZNuNu[:]
    cutsZNuNuLowDMHighMET.append("mid dPhi")
    cutsZNuNuLowDMHighMET.append("Low #Deltam")
    cutsZNuNuHighDM = cutsZNuNu[:]
    cutsZNuNuHighDM.append("mid dPhi")
    cutsZNuNuHighDM.append("High #Deltam")
    cutsLeptonLowDM = cutsLepton[:]
    cutsLeptonLowDM.append("Low #Deltam")
    cutsLeptonHighDM = cutsLepton[:]
    cutsLeptonHighDM.append("High #Deltam")
    cutsPhotonLowDM = cutsPhoton[:]
    cutsPhotonLowDM.append("Low #Deltam")
    cutsPhotonHighDM = cutsPhoton[:]
    cutsPhotonHighDM.append("High #Deltam")

    # ZNuNu MC
    # CutFlow_MC_ZNuNu_met_LowDM
    # CutFlow_MC_ZNuNu_met_LowDM_HighMET
    # CutFlow_MC_ZNuNu_met_HighDM
    
    h_map = {
        "CutFlow_ZNuNu_LowDM"         : {"mc" : "CutFlow_MC_ZNuNu_met_LowDM",         "cuts" : cutsZNuNuLowDM},
        "CutFlow_ZNuNu_LowDM_HighMET" : {"mc" : "CutFlow_MC_ZNuNu_met_LowDM_HighMET", "cuts" : cutsZNuNuLowDMHighMET},
        "CutFlow_ZNuNu_HighDM"        : {"mc" : "CutFlow_MC_ZNuNu_met_HighDM",        "cuts" : cutsZNuNuHighDM},
        "CutFlow_Electron_LowDM"  : {"data" : "CutFlow_Data_Electron_LowDM_met",  "mc" : "CutFlow_MC_Electron_LowDM_met",  "cuts" : cutsLeptonLowDM},
        "CutFlow_Electron_HighDM" : {"data" : "CutFlow_Data_Electron_HighDM_met", "mc" : "CutFlow_MC_Electron_HighDM_met", "cuts" : cutsLeptonHighDM},
        "CutFlow_Muon_LowDM"  : {"data" : "CutFlow_Data_Muon_LowDM_met",  "mc" : "CutFlow_MC_Muon_LowDM_met",  "cuts" : cutsLeptonLowDM},
        "CutFlow_Muon_HighDM" : {"data" : "CutFlow_Data_Muon_HighDM_met", "mc" : "CutFlow_MC_Muon_HighDM_met", "cuts" : cutsLeptonHighDM}
    }
    if doPhotons:
        h_map["CutFlow_Photon_LowDM"]  = {"data" : "CutFlow_Data_Photon_LowDM_met",  "mc" : "CutFlow_MC_Photon_LowDM_met",  "cuts" : cutsPhotonLowDM}
        h_map["CutFlow_Photon_HighDM"] = {"data" : "CutFlow_Data_Photon_HighDM_met", "mc" : "CutFlow_MC_Photon_HighDM_met", "cuts" : cutsPhotonHighDM}
    plot_map = {
        "CutFlow_LowDM"  : {"Electron" : "CutFlow_Electron_LowDM",  "Muon" : "CutFlow_Muon_LowDM",  "cuts" :  cutsLeptonLowDM},
        "CutFlow_HighDM" : {"Electron" : "CutFlow_Electron_HighDM", "Muon" : "CutFlow_Muon_HighDM", "cuts" :  cutsLeptonHighDM}
    }
    
    for key in h_map:
        doData = False
        if "data" in h_map[key]:
            doData = True
        plot_name = plot_dir + key
        cutList = h_map[key]["cuts"]    
        h_mc_name   = h_dir + h_map[key]["mc"]
        h_mc        = rootFile.Get(h_mc_name)
        if not h_mc:
            print "ERROR: Unable to load histogram {0}".format(h_mc_name)
            exit(1)
        if doData:
            h_data_name = h_dir + h_map[key]["data"]
            h_data      = rootFile.Get(h_data_name)
            if not h_data:
                print "ERROR: Unable to load histogram {0}".format(h_data_name)
                exit(1)
        
        # setup histograms
        #setupHist(hist, labels, title, color, y_min, y_max):
        setupHist(h_mc,   cutList, key + eraTag, blue_color, y_min, y_max)
        if doData:
            setupHist(h_data, cutList, key + eraTag, red_color,  y_min, y_max)
        
        # draw histograms
        h_mc.Draw("hist")
        if doData:
            h_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        if doData:
            legend.AddEntry(h_data, "Data", "l")
        legend.AddEntry(h_mc,   "MC",   "l")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_name + eraTag + ".pdf")
        c.SaveAs(plot_name + eraTag + ".png")

    for key in plot_map:
        plot_name = plot_dir + key
        cutList         = plot_map[key]["cuts"]
        key_Electron    = plot_map[key]["Electron"]
        key_Muon        = plot_map[key]["Muon"]
        # histograms
        h_Electron_mc   = rootFile.Get(h_dir + h_map[key_Electron]["mc"])
        h_Electron_data = rootFile.Get(h_dir + h_map[key_Electron]["data"])
        h_Muon_mc       = rootFile.Get(h_dir + h_map[key_Muon]["mc"])
        h_Muon_data     = rootFile.Get(h_dir + h_map[key_Muon]["data"])
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
        setupHist(h_Electron_mc,    cutList, key + eraTag, blue_color,   y_min, y_max)
        setupHist(h_Electron_data,  cutList, key + eraTag, red_color,    y_min, y_max)
        setupHist(h_Muon_mc,        cutList, key + eraTag, green_color,  y_min, y_max)
        setupHist(h_Muon_data,      cutList, key + eraTag, purple_color, y_min, y_max)
        setupHist(h_ratio_mc,       cutList, key + eraTag, blue_color,   0.0, 2.0)
        setupHist(h_ratio_data,     cutList, key + eraTag, red_color,    0.0, 2.0)
        setupHist(h_ratio_Electron, cutList, key + eraTag, red_color,    0.0, 2.0)
        setupHist(h_ratio_Muon,     cutList, key + eraTag, purple_color, 0.0, 2.0)
        
        # --- draw histograms --- #
        h_Electron_mc.Draw("hist")
        h_Electron_data.Draw("hist same")
        h_Muon_mc.Draw("hist same")
        h_Muon_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        legend.AddEntry(h_Electron_data, "Electron Data", "l")
        legend.AddEntry(h_Electron_mc,   "Electron MC",   "l")
        legend.AddEntry(h_Muon_data,     "Muon Data",     "l")
        legend.AddEntry(h_Muon_mc,       "Muon MC",       "l")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_name + eraTag + ".pdf")
        c.SaveAs(plot_name + eraTag + ".png")

        # --- draw histograms --- #
        h_ratio_mc.Draw("hist")
        h_ratio_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        legend.AddEntry(h_ratio_data, "Electron/Muon Data", "l")
        legend.AddEntry(h_ratio_mc,   "Electron/Muon MC",   "l")
        legend.Draw()
        
        c.SetLogy(0) # unset log y
        c.Update()
        c.SaveAs(plot_name + "_ElectronMuonRatios" + eraTag + ".pdf")
        c.SaveAs(plot_name + "_ElectronMuonRatios" + eraTag + ".png")
        
        # --- draw histograms --- #
        h_ratio_Electron.Draw("hist")
        h_ratio_Muon.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        legend.AddEntry(h_ratio_Electron, "Electron Data/MC", "l")
        legend.AddEntry(h_ratio_Muon,     "Muon Data/MC",     "l")
        legend.Draw()
        
        c.SetLogy(0) # unset log y
        c.Update()
        c.SaveAs(plot_name + "_DataMCRatios" + eraTag + ".pdf")
        c.SaveAs(plot_name + "_DataMCRatios" + eraTag + ".png")

        # log scale of same ratio plot
        # set y axis for log y
        setupHist(h_ratio_Electron, cutList, key + eraTag, red_color,    0.1, 100.0)
        setupHist(h_ratio_Muon,     cutList, key + eraTag, purple_color, 0.1, 100.0)
        
        # --- draw histograms --- #
        h_ratio_Electron.Draw("hist")
        h_ratio_Muon.Draw("hist same")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_name + "_DataMCRatios_LogScale" + eraTag + ".pdf")
        c.SaveAs(plot_name + "_DataMCRatios_LogScale" + eraTag + ".png")

def main():
    plot_dir = "more_plots"
    # set per version
    doPhotons = True
    makeCutflows("histoutput.root",    "2016",    plot_dir, doPhotons)

if __name__ == "__main__":
    main()


