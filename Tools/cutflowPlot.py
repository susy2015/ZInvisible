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

def makePlots(file_name, year, doPhotons, useHEMVeto):
    yearTag = "_" + year
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    plot_dir = "cutflows/"
    h_dir = "CutFlows/"
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
    cuts = []
    if useHEMVeto:
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
    else:
        cutsLepton = [
            "No cuts",
            "Trigger",
            "Lepton veto",
            "LL pt",
            "LL charge",
            "m_LL > 50",
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
            "JetID",
            "Event filter",
            "original MET < 250",
            "modified MET > 250",
            "HT",
            "nJets",
            "dPhi",
        ]


    # copy list instead of pointer to list
    cutsLeptonLowDM = cutsLepton[:]
    cutsLeptonLowDM.append("Low #Deltam")
    cutsLeptonHighDM = cutsLepton[:]
    cutsLeptonHighDM.append("High #Deltam")
    cutsPhotonLowDM = cutsPhoton[:]
    cutsPhotonLowDM.append("Low #Deltam")
    cutsPhotonHighDM = cutsPhoton[:]
    cutsPhotonHighDM.append("High #Deltam")
    
    h_map = {}
    if doPhotons:
        h_map = {
            "CutFlow_Electron_LowDM"  : {"data" : "CutFlow_Data_Electron_LowDM_met",  "mc" : "CutFlow_MC_Electron_LowDM_met",  "cuts" : cutsLeptonLowDM},
            "CutFlow_Electron_HighDM" : {"data" : "CutFlow_Data_Electron_HighDM_met", "mc" : "CutFlow_MC_Electron_HighDM_met", "cuts" : cutsLeptonHighDM},
            "CutFlow_Muon_LowDM"  : {"data" : "CutFlow_Data_Muon_LowDM_met",  "mc" : "CutFlow_MC_Muon_LowDM_met",  "cuts" : cutsLeptonLowDM},
            "CutFlow_Muon_HighDM" : {"data" : "CutFlow_Data_Muon_HighDM_met", "mc" : "CutFlow_MC_Muon_HighDM_met", "cuts" : cutsLeptonHighDM},
            "CutFlow_Photon_LowDM"  : {"data" : "CutFlow_Data_Photon_LowDM_met",  "mc" : "CutFlow_MC_Photon_LowDM_met",  "cuts" : cutsPhotonLowDM},
            "CutFlow_Photon_HighDM" : {"data" : "CutFlow_Data_Photon_HighDM_met", "mc" : "CutFlow_MC_Photon_HighDM_met", "cuts" : cutsPhotonHighDM}
        }
    else:
        h_map = {
            "CutFlow_Electron_LowDM"  : {"data" : "CutFlow_Data_Electron_LowDM_met",  "mc" : "CutFlow_MC_Electron_LowDM_met",  "cuts" : cutsLeptonLowDM},
            "CutFlow_Electron_HighDM" : {"data" : "CutFlow_Data_Electron_HighDM_met", "mc" : "CutFlow_MC_Electron_HighDM_met", "cuts" : cutsLeptonHighDM},
            "CutFlow_Muon_LowDM"  : {"data" : "CutFlow_Data_Muon_LowDM_met",  "mc" : "CutFlow_MC_Muon_LowDM_met",  "cuts" : cutsLeptonLowDM},
            "CutFlow_Muon_HighDM" : {"data" : "CutFlow_Data_Muon_HighDM_met", "mc" : "CutFlow_MC_Muon_HighDM_met", "cuts" : cutsLeptonHighDM}
        }
    plot_map = {
        "CutFlow_LowDM"  : {"Electron" : "CutFlow_Electron_LowDM",  "Muon" : "CutFlow_Muon_LowDM",  "cuts" :  cutsLeptonLowDM},
        "CutFlow_HighDM" : {"Electron" : "CutFlow_Electron_HighDM", "Muon" : "CutFlow_Muon_HighDM", "cuts" :  cutsLeptonHighDM}
    }
    
    for key in h_map:
        plot_name = plot_dir + key
        cutList = h_map[key]["cuts"]    
        h_mc_name   = h_dir + h_map[key]["mc"]
        h_data_name = h_dir + h_map[key]["data"]
        h_mc        = rootFile.Get(h_mc_name)
        h_data      = rootFile.Get(h_data_name)
        if not h_mc:
            print "ERROR: Unable to load histogram {0}".format(h_mc_name)
            exit(1)
        if not h_data:
            print "ERROR: Unable to load histogram {0}".format(h_data_name)
            exit(1)
        
        # setup histograms
        #setupHist(hist, labels, title, color, y_min, y_max):
        setupHist(h_mc,   cutList, key + yearTag, blue_color, y_min, y_max)
        setupHist(h_data, cutList, key + yearTag, red_color,  y_min, y_max)
        
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
        c.SaveAs(plot_name + yearTag + ".pdf")
        c.SaveAs(plot_name + yearTag + ".png")

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
        setupHist(h_Electron_mc,    cutList, key + yearTag, blue_color,   y_min, y_max)
        setupHist(h_Electron_data,  cutList, key + yearTag, red_color,    y_min, y_max)
        setupHist(h_Muon_mc,        cutList, key + yearTag, green_color,  y_min, y_max)
        setupHist(h_Muon_data,      cutList, key + yearTag, purple_color, y_min, y_max)
        setupHist(h_ratio_mc,       cutList, key + yearTag, blue_color,   0.0, 2.0)
        setupHist(h_ratio_data,     cutList, key + yearTag, red_color,    0.0, 2.0)
        setupHist(h_ratio_Electron, cutList, key + yearTag, red_color,    0.0, 2.0)
        setupHist(h_ratio_Muon,     cutList, key + yearTag, purple_color, 0.0, 2.0)
        
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
        c.SaveAs(plot_name + yearTag + ".pdf")
        c.SaveAs(plot_name + yearTag + ".png")

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
        c.SaveAs(plot_name + "_ElectronMuonRatios" + yearTag + ".pdf")
        c.SaveAs(plot_name + "_ElectronMuonRatios" + yearTag + ".png")
        
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
        c.SaveAs(plot_name + "_DataMCRatios" + yearTag + ".pdf")
        c.SaveAs(plot_name + "_DataMCRatios" + yearTag + ".png")

        # log scale of same ratio plot
        # set y axis for log y
        setupHist(h_ratio_Electron, cutList, key + yearTag, red_color,    0.1, 100.0)
        setupHist(h_ratio_Muon,     cutList, key + yearTag, purple_color, 0.1, 100.0)
        
        # --- draw histograms --- #
        h_ratio_Electron.Draw("hist")
        h_ratio_Muon.Draw("hist same")
        legend.Draw()
        
        c.SetLogy(1) # set log y
        c.Update()
        c.SaveAs(plot_name + "_DataMCRatios_LogScale" + yearTag + ".pdf")
        c.SaveAs(plot_name + "_DataMCRatios_LogScale" + yearTag + ".png")

def main():
    version = 4
    # set per version
    useHEMVeto_map = {1: False, 2:False, 3:True,  4:True}
    doPhotons_map  = {1: False, 2:False, 3:False, 4:True}
    useHEMVeto     = useHEMVeto_map[version]
    doPhotons      = doPhotons_map[version]
    if version == 1:
        # May 5, 2019 Results
        makePlots("condor/DataMC_2016_submission_2019-05-05_21-57-41/result.root", "2016", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2017_submission_2019-05-05_22-28-09/result.root", "2017", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_submission_2019-05-05_22-44-26/result.root", "2018", doPhotons, useHEMVeto)
    if version == 2:
        # May 9, 2019 Results
        makePlots("condor/DataMC_2016_submission_2019-05-09_17-19-42/result.root", "2016", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2017_submission_2019-05-09_17-16-54/result.root", "2017", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_submission_2019-05-09_17-15-04/result.root", "2018", doPhotons, useHEMVeto)
    if version == 3:
        # May 16, 2019 Results
        makePlots("condor/DataMC_2016_submission_2019-05-16_10-06-59/result.root",    "2016",    doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2017_submission_2019-05-16_10-09-29/result.root",    "2017",    doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_AB_submission_2019-05-16_10-10-30/result.root", "2018_AB", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_CD_submission_2019-05-16_10-12-04/result.root", "2018_CD", doPhotons, useHEMVeto)
    if version == 4:
        # May 17, 2019 Results
        makePlots("condor/DataMC_2016_submission_2019-05-17_18-46-29/result.root",    "2016",    doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2017_submission_2019-05-17_18-47-28/result.root",    "2017",    doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_AB_submission_2019-05-17_18-48-10/result.root", "2018_AB", doPhotons, useHEMVeto)
        makePlots("condor/DataMC_2018_CD_submission_2019-05-17_18-50-16/result.root", "2018_CD", doPhotons, useHEMVeto)

if __name__ == "__main__":
    main()


