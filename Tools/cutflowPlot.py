# cutflowPlot.py

import ROOT
from colors import getColorIndex

def labelBins(hist, labels):
    x_axis = hist.GetXaxis()
    for i, label in enumerate(labels, 1):
        #print i, label
        x_axis.SetBinLabel(i, label)

def main():
    # don't display canvases while running
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    #f_name = "quickResult_v2.root"
    f_name = "ElectronCutFlow_v4.root"
    plot_dir = "cutflows/"
    h_dir = "CutFlows/"
    c = ROOT.TCanvas("c", "c", 800, 800)
    f = ROOT.TFile(f_name)
    if not f:
        print "Unable to load the file {0}".format(f_name)
        return
    # Data: apply trigger
    data_cuts = [
        "no cuts",
        "trigger",
        "filter",
        "MET",
        "HT",
        "nJets",
        "dPhi",
        "baseline",
        "Z to ee"
    ]
    # MC: no trigger, use placeholder
    mc_cuts = [
        "no cuts",
        "no cuts",
        "filter",
        "MET",
        "HT",
        "nJets",
        "dPhi",
        "baseline",
        "Z to ee"
    ]

    h_map = {
        "CutFlow_Electron_LowDM"  : {"data" : "CutFlow_Data_Electron_LowDM_met",  "mc" : "CutFlow_MC_Electron_LowDM_met"},
        "CutFlow_Electron_HighDM" : {"data" : "CutFlow_Data_Electron_HighDM_met", "mc" : "CutFlow_MC_Electron_HighDM_met"}
    }
    
    for key in h_map:
        h_mc   = f.Get(h_dir + h_map[key]["mc"])
        h_data = f.Get(h_dir + h_map[key]["data"])
        
        # setup histograms
        labelBins(h_mc, data_cuts)
        h_mc.SetStats(ROOT.kFALSE)
        h_mc.GetYaxis().SetRangeUser(0.1, 10.0**9)
        h_mc.SetLineColor(getColorIndex("blue"))
        h_data.SetLineColor(getColorIndex("red"))
        h_mc.SetLineWidth(3)
        h_data.SetLineWidth(3)
        
        # draw histograms
        h_mc.Draw("hist")
        h_data.Draw("hist same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(0.7, 0.9, 0.7, 0.9)
        #legend = ROOT.TLegend()
        legend.AddEntry(h_mc,   "MC",  "l")
        legend.AddEntry(h_data, "Data","l")
        legend.Draw()
        
        c.SetLogy()
        c.Update()
        c.SaveAs(plot_dir + key + ".pdf")


if __name__ == "__main__":
    main()

