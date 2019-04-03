import ROOT

def labelBins(hist, labels):
    x_axis = hist.GetXaxis()
    for i, label in enumerate(labels, 1):
        #print i, label
        x_axis.SetBinLabel(i, label)

def main():
    # don't display canvases while running
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    f_name = "quickResult.root"
    #f_name = "ElectronCutFlow_v2.root"
    plot_dir = "cutflows/"
    h_dir = "CutFlows/"
    c = ROOT.TCanvas("c", "c", 800, 800)
    f = ROOT.TFile(f_name)
    if not f:
        print "Unable to load the file {0}".format(f_name)
        return
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
    mc_cuts = [
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
        "CutFlow_MC_Electron_LowDM_met" : mc_cuts,
        "CutFlow_MC_Electron_HighDM_met" : mc_cuts,
        "CutFlow_Data_Electron_LowDM_met" : data_cuts,
        "CutFlow_Data_Electron_HighDM_met" : data_cuts
    }
    for h_name in h_map:
        h = f.Get(h_dir + h_name)
        labelBins(h, h_map[h_name])
        h.GetYaxis().SetRangeUser(0.1, 10.0**9)
        h.Draw("hist")
        c.SetLogy()
        c.Update()
        c.SaveAs(plot_dir + h_name + ".pdf")


if __name__ == "__main__":
    main()

