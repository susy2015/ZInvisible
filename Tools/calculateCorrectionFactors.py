# calculateCorrectionFactors.py

import ROOT

# file: condor/histos_PhotonAndMuonControlRegionSelection_2016_11_Mar_2019_2/result.root
#
# directory: cutPhotonPt
# Low DM
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtDatadata
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPt#gamma+jetsstack
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtQCDstack
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtW(l#nu)+jetsstack
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtt#bar{t}stack
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPttWstack
# KEY: TH1D    DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtRarestack
# High DM
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPtDatadata
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPt#gamma+jetsstack
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPtQCDstack
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPtW(l#nu)+jetsstack
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPtt#bar{t}stack
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPttWstack
# KEY: TH1D    DataMC_Photon_HighDM_PhotonPt_2016cutPhotonPtcutPhotonPtRarestack

# directory: metWithPhoton
# Low DM
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonDatadata
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonQCDstack
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonW(l#nu)+jetsstack
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotont#bar{t}stack
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotontWstack
# KEY: TH1D  DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonRarestack
# High DM
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonDatadata
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonQCDstack
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonW(l#nu)+jetsstack
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotont#bar{t}stack
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotontWstack
# KEY: TH1D  DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonRarestack


class Histogram:

    def __init__(self, directory, name):
        self.directory = directory
        self.name = name
        self.histogram = 0

def main():
    plot_dir = "correctionFactors/"
    c = ROOT.TCanvas("c", "c", 800, 800)
    f_name = "condor/histos_PhotonAndMuonControlRegionSelection_2016_11_Mar_2019_2/result.root"
    f = ROOT.TFile(f_name)
    if not f:
        print "Unable to load the file {0}".format(f_name)
        return
    # map simple keys to full histogram name
    h_map = {}
    # --- cutPhotonPt --- #
    h_map["cutPhotonPt_Data_LowDM"]   = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtDatadata"          )
    h_map["cutPhotonPt_GJets_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPt#gamma+jetsstack"  )
    h_map["cutPhotonPt_QCD_LowDM"]    = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtQCDstack"          )
    h_map["cutPhotonPt_WJets_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtW(l#nu)+jetsstack" )
    h_map["cutPhotonPt_TTbar_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtt#bar{t}stack"     )
    h_map["cutPhotonPt_tW_LowDM"]     = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPttWstack"           )
    h_map["cutPhotonPt_Rare_LowDM"]   = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtRarestack"         )
    # --- metWithPhoton --- #
    h_map["metWithPhoton_Data_LowDM"]   = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonDatadata"          )
    h_map["metWithPhoton_GJets_LowDM"]  = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack"  )
    h_map["metWithPhoton_QCD_LowDM"]    = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonQCDstack"          )
    h_map["metWithPhoton_WJets_LowDM"]  = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonW(l#nu)+jetsstack" )
    h_map["metWithPhoton_TTbar_LowDM"]  = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotont#bar{t}stack"     )
    h_map["metWithPhoton_tW_LowDM"]     = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotontWstack"           )
    h_map["metWithPhoton_Rare_LowDM"]   = Histogram("metWithPhoton", "DataMC_Photon_HighDM_met_2016metWithPhotonmetWithPhotonRarestack"         )
    for h_key in h_map:
        myHist = h_map[h_key]
        myHist.key = h_key
        myHist.histogram = f.Get(myHist.directory + "/" + myHist.name)
        #print "{0} histogram: type={1}, n_bins={2}".format(myHist.key, type(myHist.histogram), myHist.histogram.GetNbinsX())
        if not myHist.histogram:
            print "Unable to load the histogram for {0} named {1}".format(myHist.key, myHist.name)

    variables = ["cutPhotonPt", "metWithPhoton"]
    for var in variables:
        # Data
        h_Data_LowDM = h_map[var + "_Data_LowDM"].histogram.Clone(var + "_Data_LowDM")
        # MC
        h_MC_LowDM   = h_map[var + "_GJets_LowDM"].histogram.Clone(var + "_GJets_LowDM")
        h_MC_LowDM.Add(h_map[var + "_QCD_LowDM"].histogram)
        h_MC_LowDM.Add(h_map[var + "_WJets_LowDM"].histogram)
        h_MC_LowDM.Add(h_map[var + "_TTbar_LowDM"].histogram)
        h_MC_LowDM.Add(h_map[var + "_tW_LowDM"].histogram)
        h_MC_LowDM.Add(h_map[var + "_Rare_LowDM"].histogram)
        # GJets
        h_GJets_LowDM = h_map[var + "_GJets_LowDM"].histogram.Clone(var + "_GJets_LowDM")
        # Background MC (without GJets)
        h_Background_LowDM   = h_map[var + "_QCD_LowDM"].histogram.Clone(var + "_QCD_LowDM")
        h_Background_LowDM.Add(h_map[var + "_WJets_LowDM"].histogram)
        h_Background_LowDM.Add(h_map[var + "_TTbar_LowDM"].histogram)
        h_Background_LowDM.Add(h_map[var + "_tW_LowDM"].histogram)
        h_Background_LowDM.Add(h_map[var + "_Rare_LowDM"].histogram)
        # Data minus background
        h_DataMinusBackground_LowDM = h_map[var + "_Data_LowDM"].histogram.Clone(var + "_Data_LowDM")
        h_DataMinusBackground_LowDM.Add(h_Background_LowDM, -1)

        h_MC_LowDM.Draw()
        h_Data_LowDM.Draw("same")
        c.Update()
        c.SaveAs(plot_dir + var + "_Data_MC.pdf")
        
        h_GJets_LowDM.Draw()
        h_DataMinusBackground_LowDM.Draw("same")
        c.Update()
        c.SaveAs(plot_dir + var + "_DataMinusBackground_GJets.pdf")

if __name__ == "__main__":
    main()




