# calculateCorrectionFactors.py

import ROOT

# file: condor/histos_PhotonAndMuonControlRegionSelection_2016_11_Mar_2019_2/result.root
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

class Histogram:

    def __init__(self, directory, name):
        self.directory = directory
        self.name = name
        self.histogram = 0

def main():
    c = ROOT.TCanvas("c", "c", 800, 800)
    f_name = "condor/histos_PhotonAndMuonControlRegionSelection_2016_11_Mar_2019_2/result.root"
    f = ROOT.TFile(f_name)
    if not f:
        print "Unable to load the file {0}".format(f_name)
        return
    # map simple keys to full histogram name
    h_map = {}
    h_map["PhotonPt_Data_LowDM"]   = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtDatadata"          )
    h_map["PhotonPt_GJets_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPt#gamma+jetsstack"  )
    h_map["PhotonPt_QCD_LowDM"]    = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtQCDstack"          )
    h_map["PhotonPt_WJets_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtW(l#nu)+jetsstack" )
    h_map["PhotonPt_TTbar_LowDM"]  = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtt#bar{t}stack"     )
    h_map["PhotonPt_tW_LowDM"]     = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPttWstack"           )
    h_map["PhotonPt_Rare"]         = Histogram("cutPhotonPt", "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtRarestack"         )
    for h_key in h_map:
        myHist = h_map[h_key]
        myHist.key = h_key
        myHist.histogram = f.Get(myHist.directory + "/" + myHist.name)
        #myHist.histogram = ROOT.TH1F(f.Get(myHist.directory + "/" + myHist.name))
        print "{0} histogram: type={1}, n_bins={2}".format(myHist.key, type(myHist.histogram), myHist.histogram.GetNbinsX())
        if not myHist.histogram:
            print "Unable to load the histogram for {0} named {1}".format(myHist.key, myHist.name)
        #myHist.histogram.Draw()
        #c.Update()



    #h_map["PhotonPt_Data_LowDM"].histogram.Draw()

    # Data
    h_PhotonPt_Data_LowDM = h_map["PhotonPt_Data_LowDM"].histogram.Clone("PhotonPt_Data_LowDM")
    # MC
    h_PhotonPt_MC_LowDM   = h_map["PhotonPt_GJets_LowDM"].histogram.Clone("PhotonPt_GJets_LowDM")
    h_PhotonPt_MC_LowDM.Add(h_map["PhotonPt_QCD_LowDM"].histogram)
    h_PhotonPt_MC_LowDM.Add(h_map["PhotonPt_WJets_LowDM"].histogram)
    h_PhotonPt_MC_LowDM.Add(h_map["PhotonPt_TTbar_LowDM"].histogram)
    h_PhotonPt_MC_LowDM.Add(h_map["PhotonPt_tW_LowDM"].histogram)
    h_PhotonPt_MC_LowDM.Add(h_map["PhotonPt_Rare"].histogram)
    # GJets
    h_PhotonPt_GJets_LowDM = h_map["PhotonPt_GJets_LowDM"].histogram.Clone("PhotonPt_GJets_LowDM")
    # Background MC (without GJets)
    h_PhotonPt_Background_LowDM   = h_map["PhotonPt_QCD_LowDM"].histogram.Clone("PhotonPt_QCD_LowDM")
    h_PhotonPt_Background_LowDM.Add(h_map["PhotonPt_WJets_LowDM"].histogram)
    h_PhotonPt_Background_LowDM.Add(h_map["PhotonPt_TTbar_LowDM"].histogram)
    h_PhotonPt_Background_LowDM.Add(h_map["PhotonPt_tW_LowDM"].histogram)
    h_PhotonPt_Background_LowDM.Add(h_map["PhotonPt_Rare"].histogram)
    # Data minus background
    h_PhotonPt_DataMinusBackground_LowDM = h_map["PhotonPt_Data_LowDM"].histogram.Clone("PhotonPt_Data_LowDM")
    h_PhotonPt_DataMinusBackground_LowDM.Add(h_PhotonPt_Background_LowDM, -1)
 


    if not h_PhotonPt_Data_LowDM:
        print "ERROR: h_PhotonPt_Data_LowDM not created"
    if not h_PhotonPt_MC_LowDM:
        print "ERROR: h_PhotonPt_MC_LowDM not created"

    print "h_PhotonPt_Data_LowDM histogram: type={0}, n_bins={1}".format(type(h_PhotonPt_Data_LowDM), h_PhotonPt_Data_LowDM.GetNbinsX())
    print "h_PhotonPt_MC_LowDM histogram: type={0}, n_bins={1}".format(type(h_PhotonPt_MC_LowDM), h_PhotonPt_MC_LowDM.GetNbinsX())
    
    h_PhotonPt_MC_LowDM.Draw()
    h_PhotonPt_Data_LowDM.Draw("same")
    c.Update()
    c.SaveAs("PhotonPt_Data_MC.pdf")
    
    h_PhotonPt_GJets_LowDM.Draw()
    h_PhotonPt_DataMinusBackground_LowDM.Draw("same")
    c.Update()
    c.SaveAs("PhotonPt_DataMinusBackground_GJets.pdf")

def test():
    myC = ROOT.TCanvas("myC", "myC", 800, 800)
    myC.cd()
    f_name = "condor/histos_PhotonAndMuonControlRegionSelection_2016_11_Mar_2019_2/result.root"
    f = ROOT.TFile(f_name)
    histogram = f.Get("cutPhotonPt" + "/" + "DataMC_Photon_LowDM_PhotonPt_2016cutPhotonPtcutPhotonPtDatadata")
    print histogram.GetNbinsX()
    histogram.Draw()
    myC.Update()
    #myC.SaveAs("myC.pdf")

if __name__ == "__main__":
    main()
    #test()




