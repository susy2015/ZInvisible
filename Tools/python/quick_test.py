# make_table.py
import ROOT
import os

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Table:
    def __init__(self, out_file):
        self.out_file = out_file
        self.jetpt_cut  - "jetpt30"
        self.regions    = ["LowDM"]
        self.selections = ["nb0", "nb1"]
        self.variable   = "nJets_drPhotonCleaned_" + self.jetpt_cut
        self.histos = {}

    def setupHistoMap(self, era):
        # histogram examples
        # DataMC_Photon_LowDM_nj_nb0_rebin_jetpt30_2016nJets_drPhotonCleaned_jetpt30nJets_drPhotonCleaned_jetpt30Datadata
        # DataMC_Photon_LowDM_nj_nb1_rebin_jetpt30_2016nJets_drPhotonCleaned_jetpt30nJets_drPhotonCleaned_jetpt30Datadata
        self.histos[era] = {}
        for region in self.regions:
            self.histos[era][region] = {}
            for selection in self.selections:
                self.histos[era][region][selection] = {
                        "Data"  : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "Datadata",
                        "GJets" : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "#gamma+jetsstack",
                        "QCD"   : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "QCDstack",
                        "WJets" : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "W(l#nu)+jetsstack",
                        "TTG"   : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "t#bar{t}#gamma+jetsstack",
                        "TTbar" : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "t#bar{t}stack",
                        "tW"    : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "tWstack",
                        "Rare"  : "DataMC_Photon_" + region + "_nj_" + selection + "_rebin_" + self.jetpt_cut + "_" + era + 2 * self.variable + "Rarestack",
                }

    def makeTable(self, in_file, era):
        print ">>> Making table for {0}".format(era)
        # check that the file exists
        if not os.path.isfile(in_file): 
            print "ERROR: The file {0} does not exist".format(in_file)
            return
        
        f = ROOT.TFile(in_file)
    
        # setup histogram map
        self.setupHistoMap(era)
        
        for region in self.regions:
            for selection in self.selections:
                h_Data  = f.Get(self.variable + "/" + self.histos[era][region][selection]["Data"])
                h_GJets = f.Get(self.variable + "/" + self.histos[era][region][selection]["GJets"])
                h_QCD   = f.Get(self.variable + "/" + self.histos[era][region][selection]["QCD"])
                h_WJets = f.Get(self.variable + "/" + self.histos[era][region][selection]["WJets"])
                h_TTG   = f.Get(self.variable + "/" + self.histos[era][region][selection]["TTG"])
                h_TTbar = f.Get(self.variable + "/" + self.histos[era][region][selection]["TTbar"])
                h_tW    = f.Get(self.variable + "/" + self.histos[era][region][selection]["tW"])
                h_Rare  = f.Get(self.variable + "/" + self.histos[era][region][selection]["Rare"])
                
                # MC_background
                h_back = h_GJets.Clone("h_back")
                h_back.Add(h_QCD)
                h_back.Add(h_WJets)
                h_back.Add(h_TTG)
                h_back.Add(h_TTbar)
                h_back.Add(h_tW)
                h_back.Add(h_Rare)
                
                h_num = h_Data.Clone("h_num")
                h_den = h_back.Clone("h_den") 
                
                # number of events for normalization
                nNum  = h_num.Integral(0, h_num.GetNbinsX() + 1)
                nDen  = h_den.Integral(0, h_den.GetNbinsX() + 1)
                ratio = float(nNum) / float(nDen)
                
                h_den_normalized = h_den.Clone("h_den_normalized")
                h_den_normalized.Scale(ratio)

                h_ratio_normalized = h_num.Clone("h_ratio_normalized")
                h_ratio_normalized.Divide(h_den_normalized)
                
                self.out_file.write("{0} {1} {2}\n".format(era, region, selection))
                for b in xrange(1, h_ratio_normalized.GetNbinsX() +1):
                    self.out_file.write("bin {0}: {1:.4f}\n".format(b, h_ratio_normalized.GetBinContent(b)))




