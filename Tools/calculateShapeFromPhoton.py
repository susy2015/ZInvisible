# calculateShapeFromPhoton.py

import os
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Shape:
    def __init__(self, verbose):
        self.verbose = verbose
        self.histos = {}
        self.regions   = ["LowDM", "HighDM"]
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "metWithPhoton"
        self.plot_dir = "met_shapes/"
        if self.plot_dir[-1] != "/":
            self.plot_dir += "/"

    def setupHistoMap(self, year):
        # histograms
        # example
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonDatadata
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack
        self.histos[year] = {}
        for region in self.regions:
            self.histos[year][region] = {
                    "Data"  : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotonDatadata",
                    "GJets" : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack",
                    "QCD"   : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotonQCDstack",
                    "WJets" : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotonW(l#nu)+jetsstack",
                    "TTbar" : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotont#bar{t}stack",
                    "tW"    : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotontWstack",
                    "Rare"  : "DataMC_Photon_" + region + "_met_2016metWithPhotonmetWithPhotonRarestack",
            }
    
    def getShape(self, file_name, era): 
        # currently the histograms are named by year (2018) and not era (2018_AB)
        # we should probalby change the histograms to use era (2018_AB)
        year = era[0:4]
        eraTag = "_" + era
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        # make directory for plots if it does not exist
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)
        
        f = ROOT.TFile(file_name)
        c = ROOT.TCanvas("c", "c", 800, 800)
        
        # setup histogram map
        self.setupHistoMap(year)
        
        for region in self.regions:
            plot_name = self.plot_dir + self.variable + "_" + region
            h_Data  = f.Get(self.variable + "/" + self.histos[year][region]["Data"])
            h_GJets = f.Get(self.variable + "/" + self.histos[year][region]["GJets"])
            h_QCD   = f.Get(self.variable + "/" + self.histos[year][region]["QCD"])
            h_WJets = f.Get(self.variable + "/" + self.histos[year][region]["WJets"])
            h_TTbar = f.Get(self.variable + "/" + self.histos[year][region]["TTbar"])
            h_tW    = f.Get(self.variable + "/" + self.histos[year][region]["tW"])
            h_Rare  = f.Get(self.variable + "/" + self.histos[year][region]["Rare"])
            h_MC = h_GJets.Clone("h_MC") 
            h_MC.Add(h_QCD)
            h_MC.Add(h_WJets)
            h_MC.Add(h_TTbar)
            h_MC.Add(h_tW)
            h_MC.Add(h_Rare)

            # number of data events
            nData = h_Data.Integral(0, h_Data.GetNbinsX() + 1)
            nMC   = h_MC.Integral(0,   h_MC.GetNbinsX() + 1)
            ratio = nData / nMC
            
            if self.verbose:
                print "{0} {1}: nData = {2:.3f}, nMC = {3:.3f}, ratio = {4:.3f}".format(era, region, nData, nMC, ratio)

            h_MC_normalized = h_MC.Clone("h_MC_normalized")
            h_MC_normalized.Scale(ratio)
        
            
            # draw histograms
            h_Data.Draw("hist")
            h_MC.Draw("hist same")
            # save histograms
            c.SetLogy(1) # set log y
            c.Update()
            c.SaveAs(plot_name + eraTag + ".pdf")
            c.SaveAs(plot_name + eraTag + ".png")
            
            # draw histograms
            h_Data.Draw("hist")
            h_MC_normalized.Draw("hist same")
            # save histograms
            c.SetLogy(1) # set log y
            c.Update()
            c.SaveAs(plot_name + "_normalized" + eraTag + ".pdf")
            c.SaveAs(plot_name + "_normalized" + eraTag + ".png")
        

def main():
    verbose = True
    S = Shape(verbose)
    S.getShape("condor/DataMC_2016_submission_2019-05-21_13-29-18/result.root", "2016")


if __name__ == "__main__":
    main()



