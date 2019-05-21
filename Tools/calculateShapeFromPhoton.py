# calculateShapeFromPhoton.py

import os
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Shape:
    def __init__(self, verbose):
        self.verbose = verbose
        self.histos = {}
        self.regions   = ["LowDM", "HighDM"]
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "metWithPhoton"

    def setupHistos(self, year):
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
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        # setup histogram map
        self.setupHistos(year)
        

def main():
    verbose = True
    S = Shape(verbose)
    S.getShape("quickResult_photonShape_2016.root", "2016")


if __name__ == "__main__":
    main()



