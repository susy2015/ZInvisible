# calculateShapeFromPhoton.py

import os
import numpy as np
import ROOT
from colors import getColorIndex

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Shape:
    def __init__(self, plot_dir, verbose):
        self.verbose = verbose
        self.histos = {}
        self.regions   = ["LowDM", "HighDM"]
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "metWithPhoton"
        self.plot_dir = plot_dir
        if self.plot_dir[-1] != "/":
            self.plot_dir += "/"
            
        # labels
        self.label_met    = "#slash{E}_{T}^{#gamma} [GeV]"
        self.label_events = "Events"
        self.label_ratio  = "Data / MC"
        # colors
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_black  = "black"

    def setupHist(self, hist, title, x_title, y_title, color, y_min, y_max):
        x_axis = hist.GetXaxis()
        y_axis = hist.GetYaxis()
        y_axis.SetRangeUser(y_min, y_max)
        hist.SetTitle(title)
        x_axis.SetTitle(x_title)
        y_axis.SetTitle(y_title)
        hist.SetStats(ROOT.kFALSE)
        hist.SetLineColor(getColorIndex(color))
        hist.SetLineWidth(3)

    def setupHistoMap(self, year):
        # histograms
        # example
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonDatadata
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack
        self.histos[year] = {}
        for region in self.regions:
            self.histos[year][region] = {
                    "Data"  : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotonDatadata",
                    "GJets" : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhoton#gamma+jetsstack",
                    "QCD"   : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotonQCDstack",
                    "WJets" : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotonW(l#nu)+jetsstack",
                    "TTbar" : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotont#bar{t}stack",
                    "tW"    : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotontWstack",
                    "Rare"  : "DataMC_Photon_" + region + "_met_" + year + "metWithPhotonmetWithPhotonRarestack",
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
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 
        
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
            
            # add MC
            h_MC = h_GJets.Clone("h_MC") 
            h_MC.Add(h_QCD)
            h_MC.Add(h_WJets)
            h_MC.Add(h_TTbar)
            h_MC.Add(h_tW)
            h_MC.Add(h_Rare)
            
            # number of events
            nData = h_Data.Integral(0, h_Data.GetNbinsX() + 1)
            nMC   = h_MC.Integral(0,   h_MC.GetNbinsX() + 1)
            ratio = nData / nMC
            
            if self.verbose:
                print "{0} {1}: nData = {2:.3f}, nMC = {3:.3f}, ratio = {4:.3f}".format(era, region, nData, nMC, ratio)

            h_MC_normalized = h_MC.Clone("h_MC_normalized")
            h_MC_normalized.Scale(ratio)
            
            # ratios
            h_ratio = h_Data.Clone("h_ratio")
            h_ratio.Divide(h_MC)
            h_ratio_normalized = h_Data.Clone("h_ratio_normalized")
            h_ratio_normalized.Divide(h_MC_normalized)
        
            # setup histograms
            #setupHist(self, hist, title, x_title, y_title, color, y_min, y_max)
            self.setupHist(h_Data,              self.variable + "_" + region + eraTag, self.label_met, "Events", self.color_red,   10.0 ** -1, 10.0 ** 6)
            self.setupHist(h_MC,                self.variable + "_" + region + eraTag, self.label_met, "Events", self.color_blue,  10.0 ** -1, 10.0 ** 6)
            self.setupHist(h_MC_normalized,     self.variable + "_" + region + eraTag, self.label_met, "Events", self.color_blue,  10.0 ** -1, 10.0 ** 6)
            self.setupHist(h_ratio,             self.variable + "_" + region + eraTag, self.label_met, "Data/MC", self.color_black, 0.0, 3.0)
            self.setupHist(h_ratio_normalized,  self.variable + "_" + region + eraTag, self.label_met, "Data/MC", self.color_black, 0.0, 3.0)
         
            # Data and MC
            
            # draw histograms
            h_Data.Draw("hist")
            h_MC.Draw("hist same")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_Data, "Data", "l")
            legend.AddEntry(h_MC,   "MC",   "l")
            legend.Draw()
            # save histograms
            c.SetLogy(1) # set log y
            c.Update()
            c.SaveAs(plot_name + eraTag + ".pdf")
            c.SaveAs(plot_name + eraTag + ".png")
            
            # Data and Normalized MC
            
            # draw histograms
            h_Data.Draw("hist")
            h_MC_normalized.Draw("hist same")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_Data,          "Data",          "l")
            legend.AddEntry(h_MC_normalized, "Normalized MC", "l")
            legend.Draw()
            # save histograms
            c.SetLogy(1) # set log y
            c.Update()
            c.SaveAs(plot_name + "_normalized" + eraTag + ".pdf")
            c.SaveAs(plot_name + "_normalized" + eraTag + ".png")
            
            # Data/MC

            # draw histograms
            h_ratio.Draw("hist")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_ratio, "Data/MC", "l")
            legend.Draw()
            # save histograms
            c.SetLogy(0) # unset log y
            c.Update()
            c.SaveAs(plot_name + "_ratio" + eraTag + ".pdf")
            c.SaveAs(plot_name + "_ratio" + eraTag + ".png")
            
            # Data/(Normalized MC)
            
            # draw histograms
            h_ratio_normalized.Draw("hist")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_ratio_normalized, "Data/(Normalized MC)", "l")
            legend.Draw()
            # save histograms
            c.SetLogy(0) # unset log y
            c.Update()
            c.SaveAs(plot_name + "_ratio_normalized" + eraTag + ".pdf")
            c.SaveAs(plot_name + "_ratio_normalized" + eraTag + ".png")
        

def main():
    plot_dir = "more_plots"
    verbose = True
    S = Shape(plot_dir, verbose)
    S.getShape("condor/DataMC_2016_submission_2019-05-21_13-29-18/result.root", "2016")
    S.getShape("condor/DataMC_2017_submission_2019-05-21_13-32-05/result.root", "2017")
    S.getShape("condor/DataMC_2018_AB_submission_2019-05-21_13-32-36/result.root", "2018_AB")
    S.getShape("condor/DataMC_2018_CD_submission_2019-05-21_13-33-23/result.root", "2018_CD")


if __name__ == "__main__":
    main()



