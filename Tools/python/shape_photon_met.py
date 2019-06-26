# shape_photon_met.py

import os
import numpy as np
import ROOT
from colors import getColorIndex
from tools import setupHist

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Shape:
    def __init__(self, plot_dir, draw, verbose):
        self.draw = draw
        self.verbose = verbose
        self.histos = {}
        self.regions   = ["LowDM", "HighDM"]
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "metWithPhoton"
        self.plot_dir = plot_dir
        if self.plot_dir[-1] != "/":
            self.plot_dir += "/"
        self.ratio_map = {}
        self.shape_map = {}
            
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

    def setupHistoMap(self, era):
        # histograms
        # example
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotonDatadata
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhoton#gamma+jetsstack
        # DataMC_Photon_LowDM_met_2016metWithPhotonmetWithPhotont#bar{t}#gamma+jetsstack
        self.histos[era] = {}
        for region in self.regions:
            self.histos[era][region] = {
                    "Data"  : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotonDatadata",
                    "GJets" : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhoton#gamma+jetsstack",
                    "QCD"   : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotonQCDstack",
                    "WJets" : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotonW(l#nu)+jetsstack",
                    "TTG"   : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotont#bar{t}#gamma+jetsstack",
                    "TTbar" : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotont#bar{t}stack",
                    "tW"    : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotontWstack",
                    "Rare"  : "DataMC_Photon_" + region + "_met_" + era + "metWithPhotonmetWithPhotonRarestack",
            }
    
    def getShape(self, file_name, era): 
        self.ratio_map[era] = {}
        self.shape_map[era] = {}
        draw_option = "hist error"
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

        xbin_250 = np.array([0.0, 250.0, 1000.0])
        xbin_300 = np.array([0.0, 250.0, 300.0, 1000.0])
        xbin_400 = np.array([0.0, 250.0, 400.0, 1000.0])
        
        # setup histogram map
        self.setupHistoMap(era)
        
        for region in self.regions:
            self.ratio_map[era][region] = {}
            self.shape_map[era][region] = {}
            
            plot_name = self.plot_dir + self.variable + "_" + region
            
            h_Data  = f.Get(self.variable + "/" + self.histos[era][region]["Data"])
            h_GJets = f.Get(self.variable + "/" + self.histos[era][region]["GJets"])
            h_QCD   = f.Get(self.variable + "/" + self.histos[era][region]["QCD"])
            h_WJets = f.Get(self.variable + "/" + self.histos[era][region]["WJets"])
            h_TTG   = f.Get(self.variable + "/" + self.histos[era][region]["TTG"])
            h_TTbar = f.Get(self.variable + "/" + self.histos[era][region]["TTbar"])
            h_tW    = f.Get(self.variable + "/" + self.histos[era][region]["tW"])
            h_Rare  = f.Get(self.variable + "/" + self.histos[era][region]["Rare"])
            
            # MC_background
            h_back = h_QCD.Clone("h_back")
            h_back.Add(h_WJets)
            h_back.Add(h_TTG)
            h_back.Add(h_TTbar)
            h_back.Add(h_tW)
            h_back.Add(h_Rare)

            # numerator = Data - MC_background
            h_num = h_Data.Clone("h_num")
            h_num.Add(h_back, -1)
             
            # denominator = MC_signal
            h_den = h_GJets.Clone("h_den") 
            
            # number of events for normalization
            nNum  = h_num.Integral(0, h_num.GetNbinsX() + 1)
            nDen  = h_den.Integral(0, h_den.GetNbinsX() + 1)
            ratio = float(nNum) / float(nDen)
            
            if self.verbose:
                print "{0} {1}: nNum = {2:.3f}, nDen = {3:.3f}, ratio = {4:.3f}".format(era, region, nNum, nDen, ratio)

            h_den_normalized = h_den.Clone("h_den_normalized")
            h_den_normalized.Scale(ratio)
            
            # rebin
            h_num.Rebin(2)
            h_den.Rebin(2)
            h_den_normalized.Rebin(2)
            h_num_250 = h_num.Rebin(2, "h_num_250", xbin_250)
            h_num_300 = h_num.Rebin(3, "h_num_300", xbin_300)
            h_num_400 = h_num.Rebin(3, "h_num_400", xbin_400)
            h_den_250 = h_den.Rebin(2, "h_den_250", xbin_250)
            h_den_300 = h_den.Rebin(3, "h_den_300", xbin_300)
            h_den_400 = h_den.Rebin(3, "h_den_400", xbin_400)
            h_den_normalized_250 = h_den_normalized.Rebin(2, "h_den_normalized_250", xbin_250)
            h_den_normalized_300 = h_den_normalized.Rebin(3, "h_den_normalized_300", xbin_300)
            h_den_normalized_400 = h_den_normalized.Rebin(3, "h_den_normalized_400", xbin_400)
            h_map = {
                        "standard" : {"num":h_num,     "den":h_den,     "den_norm":h_den_normalized},
                        "250"      : {"num":h_num_250, "den":h_den_250, "den_norm":h_den_normalized_250},
                        "300"      : {"num":h_num_300, "den":h_den_300, "den_norm":h_den_normalized_300},
                        "400"      : {"num":h_num_400, "den":h_den_400, "den_norm":h_den_normalized_400}
            }
            
            for key in h_map:
                keyTag = "_" + key
                h_num            = h_map[key]["num"]
                h_den            = h_map[key]["den"]
                h_den_normalized = h_map[key]["den_norm"]
            
                # ratios
                h_ratio = h_num.Clone("h_ratio")
                h_ratio.Divide(h_den)
                h_ratio_normalized = h_num.Clone("h_ratio_normalized")
                h_ratio_normalized.Divide(h_den_normalized)
        
                # setup histograms
                #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                setupHist(h_num,               self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_red,   10.0 ** -1, 10.0 ** 6)
                setupHist(h_den,               self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_blue,  10.0 ** -1, 10.0 ** 6)
                setupHist(h_den_normalized,    self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_blue,  10.0 ** -1, 10.0 ** 6)
                setupHist(h_ratio,             self.variable + "_" + region + eraTag, self.label_met, "(Data - Back.)/Sig.",         self.color_black, 0.0, 3.0)
                setupHist(h_ratio_normalized,  self.variable + "_" + region + eraTag, self.label_met, "(Data - Back.)/(Norm. Sig.)", self.color_black, 0.0, 3.0)
         
                # map for normalized ratios
                self.ratio_map[era][region][key] = h_ratio_normalized
                
                ###################
                # Draw Histograms #
                ###################
                if self.draw:

                    # Data and MC
                    
                    # draw histograms
                    h_den.Draw(draw_option)
                    h_num.Draw(draw_option + "same")
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    legend.AddEntry(h_num, "Data - Background", "l")
                    legend.AddEntry(h_den, "Signal",            "l")
                    legend.Draw()
                    # save histograms
                    c.SetLogy(1) # set log y
                    c.Update()
                    c.SaveAs(plot_name + keyTag + eraTag + ".pdf")
                    c.SaveAs(plot_name + keyTag + eraTag + ".png")
                    
                    # Data and Normalized MC
                    
                    # draw histograms
                    h_den_normalized.Draw(draw_option)
                    h_num.Draw(draw_option + "same")
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    legend.AddEntry(h_num,            "Data - Background", "l")
                    legend.AddEntry(h_den_normalized, "Normalized Signal", "l")
                    legend.Draw()
                    # save histograms
                    c.SetLogy(1) # set log y
                    c.Update()
                    c.SaveAs(plot_name + "_normalized" + keyTag + eraTag + ".pdf")
                    c.SaveAs(plot_name + "_normalized" + keyTag + eraTag + ".png")
                    
                    # Data/MC

                    # draw histograms
                    h_ratio.Draw(draw_option)
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    legend.AddEntry(h_ratio, "(Data - Back.)/Sig.", "l")
                    legend.Draw()
                    # save histograms
                    c.SetLogy(0) # unset log y
                    c.Update()
                    c.SaveAs(plot_name + "_ratio" + keyTag + eraTag + ".pdf")
                    c.SaveAs(plot_name + "_ratio" + keyTag + eraTag + ".png")
                    
                    # Data/(Normalized MC)
                    
                    # draw histograms
                    h_ratio_normalized.Draw(draw_option)
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    legend.AddEntry(h_ratio_normalized, "(Data - Back.)/(Norm. Sig.)", "l")
                    legend.Draw()
                    # save histograms
                    c.SetLogy(0) # unset log y
                    c.Update()
                    c.SaveAs(plot_name + "_ratio_normalized" + keyTag + eraTag + ".pdf")
                    c.SaveAs(plot_name + "_ratio_normalized" + keyTag + eraTag + ".png")
            
        # values for search bins
        # TODO: save shape values for search bins

        # values for validation  bins
        # LowDM
        self.shape_map[era]["LowDM"]["met_250to300"]       = self.ratio_map[era]["LowDM"]["300"].GetBinContent(2)
        self.shape_map[era]["LowDM"]["met_250to300_error"] = self.ratio_map[era]["LowDM"]["300"].GetBinError(2)
        self.shape_map[era]["LowDM"]["met_300toINF"]       = self.ratio_map[era]["LowDM"]["300"].GetBinContent(3)
        self.shape_map[era]["LowDM"]["met_300toINF_error"] = self.ratio_map[era]["LowDM"]["300"].GetBinError(3)
        self.shape_map[era]["LowDM"]["met_250to400"]       = self.ratio_map[era]["LowDM"]["400"].GetBinContent(2)
        self.shape_map[era]["LowDM"]["met_250to400_error"] = self.ratio_map[era]["LowDM"]["400"].GetBinError(2)
        self.shape_map[era]["LowDM"]["met_400toINF"]       = self.ratio_map[era]["LowDM"]["400"].GetBinContent(3)
        self.shape_map[era]["LowDM"]["met_400toINF_error"] = self.ratio_map[era]["LowDM"]["400"].GetBinError(3)
        self.shape_map[era]["LowDM"]["met_250toINF"]       = self.ratio_map[era]["LowDM"]["250"].GetBinContent(2)
        self.shape_map[era]["LowDM"]["met_250toINF_error"] = self.ratio_map[era]["LowDM"]["250"].GetBinError(2)
        # HighDM
        self.shape_map[era]["HighDM"]["met_250to400"]       = self.ratio_map[era]["HighDM"]["400"].GetBinContent(2)
        self.shape_map[era]["HighDM"]["met_250to400_error"] = self.ratio_map[era]["HighDM"]["400"].GetBinError(2)
        self.shape_map[era]["HighDM"]["met_400toINF"]       = self.ratio_map[era]["HighDM"]["400"].GetBinContent(3)
        self.shape_map[era]["HighDM"]["met_400toINF_error"] = self.ratio_map[era]["HighDM"]["400"].GetBinError(3)
        

def main():
    plot_dir = "more_plots"
    draw = False
    verbose = False
    S = Shape(plot_dir, draw, verbose)
    S.getShape("condor/DataMC_2016_submission_2019-05-21_13-29-18/result.root", "2016")
    S.getShape("condor/DataMC_2017_submission_2019-05-21_13-32-05/result.root", "2017")
    S.getShape("condor/DataMC_2018_AB_submission_2019-05-21_13-32-36/result.root", "2018_AB")
    S.getShape("condor/DataMC_2018_CD_submission_2019-05-21_13-33-23/result.root", "2018_CD")


if __name__ == "__main__":
    main()



