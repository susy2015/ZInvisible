# shape_photon_met.py
import json
import os
import numpy as np
import ROOT
from colors import getColorIndex
from tools import setupHist, getMETBinEdges, getSelections, removeCuts, stringifyMap

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Shape:
    def __init__(self, plot_dir, draw, verbose):
        self.draw = draw
        self.verbose = verbose
        self.plot_dir = plot_dir
        self.histos = {}
        self.ratio_map = {}
        self.shape_map = {}
        # variable is also TDirectoryFile that holds histograms 
        self.variable   = "metWithPhoton"
        self.bin_types  = ["validation", "search"]
        self.regions    = ["LowDM", "HighDM"]
        self.bin_maps = {}
        with open("validation_bins.json", "r") as j:
            self.bin_maps["validation"] = stringifyMap(json.load(j))
        with open("search_bins.json", "r") as j:
            self.bin_maps["search"] = stringifyMap(json.load(j))
        
        # Note: some selections are repeated, and there can be different MET binning for the same selection
        # get selections from json file
        self.selections = {}
        for bin_type in self.bin_types:
            self.selections[bin_type] = getSelections(self.bin_maps[bin_type], bin_type, "NSV")
        
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
        # histogram examples
        # example; note that the variable is written twice
        # with selection
        # DataMC_Photon_LowDM_met_NBeq0_NJge6_jetpt20_2016metWithPhotonmetWithPhotonDatadata
        # without selection
        # DataMC_Photon_LowDM_met_jetpt20_2016metWithPhotonmetWithPhotonDatadata
        eraTag = "_" + era
        self.histos[era] = {}
        for bin_type in self.bin_types:
            self.histos[era][bin_type] = {}
            for region in self.regions:
                self.histos[era][bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    selectionTag = "_" + selection + "_jetpt20"
                    self.histos[era][bin_type][region][selection] = {
                            "Data"  : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "Datadata",
                            "GJets" : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "#gamma+jetsstack",
                            "QCD"   : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "QCDstack",
                            "WJets" : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "W(l#nu)+jetsstack",
                            "TTG"   : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "t#bar{t}#gamma+jetsstack",
                            "TTbar" : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "t#bar{t}stack",
                            "tW"    : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "tWstack",
                            "Rare"  : "DataMC_Photon_" + region + "_met" + selectionTag + eraTag + 2 * self.variable + "Rarestack",
                    }
    
    def getShape(self, file_name, era): 
        skip = True
        draw_option = "hist error"
        eraTag = "_" + era
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        if self.verbose:
            print "file_name: {0}".format(file_name)
        f = ROOT.TFile(file_name)
        c = ROOT.TCanvas("c", "c", 800, 800)

        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 
        
        # setup histogram map
        self.setupHistoMap(era)
        
        self.ratio_map[era] = {}
        self.shape_map[era] = {}
        
        for bin_type in self.bin_types:
            self.ratio_map[era][bin_type] = {}
            self.shape_map[era][bin_type] = {}
            for region in self.regions:
                self.ratio_map[era][bin_type][region] = {}
                self.shape_map[era][bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    self.shape_map[era][bin_type][region][selection] = {}
                    plot_name = self.plot_dir + self.variable + "_" + region
                    if self.verbose:
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]
                    
                    #WARNING: strings loaded from json file have type 'unicode'
                    # ROOT cannot load histograms using unicode input: use type 'str'
                    h_Data  = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]  ) )
                    h_GJets = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["GJets"] ) )
                    h_QCD   = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]   ) )
                    h_WJets = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["WJets"] ) )
                    h_TTG   = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTG"]   ) )
                    h_TTbar = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTbar"] ) )
                    h_tW    = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["tW"]    ) )
                    h_Rare  = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Rare"]  ) )
                    
                    # check if histograms load
                    if not h_Data:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"])
                    if not h_QCD:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"])
                    
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
                    
                    # rebin in MET
                    # variable binning for background prediction
                    # Note: some selections are repeated, and there can be different MET binning for the same selection
                    # Loop over different met binning for given selection
                    # even if selection and MET binnings are repeated, it is ok overwrite them (as only selection and met binning determine value)
                    
                    # contstant binning for plots
                    h_num.Rebin(2)
                    h_den.Rebin(2)
                    h_den_normalized.Rebin(2)
                    
                    # ratios
                    h_ratio = h_num.Clone("h_ratio")
                    h_ratio.Divide(h_den)
                    h_ratio_normalized = h_num.Clone("h_ratio_normalized")
                    h_ratio_normalized.Divide(h_den_normalized)
                         
                    # map for normalized ratios
                    self.ratio_map[era][bin_type][region][selection] = h_ratio_normalized
                     
                    met_dict  = getMETBinEdges(self.bin_maps[bin_type], selection)
                    met_names = met_dict["names"]
                    met_xbins = met_dict["xbins"]
                    if len(met_names) != len(met_names):
                        print "ERROR: met_names and met_xbins do not have the same length"
                        return
                    
                    for i in xrange(len(met_names)):
                        names   = met_names[i]
                        xbins   = met_xbins[i]
                        if len(names) != len(xbins) - 1:
                            print "ERROR: length of names should be 1 less than length of xbins, but it is not"
                            return
                        n_bins  = len(names)
                        h_num_rebinned              = h_num.Rebin(n_bins, "h_num_rebinned", xbins)
                        h_den_rebinned              = h_den.Rebin(n_bins, "h_den_rebinned", xbins)
                        h_den_rebinned_normalized   = h_den_normalized.Rebin(n_bins, "h_den_rebinned_normalized", xbins)
                        
                        # ratios
                        h_ratio_rebinned = h_num_rebinned.Clone("h_ratio_rebinned")
                        h_ratio_rebinned.Divide(h_den_rebinned)
                        h_ratio_rebinned_normalized = h_num_rebinned.Clone("h_ratio_rebinned_normalized")
                        h_ratio_rebinned_normalized.Divide(h_den_rebinned_normalized)
                        
                        for j in xrange(n_bins):
                            name = names[j]
                            self.shape_map[era][bin_type][region][selection][name]            = h_ratio_rebinned_normalized.GetBinContent(j + 1)
                            self.shape_map[era][bin_type][region][selection][name + "_error"] = h_ratio_rebinned_normalized.GetBinError(j + 1)
                            if self.verbose:
                                print "setting value in shape_map: {0} {1} {2} {3} {4}".format(era, bin_type, region, selection, name)

                    # TODO: update "skip" section
                    if not skip:
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
                            self.ratio_map[era][bin_type][region][selection][key] = h_ratio_normalized
                            
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

def main():
    json_file = "runs/run_2019-07-17.json"
    eras = ["2016", "2017", "2018_AB", "2018_CD"]
    #eras = ["2016", "2017", "2018_PreHEM", "2018_PostHEM"]
    plot_dir = "more_plots"
    if plot_dir[-1] != "/":
        plot_dir += "/"
    draw = True
    verbose = False
    S = Shape(plot_dir, draw, verbose)
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            S.getShape(result_file, era)


if __name__ == "__main__":
    main()



