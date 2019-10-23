# shape_photon_met.py
import ROOT
import copy
import json
import os
import numpy as np
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
        self.cr_unit_histos ={}
        self.cr_unit_histos_summed ={}
        self.shape_map = {}
        # variable is also TDirectoryFile that holds histograms 
        self.variable   = "metWithPhoton"
        self.bin_types  = ["validation", "search"]
        self.regions    = ["LowDM", "HighDM"]
        self.bin_maps = {}
        with open("validation_bins_v3.json", "r") as j:
            self.bin_maps["validation"] = stringifyMap(json.load(j))
        with open("search_bins_v4.json", "r") as j:
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

    # here Tag variables should begin with underscore: e.g. _met, _2016, etc.
    def getSimpleMap(self, region, nameTag, selectionTag, eraTag, variable):
        temp_map = {
                        "Data"              : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "Datadata",
                        "GJets"             : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "#gamma+jetsstack",
                        "QCD_Fragmented"    : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "QCD Fragmentedstack",
                        "QCD_Fake"          : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "QCD Fakestack",
                        "WJets"             : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "W(l#nu)+jetsstack",
                        "TTG"               : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "t#bar{t}#gamma+jetsstack",
                        "TTbar"             : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "t#bar{t}stack",
                        "tW"                : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "tWstack",
                        "Rare"              : "DataMC_Photon_" + region + nameTag + selectionTag + eraTag + 2 * variable + "Rarestack",
        }
        return temp_map

    def getHistoMap(self, name, variable, era):
        # histogram examples
        # example; note that the variable is written twice
        # with selection
        # DataMC_Photon_LowDM_met_NBeq0_NJge6_jetpt30_2016metWithPhotonmetWithPhotonDatadata
        # without selection
        # DataMC_Photon_LowDM_met_jetpt30_2016metWithPhotonmetWithPhotonDatadata
        nameTag = "_" + name
        eraTag = "_" + era
        temp_map = {}
        for bin_type in self.bin_types:
            temp_map[bin_type] = {}
            for region in self.regions:
                temp_map[bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    selectionTag = "_" + selection + "_jetpt30"
                    temp_map[bin_type][region][selection] = self.getSimpleMap(region, nameTag, selectionTag, eraTag, variable)

        return temp_map
    
    def getCRUnitMap(self, region, name, variable, era):
        # histogram examples
        # example; note that the variable is written twice
        # CR units: different variable, no selections
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30Datadata
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30#gamma+jetsstack
        nameTag = "_" + name
        eraTag = "_" + era
        selectionTag = "_jetpt30"
        temp_map = self.getSimpleMap(region, nameTag, selectionTag, eraTag, variable)
        return temp_map
    
    def getShape(self, file_name, era): 
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
        
        # setup histogram map (standard shape histograms)
        self.histos[era] = self.getHistoMap("met", self.variable, era)
        
        # CR units
        # nCRUnitLowDM_drPhotonCleaned_jetpt30
        # nCRUnitHighDM_drPhotonCleaned_jetpt30 
        self.cr_unit_histos[era] = {}
        self.cr_unit_histos[era]["LowDM"]  = self.getCRUnitMap("LowDM",  "nCRUnitLowDM",  "nCRUnitLowDM_drPhotonCleaned_jetpt30",  era)
        self.cr_unit_histos[era]["HighDM"] = self.getCRUnitMap("HighDM", "nCRUnitHighDM", "nCRUnitHighDM_drPhotonCleaned_jetpt30", era)
        
        self.shape_map[era] = {}
        
        # standard shape histograms
        for bin_type in self.bin_types:
            self.shape_map[era][bin_type] = {}
            for region in self.regions:
                self.shape_map[era][bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    self.shape_map[era][bin_type][region][selection] = {}
                    plot_name = self.plot_dir + self.variable + "_" + region
                    if self.verbose:
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"]
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]
                    
                    #WARNING: strings loaded from json file have type 'unicode'
                    # ROOT cannot load histograms using unicode input: use type 'str'
                    h_Data              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]            ) )
                    h_GJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["GJets"]           ) )
                    h_QCD_Fragmented    = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"]  ) )
                    h_QCD_Fake          = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]        ) )
                    h_WJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["WJets"]           ) )
                    h_TTG               = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTG"]             ) )
                    h_TTbar             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTbar"]           ) )
                    h_tW                = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["tW"]              ) )
                    h_Rare              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Rare"]            ) )
                    
                    # check if histograms load
                    if not h_Data:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"])
                    if not h_QCD_Fragmented:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"])
                    if not h_QCD_Fake:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"])
                    
                    # MC_background
                    h_back = h_QCD_Fragmented.Clone("h_back")
                    h_back.Add(h_QCD_Fake)
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
                    
                    #if self.verbose:
                    #    print "{0} {1}: nNum = {2:.3f}, nDen = {3:.3f}, ratio = {4:.3f}".format(era, region, nNum, nDen, ratio)

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
                         
                    met_dict  = getMETBinEdges(self.bin_maps[bin_type], selection)
                    #list of list of bin labels
                    met_names = met_dict["names"]
                    # list of list of bin edges
                    met_xbins = met_dict["xbins"]
                    
                    if len(met_names) != len(met_names):
                        print "ERROR: met_names and met_xbins do not have the same length"
                        return
                    
                    for i in xrange(len(met_names)):
                        # list of bin labels
                        names   = met_names[i]
                        # list of bin edges
                        xbins   = met_xbins[i]
                        if self.verbose:
                            print "Rebin shape factor: selection = {0}, rebin index = {1}, names = {2}, xbins = {3}".format(selection, i + 1, names, xbins)
                        if len(names) != len(xbins) - 1:
                            print "ERROR: length of names should be 1 less than length of xbins, but it is not"
                            #return
                        n_bins  = len(names)
                        h_num_rebinned              = h_num.Rebin(n_bins, "h_num_rebinned", xbins)
                        h_den_rebinned              = h_den.Rebin(n_bins, "h_den_rebinned", xbins)
                        h_den_rebinned_normalized   = h_den_normalized.Rebin(n_bins, "h_den_rebinned_normalized", xbins)
                        
                        # ratios
                        h_ratio_rebinned = h_num_rebinned.Clone("h_ratio_rebinned")
                        h_ratio_rebinned.Divide(h_den_rebinned)
                        h_ratio_rebinned_normalized = h_num_rebinned.Clone("h_ratio_rebinned_normalized")
                        h_ratio_rebinned_normalized.Divide(h_den_rebinned_normalized)
                        
                        # save shape factors to map
                        for j in xrange(n_bins):
                            name = names[j]
                            self.shape_map[era][bin_type][region][selection][name]            = h_ratio_rebinned_normalized.GetBinContent(j + 1)
                            self.shape_map[era][bin_type][region][selection][name + "_error"] = h_ratio_rebinned_normalized.GetBinError(j + 1)
                            #if self.verbose:
                            #    print "setting value in shape_map: {0} {1} {2} {3} {4}".format(era, bin_type, region, selection, name)

                        ###################
                        # Draw Histograms #
                        ###################
                        if self.draw:
                            selectionTag = "_{0}_rebin{1}".format(selection, i + 1)
                            # setup histograms
                            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                            setupHist(h_num_rebinned,               self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_red,   10.0 ** -1, 10.0 ** 6)
                            setupHist(h_den_rebinned,               self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_blue,  10.0 ** -1, 10.0 ** 6)
                            setupHist(h_den_rebinned_normalized,    self.variable + "_" + region + eraTag, self.label_met, "Events",  self.color_blue,  10.0 ** -1, 10.0 ** 6)
                            setupHist(h_ratio_rebinned,             self.variable + "_" + region + eraTag, self.label_met, "(Data - Back.)/Sig.",         self.color_black, 0.0, 3.0)
                            setupHist(h_ratio_rebinned_normalized,  self.variable + "_" + region + eraTag, self.label_met, "(Data - Back.)/(Norm. Sig.)", self.color_black, 0.0, 3.0)

                            # Data and MC
                            
                            # draw histograms
                            h_den_rebinned.Draw(draw_option)
                            h_num_rebinned.Draw(draw_option + "same")
                            # legend: TLegend(x1,y1,x2,y2)
                            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                            legend.AddEntry(h_num_rebinned, "Data - Background", "l")
                            legend.AddEntry(h_den_rebinned, "Signal",            "l")
                            legend.Draw()
                            # save histograms
                            c.SetLogy(1) # set log y
                            c.Update()
                            c.SaveAs(plot_name + selectionTag + eraTag + ".pdf")
                            c.SaveAs(plot_name + selectionTag + eraTag + ".png")
                            
                            # Data and Normalized MC
                            
                            # draw histograms
                            h_den_rebinned_normalized.Draw(draw_option)
                            h_num_rebinned.Draw(draw_option + "same")
                            # legend: TLegend(x1,y1,x2,y2)
                            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                            legend.AddEntry(h_num_rebinned,            "Data - Background", "l")
                            legend.AddEntry(h_den_rebinned_normalized, "Normalized Signal", "l")
                            legend.Draw()
                            # save histograms
                            c.SetLogy(1) # set log y
                            c.Update()
                            c.SaveAs(plot_name + "_normalized" + selectionTag + eraTag + ".pdf")
                            c.SaveAs(plot_name + "_normalized" + selectionTag + eraTag + ".png")
                            
                            # Data/MC

                            # draw histograms
                            h_ratio_rebinned.Draw(draw_option)
                            # legend: TLegend(x1,y1,x2,y2)
                            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                            legend.AddEntry(h_ratio_rebinned, "(Data - Back.)/Sig.", "l")
                            legend.Draw()
                            # save histograms
                            c.SetLogy(0) # unset log y
                            c.Update()
                            c.SaveAs(plot_name + "_ratio" + selectionTag + eraTag + ".pdf")
                            c.SaveAs(plot_name + "_ratio" + selectionTag + eraTag + ".png")
                            
                            # Data/(Normalized MC)
                            
                            # draw histograms
                            h_ratio_rebinned_normalized.Draw(draw_option)
                            # legend: TLegend(x1,y1,x2,y2)
                            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                            legend.AddEntry(h_ratio_rebinned_normalized, "(Data - Back.)/(Norm. Sig.)", "l")
                            legend.Draw()
                            # save histograms
                            c.SetLogy(0) # unset log y
                            c.Update()
                            c.SaveAs(plot_name + "_ratio_normalized" + selectionTag + eraTag + ".pdf")
                            c.SaveAs(plot_name + "_ratio_normalized" + selectionTag + eraTag + ".png")


        h_map = {}
        
        # CR unit shape histograms
        for region in self.regions:
            h_map[region] = {}

            # nCRUnitLowDM_drPhotonCleaned_jetpt30
            # nCRUnitHighDM_drPhotonCleaned_jetpt30 
            variable = "nCRUnit" + region + "_drPhotonCleaned_jetpt30"
            #WARNING: strings loaded from json file have type 'unicode'
            # ROOT cannot load histograms using unicode input: use type 'str'
            samples = ["Data", "GJets", "QCD_Fragmented", "QCD_Fake", "WJets", "TTG", "TTbar", "tW", "Rare"]
            print "Shape factor CR units; Loading {0} histograms".format(region)
            for sample in samples:
                hist_name = str(variable + "/" + self.cr_unit_histos[era][region][sample])
                print "\t{0}".format(hist_name) 
            
            h_Data              = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["Data"]              ) )
            h_GJets             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["GJets"]             ) )
            h_QCD_Fragmented    = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["QCD_Fragmented"]    ) )
            h_QCD_Fake          = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["QCD_Fake"]          ) )
            h_WJets             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["WJets"]             ) )
            h_TTG               = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["TTG"]               ) )
            h_TTbar             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["TTbar"]             ) )
            h_tW                = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["tW"]                ) )
            h_Rare              = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["Rare"]              ) )
            
            # MC_background
            h_back = h_QCD_Fragmented.Clone("h_back")
            h_back.Add(h_QCD_Fake)
            h_back.Add(h_WJets)
            h_back.Add(h_TTG)
            h_back.Add(h_TTbar)
            h_back.Add(h_tW)
            h_back.Add(h_Rare)

            # Data, MC_background, MC_gjets
            h_map[region]["data"]         = h_Data
            h_map[region]["mc_back"]      = h_back
            h_map[region]["mc_gjets"]     = h_GJets
        
        # WARNING
        # - histograms will be deleted when TFile is closed
        # - histograms need to be copied to use them later on 
        self.cr_unit_histos_summed[era] = copy.deepcopy(h_map)




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



