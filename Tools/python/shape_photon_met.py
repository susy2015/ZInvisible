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
    def __init__(self, plot_dir, draw, doUnits, verbose):
        self.plot_dir               = plot_dir
        self.draw                   = draw
        self.doUnits                = doUnits
        self.verbose                = verbose
        self.splitQCD               = False
        self.systTag                = ""
        self.histos                 = {}
        self.cr_unit_histos         = {}
        self.cr_unit_histos_summed  = {}
        self.shape_map              = {}
        self.ratio_map              = {}
        self.ratio_rebinned_map     = {}
        self.eras = []
        # variable is also TDirectoryFile that holds histograms 
        self.variable   = "metWithPhoton"
        self.bin_types  = ["validation", "validationMetStudy", "search"]
        self.regions    = ["LowDM", "HighDM"]
        self.bin_maps = {}
        with open("validation_bins_v3.json", "r") as j:
            self.bin_maps["validation"] = stringifyMap(json.load(j))
        with open("validation_bins_metStudy.json", "r") as j:
            self.bin_maps["validationMetStudy"] = stringifyMap(json.load(j))
        with open("search_bins_v4.json", "r") as j:
            self.bin_maps["search"] = stringifyMap(json.load(j))
        
        # Note: some selections are repeated, and there can be different MET binning for the same selection
        # get selections from json file
        self.selections = {}
        for bin_type in self.bin_types:
            self.selections[bin_type] = getSelections(self.bin_maps[bin_type], bin_type, ["NSV", "MET"])
        
        # labels
        self.label_met    = "#slash{E}_{T}^{#gamma} [GeV]"
        self.label_events = "Events"
        self.label_ratio  = "Data / MC"
        # colors
        self.color_black  = "black"
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_list   = ["pinkish red", "tangerine", "emerald", "dark sky blue", "pinky purple"]

    # here Tag variables should begin with underscore: e.g. _met, _2016, etc.
    def getSimpleMap(self, region, nameTag, dataSelectionTag, mcSelectionTag, eraTag, variable):
        # testing
        #print "In getSimpleMap():" + "DataMC_Photon_" + region + nameTag + dataSelectionTag + 2 * variable + "Datadata"
        #print "In getSimpleMap():" + "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "#gamma+jetsstack"
        if self.splitQCD:
            temp_map = {
                            "Data"              : "DataMC_Photon_" + region + nameTag + dataSelectionTag + 2 * variable + "Datadata",
                            "GJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "#gamma+jetsstack",
                            "QCD_Fragmented"    : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "QCD Fragmentedstack",
                            "QCD_Fake"          : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "QCD Fakestack",
                            "WJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "W(l#nu)+jetsstack",
                            "TTG"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "t#bar{t}#gamma+jetsstack",
                            "TTbar"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "t#bar{t}stack",
                            "tW"                : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "tWstack",
                            "Rare"              : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "Rarestack",
            }
        else:
            temp_map = {
                            "Data"              : "DataMC_Photon_" + region + nameTag + dataSelectionTag + 2 * variable + "Datadata",
                            "GJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "#gamma+jetsstack",
                            "QCD"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "QCDstack",
                            "WJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "W(l#nu)+jetsstack",
                            "TTG"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "t#bar{t}#gamma+jetsstack",
                            "TTbar"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "t#bar{t}stack",
                            "tW"                : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "tWstack",
                            "Rare"              : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variable + "Rarestack",
            }
        return temp_map

    def getHistoMap(self, name, variable, era):
        # histogram examples
        # example; note that the variable is written twice
        # with selection
        # DataMC_Photon_LowDM_met_NBeq0_NJge6_jetpt30_2016metWithPhotonmetWithPhotonDatadata
        # without selection
        # DataMC_Photon_LowDM_met_jetpt30_2016metWithPhotonmetWithPhotonDatadata
        # systematics
        # KEY: TH1D    DataMC_Photon_LowDM_met_NBge2_NJge7_prefire_syst_up_jetpt30metWithPhotonmetWithPhoton#gamma+jetsstack;1 metWithPhoton
        # KEY: TH1D    DataMC_Photon_HighDM_met_NBge2_NJge7_pileup_syst_up_jetpt30metWithPhotonmetWithPhoton#gamma+jetsstack;1 metWithPhoton
        nameTag = "_" + name
        eraTag = "_" + era
        temp_map = {}
        for bin_type in self.bin_types:
            temp_map[bin_type] = {}
            for region in self.regions:
                temp_map[bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    # apply syst to MC only
                    dataSelectionTag = "_" + selection + "_jetpt30"
                    mcSelectionTag   = "_" + selection + self.systTag + "_jetpt30"
                    temp_map[bin_type][region][selection] = self.getSimpleMap(region, nameTag, dataSelectionTag, mcSelectionTag, eraTag, variable)

        return temp_map
    
    def getCRUnitMap(self, region, name, variable, era):
        # histogram examples
        # example; note that the variable is written twice
        # CR units: different variable, no selections
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30Datadata
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30#gamma+jetsstack
        nameTag = "_" + name
        eraTag = "_" + era
        # apply syst to MC only
        dataSelectionTag = "_jetpt30"
        mcSelectionTag   = self.systTag + "_jetpt30"
        temp_map = self.getSimpleMap(region, nameTag, dataSelectionTag, mcSelectionTag, eraTag, variable)
        return temp_map
    
    def getShape(self, file_name, era, systTag = ""): 
        self.systTag = systTag
        self.eras.append(era)
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
        
        if self.doUnits:
            # CR units
            # nCRUnitLowDM_drPhotonCleaned_jetpt30
            # nCRUnitHighDM_drPhotonCleaned_jetpt30 
            self.cr_unit_histos[era] = {}
            self.cr_unit_histos[era]["LowDM"]  = self.getCRUnitMap("LowDM",  "nCRUnitLowDM",  "nCRUnitLowDM_drPhotonCleaned_jetpt30",  era)
            self.cr_unit_histos[era]["HighDM"] = self.getCRUnitMap("HighDM", "nCRUnitHighDM", "nCRUnitHighDM_drPhotonCleaned_jetpt30", era)
        
        self.shape_map[era] = {}
        ratio_map           = {}
        ratio_rebinned_map  = {}
        
        # standard shape histograms
        for bin_type in self.bin_types:
            self.shape_map[era][bin_type] = {}
            ratio_map[bin_type] = {}
            ratio_rebinned_map[bin_type] = {}
            for region in self.regions:
                self.shape_map[era][bin_type][region] = {}
                ratio_map[bin_type][region] = {}
                ratio_rebinned_map[bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    self.shape_map[era][bin_type][region][selection] = {}
                    ratio_rebinned_map[bin_type][region][selection] = {}
                    plot_name = self.plot_dir + self.variable + "_" + region
                    if self.verbose:
                        print self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]
                        if self.splitQCD:
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"]
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]
                        else:
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]
                    
                    #WARNING: strings loaded from json file have type 'unicode'
                    # ROOT cannot load histograms using unicode input: use type 'str'
                    h_Data              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]            ) )
                    h_GJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["GJets"]           ) )
                    if self.splitQCD:
                        h_QCD_Fragmented    = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"]  ) )
                        h_QCD_Fake          = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]        ) )
                    else:
                        h_QCD               = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]             ) )
                    h_WJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["WJets"]           ) )
                    h_TTG               = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTG"]             ) )
                    h_TTbar             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTbar"]           ) )
                    h_tW                = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["tW"]              ) )
                    h_Rare              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Rare"]            ) )
                    
                    # check if histograms load
                    if not h_Data:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"])
                    if self.splitQCD:
                        if not h_QCD_Fragmented:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmented"])
                        if not h_QCD_Fake:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"])
                    else:
                        if not h_QCD:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"])
                    
                    # MC_background
                    # combine all MC in denominator
                    h_mc = h_GJets.Clone("h_mc") 
                    if self.splitQCD:
                        h_mc.Add(h_QCD_Fragmented)
                        h_mc.Add(h_QCD_Fake)
                    else:
                        h_mc.Add(h_QCD)
                    h_mc.Add(h_WJets)
                    h_mc.Add(h_TTG)
                    h_mc.Add(h_TTbar)
                    h_mc.Add(h_tW)
                    h_mc.Add(h_Rare)

                    # numerator = Data
                    h_num = h_Data.Clone("h_num")
                     
                    # denominator = MC
                    h_den = h_mc.Clone("h_den") 
                    
                    ########################################
                    # Apply and Save Data/MC normalization #
                    ########################################
                    
                    # number of events for normalization
                    nNum  = h_num.Integral(0, h_num.GetNbinsX() + 1)
                    nDen  = h_den.Integral(0, h_den.GetNbinsX() + 1)
                    DataMCNorm = float(nNum) / float(nDen)
                    
                    if self.verbose:
                        print "{0} {1}: nNum = {2:.3f}, nDen = {3:.3f}, DataMCNorm = {4:.3f}".format(era, region, nNum, nDen, DataMCNorm)

                    print "{0} {1} {2} {3}: nNum = {4:.3f}, nDen = {5:.3f}, DataMCNorm = {6:.3f}".format(era, bin_type, region, selection, nNum, nDen, DataMCNorm)
                    h_den_normalized = h_den.Clone("h_den_normalized")
                    h_den_normalized.Scale(DataMCNorm)
                        
                    self.shape_map[era][bin_type][region][selection]["data_mc_norm"] = DataMCNorm
                    
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


                    ################################
                    # Save shape histograms to map #
                    ################################
                    ratio_map[bin_type][region][selection] = h_ratio_normalized
                         
                    met_dict  = getMETBinEdges(self.bin_maps[bin_type], selection)
                    #list of list of bin labels
                    met_names = met_dict["names"]
                    # list of list of bin edges
                    met_xbins = met_dict["xbins"]
                    
                    if len(met_names) != len(met_names):
                        print "ERROR: met_names and met_xbins do not have the same length"
                        return
                    
                    for i in xrange(len(met_names)):
                        rebinKey = "rebin{0}".format(i + 1)
                        # list of bin labels
                        names   = met_names[i]
                        # list of bin edges
                        xbins   = met_xbins[i]
                        if self.verbose:
                            print "Rebin shape factor: selection = {0}, rebin index = {1}, names = {2}, xbins = {3}".format(selection, i + 1, names, xbins)
                        #if len(names) != len(xbins) - 1:
                        #    print "ERROR: length of names should be 1 less than length of xbins, but it is not"
                        #    print "bin_type: {0}, selection: {1}, n_names: {2}, n_xbins: {3}, names: {4}, xbins:= {5}".format(bin_type, selection, len(names), len(xbins), names, xbins)
                        #    #return
                        n_bins  = len(names)
                        h_num_rebinned              = h_num.Rebin(n_bins, "h_num_rebinned", xbins)
                        h_den_rebinned              = h_den.Rebin(n_bins, "h_den_rebinned", xbins)
                        h_den_rebinned_normalized   = h_den_normalized.Rebin(n_bins, "h_den_rebinned_normalized", xbins)
                        
                        # ratios
                        h_ratio_rebinned = h_num_rebinned.Clone("h_ratio_rebinned")
                        h_ratio_rebinned.Divide(h_den_rebinned)
                        h_ratio_rebinned_normalized = h_num_rebinned.Clone("h_ratio_rebinned_normalized")
                        h_ratio_rebinned_normalized.Divide(h_den_rebinned_normalized)
                        
                        ################################
                        # Save shape histograms to map #
                        ################################
                        ratio_rebinned_map[bin_type][region][selection][rebinKey] = h_ratio_rebinned_normalized
                        
                        #############################
                        # Save shape factors to map #
                        #############################
                        
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
                            # put rebin number in name to distinguish different binnings
                            selectionTag = "_{0}_{1}".format(selection, rebinKey)
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
                            h_num_rebinned.Draw(draw_option + " same")
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
                            h_num_rebinned.Draw(draw_option + " same")
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

        # WARNING
        # - histograms will be deleted when TFile is closed
        # - histograms need to be copied to use them later on 
        self.ratio_map[era]             = copy.deepcopy(ratio_map)
        self.ratio_rebinned_map[era]    = copy.deepcopy(ratio_rebinned_map)

        # Unit bins
        
        if self.doUnits:
            h_map = {}
            # CR unit shape histograms
            for region in self.regions:
                h_map[region] = {}

                # nCRUnitLowDM_drPhotonCleaned_jetpt30
                # nCRUnitHighDM_drPhotonCleaned_jetpt30 
                variable = "nCRUnit" + region + "_drPhotonCleaned_jetpt30"
                #WARNING: strings loaded from json file have type 'unicode'
                # ROOT cannot load histograms using unicode input: use type 'str'
                if self.splitQCD:
                    samples = ["Data", "GJets", "QCD_Fragmented", "QCD_Fake", "WJets", "TTG", "TTbar", "tW", "Rare"]
                else:
                    samples = ["Data", "GJets", "QCD", "WJets", "TTG", "TTbar", "tW", "Rare"]
                #print "Shape factor CR units; Loading {0} histograms".format(region)
                for sample in samples:
                    hist_name = str(variable + "/" + self.cr_unit_histos[era][region][sample])
                    #print "\t{0}".format(hist_name) 
                
                h_Data              = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["Data"]              ) )
                h_GJets             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["GJets"]             ) )
                if self.splitQCD:
                    h_QCD_Fragmented    = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["QCD_Fragmented"]    ) )
                    h_QCD_Fake          = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["QCD_Fake"]          ) )
                else:
                    h_QCD               = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["QCD"]               ) )
                h_WJets             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["WJets"]             ) )
                h_TTG               = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["TTG"]               ) )
                h_TTbar             = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["TTbar"]             ) )
                h_tW                = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["tW"]                ) )
                h_Rare              = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["Rare"]              ) )
                
                # MC_background
                if self.splitQCD:
                    h_back = h_QCD_Fragmented.Clone("h_back")
                    h_back.Add(h_QCD_Fake)
                else:
                    h_back = h_QCD.Clone("h_back")
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


    def makeComparison(self, bin_type):
        draw_option = "hist error"
        
        ###################
        # Draw Histograms #
        ###################

        # draw histograms
        c = ROOT.TCanvas("c", "c", 800, 800)
        c.SetGrid()
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 

        # TODO: update to loop over eras or years... use list to accept any number of eras/years

        for region in self.regions:
            for selection in self.selections[bin_type][region]:
                for rebin in self.ratio_rebinned_map["2016"][bin_type][region][selection]: 
                    title = "Shape for {0} bins, {1}, {2}, {3}".format(bin_type, region, selection, rebin)
                    x_title = "MET (GeV)" 
                    y_title = "Shape #left(S_{#gamma}#right)"
                    y_min = -1.0
                    y_max = 3.0
                    
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    
                    # map for histograms
                    hist_map = {}
                    
                    for i, era in enumerate(self.eras):
                        color = self.color_list[i]
                        if era == "Run2":
                            color = self.color_black
                        hist_map[era] = self.ratio_rebinned_map[era][bin_type][region][selection][rebin]
                        #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                        setupHist(hist_map[era],   title, x_title, y_title, color,     y_min, y_max)
                        # draw
                        if i == 0:
                            hist_map[era].Draw(draw_option)
                        else:
                            hist_map[era].Draw(draw_option + " same")
                        legend.AddEntry(hist_map[era],     era,             "l")
                    
                    legend.Draw()
                    # save histograms
                    plot_name = "{0}Shape_{1}_{2}_{3}_{4}".format(self.plot_dir, bin_type, region, selection, rebin)
                    c.Update()
                    c.SaveAs(plot_name + ".pdf")
                    c.SaveAs(plot_name + ".png")
                    

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



