# shape_photon_met.py
import ROOT
import copy
import json
import os
import numpy as np
from tools import setupHist, getMETBinEdges, getSelections, removeCuts, stringifyMap, normalize, getNormalizedRatio, getTexSelection

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# title size
ROOT.gStyle.SetTitleSize(0.15, "t")
# how to set pad title size (for histogram):
# https://root-forum.cern.ch/t/setting-histogram-title-size-in-root/4468
# https://root.cern.ch/doc/master/classTStyle.html#a92b426badbae2e2d8dcceb38a8b93912

class Shape:
    def __init__(self, plot_dir, draw, doUnits, verbose):
        self.plot_dir               = plot_dir
        self.draw                   = draw
        self.doUnits                = doUnits
        self.verbose                = verbose
        self.splitQCD               = False
        self.systTag                = ""
        self.varTag                 = ""
        self.histos                 = {}
        self.cr_unit_histos         = {}
        self.cr_unit_histos_summed  = {}
        self.shape_map              = {}
        self.ratio_map              = {}
        self.ratio_rebinned_map     = {}
        self.eras = []
        # variable is also TDirectoryFile that holds histograms 
        self.variable    = "metWithPhoton"
        self.bin_types   = ["validation", "validationMetStudy", "search"]
        self.regions     = ["LowDM", "HighDM"]
        self.regions_tex = {
                             "LowDM"  : "Low $\Delta m$",
                             "HighDM" : "High $\Delta m$"
                           }
        self.bin_maps    = {}
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
        #self.label_met    = "#slash{E}_{T}^{#gamma} [GeV]"
        self.label_met    = "Modified p_{T}^{miss} [GeV]"
        self.label_events = "Events"
        self.label_ratio  = "Data / MC"
        self.labels = {
                        "met"   : "Modified p_{T}^{miss} [GeV]", 
                        "ht"    : "H_{T} [GeV]",
                        "nj"    : "N_{j}",
                        "nb"    : "N_{b}",
                        "nmt"   : "N_{merged tops}",
                        "nrt"   : "N_{resolved tops}",
                        "nw"    : "N_{W}"
                      }
        # colors
        self.color_black  = "black"
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        #self.color_list   = ["pinkish red", "tangerine", "emerald", "dark sky blue", "pinky purple", "reddish orange", "dark peach", "teal green", "blurple", "periwinkle blue", "apricot", "electric purple", "bright orange", "marigold", "ocean"]
        # rainbow ordered
        #self.color_list   = ["bubblegum pink", "hot pink", "red", "salmon", "red orange", "orange", "goldenrod", "light olive green", "medium green", "turquoise", "aqua blue", "bright blue", "dark sky blue", "blurple", "purple", "pastel purple"]
        # for shape study
        #self.color_list   = ["hot pink", "bright magenta", "red", "orange", "goldenrod", "light olive green", "aqua blue", "bright blue", "violet", "blurple", "purple", "pastel purple"]
        # colors matching data/MC plots
        self.color_list   = ["#0373E6", "#66CC33", "#FF9900", "#C843C8", "#6699FF", "#CC3333", "#339966", "#0533FF"]

    # here Tag variables should begin with underscore: e.g. _met, _2016, etc.
    def getSimpleMap(self, region, nameTag, dataSelectionTag, mcSelectionTag, variable, varTag = "", splitQCD=False):
        # only use varTag for MC, not data
        # only use varTag if variable name changes (e.g. nCRUnit, but not metWithPhoton)
        variableWithTag = variable + varTag
        if splitQCD:
            temp_map = {
                            "Data"              : "DataMC_Photon_" + region + nameTag + dataSelectionTag + 2 * variable        + "Datadata",
                            "GJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "#gamma+jetsstack",
                            "QCD_Direct"        : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "QCD Directstack",
                            "QCD_Fragmentation" : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "QCD Fragmentationstack",
                            "QCD_NonPrompt"     : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "QCD NonPromptstack",
                            "QCD_Fake"          : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "QCD Fakestack",
                            "WJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "W(l#nu)+jetsstack",
                            "TTG"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "t#bar{t}#gamma+jetsstack",
                            "tW"                : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "tWstack",
                            "Rare"              : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "Rarestack",
            }
        else:
            temp_map = {
                            "Data"              : "DataMC_Photon_" + region + nameTag + dataSelectionTag + 2 * variable + "Datadata",
                            "GJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "#gamma+jetsstack",
                            "QCD"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "QCDstack",
                            "WJets"             : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "W(l#nu)+jetsstack",
                            "TTG"               : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "t#bar{t}#gamma+jetsstack",
                            "tW"                : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "tWstack",
                            "Rare"              : "DataMC_Photon_" + region + nameTag + mcSelectionTag   + 2 * variableWithTag + "Rarestack",
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
        temp_map = {}
        for bin_type in self.bin_types:
            temp_map[bin_type] = {}
            for region in self.regions:
                temp_map[bin_type][region] = {}
                for selection in self.selections[bin_type][region]: 
                    # apply syst to MC only
                    dataSelectionTag = "_" + selection + "_jetpt30"
                    mcSelectionTag   = "_" + selection + self.systTag + "_jetpt30"
                    # do not use varTag for MET histogram systematics
                    temp_map[bin_type][region][selection] = self.getSimpleMap(region, nameTag, dataSelectionTag, mcSelectionTag, variable, splitQCD=self.splitQCD)

        return temp_map
    
    def getCRUnitMap(self, region, name, variable, era):
        # histogram examples
        # example; note that the variable is written twice
        # CR units: different variable, no selections
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30Datadata
        # DataMC_Photon_LowDM_nCRUnitLowDM_jetpt30_2016nCRUnitLowDM_drPhotonCleaned_jetpt30nCRUnitLowDM_drPhotonCleaned_jetpt30#gamma+jetsstack
        nameTag = "_" + name
        # apply syst to MC only
        dataSelectionTag = "_jetpt30"
        mcSelectionTag   = self.systTag + "_jetpt30"
        # use varTag for CR unit systematics
        temp_map = self.getSimpleMap(region, nameTag, dataSelectionTag, mcSelectionTag, variable, self.varTag, splitQCD=self.splitQCD)
        return temp_map
    
    def getShape(self, file_name, era, systTag = "", varTag=""): 
        self.systTag = systTag
        self.varTag = varTag
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
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Direct"]
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmentation"]
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_NonPrompt"]
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]
                        else:
                            print self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]
                    
                    # WARNING: strings loaded from json file have type 'unicode'
                    # ROOT cannot load histograms using unicode input: use type 'str'
                    h_Data              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"]            ) )
                    h_GJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["GJets"]           ) )
                    if self.splitQCD:
                        h_QCD_Direct        = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Direct"]          ) )
                        h_QCD_Fragmentation = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmentation"]   ) )
                        h_QCD_NonPrompt     = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_NonPrompt"]       ) )
                        h_QCD_Fake          = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"]            ) )
                    else:
                        h_QCD               = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"]             ) )
                    h_WJets             = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["WJets"]           ) )
                    h_TTG               = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["TTG"]             ) )
                    h_tW                = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["tW"]              ) )
                    h_Rare              = f.Get( str(self.variable + "/" + self.histos[era][bin_type][region][selection]["Rare"]            ) )
                    
                    # check if histograms load
                    if not h_Data:
                        print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["Data"])
                    if self.splitQCD:
                        if not h_QCD_Direct:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Direct"])
                        if not h_QCD_Fragmentation:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fragmentation"])
                        if not h_QCD_NonPrompt:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_NonPrompt"])
                        if not h_QCD_Fake:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD_Fake"])
                    else:
                        if not h_QCD:
                            print "ERROR: unable to load histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][region][selection]["QCD"])
                    
                    # MC_background
                    # combine all MC in denominator
                    h_mc = h_GJets.Clone("h_mc") 
                    if self.splitQCD:
                        h_mc.Add(h_QCD_Direct)
                        h_mc.Add(h_QCD_Fragmentation)
                        h_mc.Add(h_QCD_NonPrompt)
                        h_mc.Add(h_QCD_Fake)
                    else:
                        h_mc.Add(h_QCD)
                    h_mc.Add(h_WJets)
                    h_mc.Add(h_TTG)
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
                    photonDataMCNorm = float(nNum) / float(nDen)
                    
                    if self.verbose:
                        print "{0} {1} {2} {3}: nNum = {4:.3f}, nDen = {5:.3f}, photonDataMCNorm = {6:.3f}".format(era, bin_type, region, selection, nNum, nDen, photonDataMCNorm)

                    h_den_normalized = h_den.Clone("h_den_normalized")
                    h_den_normalized.Scale(photonDataMCNorm)
                        
                    self.shape_map[era][bin_type][region][selection]["total_data"]              = nNum
                    self.shape_map[era][bin_type][region][selection]["total_mc"]                = nDen
                    self.shape_map[era][bin_type][region][selection]["photon_data_mc_norm"]     = photonDataMCNorm
                    self.shape_map[era][bin_type][region][selection]["total_data_tex"]          = "${0:d}$".format(int(nNum))
                    self.shape_map[era][bin_type][region][selection]["total_mc_tex"]            = "${0:.3f}$".format(nDen)
                    self.shape_map[era][bin_type][region][selection]["photon_data_mc_norm_tex"] = "${0:.3f}$".format(photonDataMCNorm)
                    
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
                # for JEC systematic:
                # nCRUnitLowDM_drPhotonCleaned_jetpt30_jesTotalUp
                # nCRUnitLowDM_drPhotonCleaned_jetpt30_jesTotalDown
                # nCRUnitHighDM_drPhotonCleaned_jetpt30_jesTotalUp
                # nCRUnitHighDM_drPhotonCleaned_jetpt30_jesTotalDown
                
                variable = "nCRUnit" + region + "_drPhotonCleaned_jetpt30"
                variableWithTag = variable + self.varTag
                
                # WARNING: strings loaded from json file have type 'unicode'
                # ROOT cannot load histograms using unicode input: use type 'str'
                # Note: systematic (and thus varTag) should be apply only to MC, not to data
                
                h_Data                  = f.Get( str(variable + "/" + self.cr_unit_histos[era][region]["Data"]              ) )
                h_GJets                 = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["GJets"]             ) )
                if self.splitQCD:
                    h_QCD_Direct        = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["QCD_Direct"]         ) )
                    h_QCD_Fragmentation = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["QCD_Fragmentation"]  ) )
                    h_QCD_NonPrompt     = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["QCD_NonPrompt"]      ) )
                    h_QCD_Fake          = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["QCD_Fake"]           ) )
                else:
                    h_QCD               = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["QCD"]     ) )
                h_WJets                 = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["WJets"]   ) )
                h_TTG                   = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["TTG"]     ) )
                h_tW                    = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["tW"]      ) )
                h_Rare                  = f.Get( str(variableWithTag + "/" + self.cr_unit_histos[era][region]["Rare"]    ) )
                
                # MC_background
                if self.splitQCD:
                    h_back = h_QCD_Direct.Clone("h_back")
                    h_back.Add(h_QCD_Fragmentation)
                    h_back.Add(h_QCD_NonPrompt)
                    h_back.Add(h_QCD_Fake)
                else:
                    h_back = h_QCD.Clone("h_back")
                h_back.Add(h_WJets)
                h_back.Add(h_TTG)
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


    def writeLine(self, line):
        self.output_file.write(line + "\n")
    
    # Run 2 Q values in table for analysis note 
    def makeTable(self, output_name, makeDoc=False):
        # values for table are stored as follows:
        # self.shape_map[era][bin_type][region][selection]["total_data"]          = nNum
        # self.shape_map[era][bin_type][region][selection]["total_mc"]            = nDen
        # self.shape_map[era][bin_type][region][selection]["photon_data_mc_norm"] = photonDataMCNorm
        era      = "Run2"
        bin_type = "search"
        channelsForTable    = ["total_data_tex", "total_mc_tex", "photon_data_mc_norm_tex"]
        header              = "$\Nb$ & $\Nj$ & $N_{\\text{data}}$ & $N_{\\text{MC}}$ & $Q$ \\\\"
        caption  = "Summary of the different regions in which the photon normalization $Q$ is applied."
        caption += " The total number of data and MC events for each selection are shown, as well as the ratio $Q$."
        with open(output_name, "w+") as f:
            self.output_file = f
            
            if makeDoc:
                # begin document
                self.writeLine("\documentclass{article}")
                self.writeLine("\usepackage[utf8]{inputenc}")
                self.writeLine("\usepackage{geometry}")
                self.writeLine("\usepackage{longtable}")
                self.writeLine("\\usepackage{xspace}")
                self.writeLine("\\usepackage{amsmath}")
                self.writeLine("\\usepackage{graphicx}")
                self.writeLine("\\usepackage{cancel}")
                self.writeLine("\\usepackage{amsmath}")
                self.writeLine("\geometry{margin=1in}")
                self.writeLine("\\input{VariableNames.tex}")
                self.writeLine("\\begin{document}")
            
            # begin table
            self.writeLine("\\begin{table}")
            self.writeLine("\\begin{center}")
            self.writeLine("\\caption{%s}" % caption)
            self.writeLine("\\label{tab:Qregions}")
            self.writeLine("\\begin{tabular}{%s}" % ( "c" * (2 + len(channelsForTable)) ) )
            self.writeLine(header)
            self.writeLine("\\hline")
            self.writeLine("\\multicolumn{2}{c}{low \dm normalization regions} \\\\")
            self.writeLine("\\hline")
            self.writeLine("0         &   $\leq5$       & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["LowDM"]["NBeq0_NJle5"][channel]) for channel in channelsForTable)) )
            self.writeLine("0         &   $\geq6$       & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["LowDM"]["NBeq0_NJge6"][channel]) for channel in channelsForTable)) )
            self.writeLine("1         &   --            & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["LowDM"]["NBeq1"][channel]      ) for channel in channelsForTable)) )
            self.writeLine("$\geq2$   &   --            & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["LowDM"]["NBge2"][channel]      ) for channel in channelsForTable)) )
            self.writeLine("$\geq2$   &   $\geq7$       & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["LowDM"]["NBge2_NJge7"][channel]) for channel in channelsForTable)) )
            self.writeLine("\\hline")
            self.writeLine("\\multicolumn{2}{c}{high \dm normalization regions} \\\\")
            self.writeLine("\\hline")
            self.writeLine("1         &   --            & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["HighDM"]["NBeq1"][channel]      )    for channel in channelsForTable)) )
            self.writeLine("1         &   $\geq7$       & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["HighDM"]["NBeq1_NJge7"][channel])    for channel in channelsForTable)) )
            #self.writeLine("2         &   --            & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["HighDM"]["NBeq2"][channel]      )    for channel in channelsForTable)) )
            self.writeLine("$\geq2$   &   --            & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["HighDM"]["NBge2"][channel]      )    for channel in channelsForTable)) )
            self.writeLine("$\geq2$   &   $\geq7$       & %s \\\\" % (" & ".join("{0}".format(self.shape_map[era][bin_type]["HighDM"]["NBge2_NJge7"][channel])    for channel in channelsForTable)) )
            self.writeLine("\\hline")
            # end table
            self.writeLine("\\end{tabular}")
            self.writeLine("\\end{center}")
            self.writeLine("\\end{table}")
            
            if makeDoc:
                # end document
                self.writeLine("\\end{document}")
    
    def makeComparison(self, bin_type):
        draw_option = "hist error"
        
        ###################
        # Draw Histograms #
        ###################

        # draw histograms
        c = ROOT.TCanvas("c", "c", 800, 800)
        pad = c.cd(1)
        leftMargin      = 0.150
        rightMargin     = 0.075
        topMargin       = 0.180
        bottomMargin    = 0.150
        # set ticks on all sides of plot
        pad.SetTickx()
        pad.SetTicky()
        pad.SetLeftMargin(leftMargin)
        pad.SetRightMargin(rightMargin)
        pad.SetTopMargin(topMargin)
        pad.SetBottomMargin(bottomMargin)
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_width  = 0.3
        legend_height = 0.3
        # legend in left corner
        legend_x1 = 0.025 + leftMargin
        legend_x2 = 0.025 + leftMargin + legend_width
        legend_y1 = 0.975 - topMargin  - legend_height
        legend_y2 = 0.975 - topMargin
        # legend in right corner
        #legend_x1 = 1.0 - rightMargin - legend_width
        #legend_x2 = 1.0 - rightMargin
        #legend_y1 = 1.0 - topMargin   - legend_height
        #legend_y2 = 1.0 - topMargin

        # use 3 years and Run2 for eras
        eras = ["2016", "2017", "2018", "Run2"]

        for region in self.regions:
            for selection in self.selections[bin_type][region]:
                for rebin in self.ratio_rebinned_map["2016"][bin_type][region][selection]: 
                    region_tex          = self.regions_tex[region] 
                    region_root_tex     = region_tex.replace("\\", "#")
                    region_root_tex     = region_root_tex.replace("$", "")
                    selections_tex      = getTexSelection(selection)
                    selections_root_tex = selections_tex.replace("\\", "#")
                    selections_root_tex = selections_root_tex.replace("$", "")
                    #title   = "Shape for {0} bins, {1}, {2}, {3}".format(bin_type, region, selection, rebin)
                    title   = "Shapes for {0}, {1}".format(region_root_tex, selections_root_tex)
                    x_title = self.labels["met"]
                    y_title = "Shape #left(S_{#gamma}#right)"
                    y_min   = 0.0
                    y_max   = 5.0
                    
                    # legend: TLegend(x1,y1,x2,y2)
                    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                    legend.SetFillStyle(0)
                    legend.SetBorderSize(0)
                    legend.SetLineWidth(1)
                    legend.SetNColumns(1)
                    legend.SetTextFont(42)
                    
                    # map for histograms
                    hist_map = {}
                    
                    for i, era in enumerate(eras):
                        color = self.color_list[i]
                        if era == "Run2":
                            color = self.color_black
                        hist_map[era] = self.ratio_rebinned_map[era][bin_type][region][selection][rebin]
                        #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                        setupHist(hist_map[era],   title, x_title, y_title, color,     y_min, y_max)
                        hist_map[era].GetXaxis().SetTitleOffset(1.5)
                        hist_map[era].GetXaxis().SetNdivisions(5, 5, 0, True)
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
    
    # ----------------------------------------------------- #
    # - study shape of data and different MC in photon CR - #
    # ----------------------------------------------------- #
    def studyShapes(self, file_name, region, era, variable, nameTag, varName, xbins, n_bins):
        
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        f = ROOT.TFile(file_name)
        
        draw_option = "hist"
        data_style  = "E"
        
        ###################
        # Draw Histograms #
        ###################

        # draw histograms
        c = ROOT.TCanvas("c", "c", 800, 800)
        pad = c.cd(1)
        # set ticks on all sides of plot
        pad.SetTickx()
        pad.SetTicky()
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.10)
        pad.SetTopMargin(0.05)
        pad.SetBottomMargin(0.15)
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.45
        legend_x2 = 0.90
        legend_y1 = 0.60
        legend_y2 = 0.90
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetLineWidth(1)
        legend.SetNColumns(1)
        legend.SetTextFont(42)
        
        dataSelectionTag = ""
        mcSelectionTag   = ""
        # getSimpleMap(self, region, nameTag, dataSelectionTag, mcSelectionTag, variable, varTag = "", splitQCD=self.splitQCD):
        hNames = self.getSimpleMap(region, nameTag, dataSelectionTag, mcSelectionTag, variable, varTag = "", splitQCD=True)
        hMap = {}
        hMap["Data"]                = f.Get( str(variable + "/" + hNames["Data"]              ) )
        hMap["GJets"]               = f.Get( str(variable + "/" + hNames["GJets"]             ) )
        hMap["QCD_Direct"]          = f.Get( str(variable + "/" + hNames["QCD_Direct"]        ) )
        hMap["QCD_Fragmentation"]   = f.Get( str(variable + "/" + hNames["QCD_Fragmentation"] ) )
        hMap["QCD_NonPrompt"]       = f.Get( str(variable + "/" + hNames["QCD_NonPrompt"]     ) )
        hMap["QCD_Fake"]            = f.Get( str(variable + "/" + hNames["QCD_Fake"]          ) )
        hMap["WJets"]               = f.Get( str(variable + "/" + hNames["WJets"]             ) )
        hMap["TTG"]                 = f.Get( str(variable + "/" + hNames["TTG"]               ) )
        hMap["tW"]                  = f.Get( str(variable + "/" + hNames["tW"]                ) )
        hMap["Rare"]                = f.Get( str(variable + "/" + hNames["Rare"]              ) )
        # add QCD Direct to QCD Fragmentation
        hMap["QCD_Fragmentation"].Add(hMap["QCD_Direct"])
        
        # ---------------------- #
        # --- Compare shapes --- #
        # ---------------------- #

        # use list to define order
        # Data and all MC
        #hList = ["Data", "GJets", "QCD_Fragmentation", "QCD_NonPrompt", "QCD_Fake", "WJets", "TTG", "tW", "Rare"]
        # only Data, GJets, QCD
        # order to draw
        hListDraw   = ["GJets", "QCD_Fragmentation", "QCD_NonPrompt", "QCD_Fake", "Data"]
        # order for legend
        hListLegend = ["Data", "GJets", "QCD_Fragmentation", "QCD_NonPrompt", "QCD_Fake"]

        hNew = {} 
        
        # draw
        for i, key in enumerate(hListDraw):
            hOriginal = hMap[key]
            # rebin 
            h = hOriginal.Rebin(n_bins, "h_" + key + "_rebinned", xbins)
            # normalize each histogram so that we can compare shapes
            normalize(h)
            if key == "Data":
                color = self.color_black
            else:
                color = self.color_list[i]
            # turn off title
            #title   = "{0} in {1} for {2}".format(varName, region, era)
            title   = ""
            x_title = varName
            if varName.lower() in self.labels:
                x_title = self.labels[varName.lower()]
            y_title = "Events (normalized)"
            y_min   = 0.0
            y_max   = 1.0
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h, title, x_title, y_title, color, y_min, y_max)
            h.GetXaxis().SetTitleOffset(1.5)
            h.GetXaxis().SetNdivisions(5, 5, 0, True)
            # change style for data
            if key == "Data":
                h.SetMarkerStyle(ROOT.kFullCircle)
                h.SetMarkerSize(1.25)
            # draw
            if i == 0:
                if key == "Data":
                    h.Draw(data_style)
                else:
                    h.Draw(draw_option)
            else:
                if key == "Data":
                    h.Draw(data_style + " same")
                else:
                    h.Draw(draw_option + " same")
            # save new histogram to map
            hNew[key] = h
        
        # legend
        for i, key in enumerate(hListLegend):
            h = hNew[key]
            # change style for data
            if key == "Data":
                legend.AddEntry(h,  key,  "pe")
            else:
                legend.AddEntry(h,  key,  "l")
                    
                    
        legend.Draw()
        # save histograms
        plot_name = "{0}StudyShapes_{1}_{2}_{3}".format(self.plot_dir, region, variable, era)
        c.Update()
        c.SaveAs(plot_name + ".pdf")
        c.SaveAs(plot_name + ".png")
        
        # ----------------------------------- #
        # --- Show affect on shape factor --- #
        # ----------------------------------- #
        
        # clear canvas
        c.Clear()
        # split canvas to show ratios
        c.Divide(1, 2)
        
        h_num_original = hMap["Data"].Clone("h_num_original")
        h_mc           = hMap["GJets"].Clone("h_mc") 
        h_mc.Add(hMap["QCD_Fragmentation"])
        h_mc.Add(hMap["QCD_NonPrompt"])
        h_mc.Add(hMap["QCD_Fake"])
        h_mc.Add(hMap["WJets"])
        h_mc.Add(hMap["TTG"])
        h_mc.Add(hMap["tW"])
        h_mc.Add(hMap["Rare"])
        
        # vary QCD components up and down
        varyList = ["QCD_Fragmentation", "QCD_NonPrompt", "QCD_Fake"]
        for i, key in enumerate(varyList): 
            # get histos
            h_den_nominal_original   = h_mc.Clone("h_den_nominal_original")
            h_den_up_original        = h_mc.Clone("h_den_up_original")
            h_den_down_original      = h_mc.Clone("h_den_down_original")
            # 50% variation up and down
            h_den_up_original.Add(hMap[key],      0.5)
            h_den_down_original.Add(hMap[key],   -0.5)
            # WARNING: do not rebin ratios; first rebin, then get ratio
            # rebin 
            h_num         = h_num_original.Rebin(n_bins,            "h_num_rebinned",           xbins)
            h_den_nominal = h_den_nominal_original.Rebin(n_bins,    "h_den_nominal_rebinned",   xbins)
            h_den_up      = h_den_up_original.Rebin(n_bins,         "h_den_up_rebinned",        xbins)
            h_den_down    = h_den_down_original.Rebin(n_bins,       "h_den_down_rebinned",      xbins)
            # get shapes by normalizing denominator to numerator and taking the ratio 
            h_shape_nominal = getNormalizedRatio(h_num, h_den_nominal) 
            h_shape_up      = getNormalizedRatio(h_num, h_den_up) 
            h_shape_down    = getNormalizedRatio(h_num, h_den_down) 
            h_ratio_up      = h_shape_up.Clone("h_ratio_up")
            h_ratio_down    = h_shape_down.Clone("h_ratio_down")
            h_ratio_up.Divide(h_shape_nominal)
            h_ratio_down.Divide(h_shape_nominal)
            
            #title   = "{0} with {1} varied 50% in {2} for {3}".format(varName, key, region, era)
            title   = "{0} varied 50% in {1}".format(key, region)
            y_title = "Data/Sim."
            y_min   = 0.5
            y_max   = 1.5
            
            # turn off x-axis titles for upper plot
            # turn off main title for lower plot
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h_shape_nominal, title, "",      y_title,             self.color_black, y_min, y_max, True)
            setupHist(h_shape_up,      title, "",      y_title,             self.color_red,   y_min, y_max, True)
            setupHist(h_shape_down,    title, "",      y_title,             self.color_blue,  y_min, y_max, True)
            setupHist(h_ratio_up,      "",    x_title, "Varied/Nominal",    self.color_red,   0.94,  1.06,  True)
            setupHist(h_ratio_down,    "",    x_title, "Varied/Nominal",    self.color_blue,  0.94,  1.06,  True)
            
            # label and title formatting
            labelSize           = 0.08 
            titleSize           = 0.08
            titleOffsetXaxis    = 1.20
            titleOffsetYaxis    = 0.80
            
            # upper plot
            #h_shape_nominal.SetTitleSize(0.075, "t")
            h_shape_nominal.GetXaxis().SetLabelSize(0) # turn off x-axis labels for upper plot
            h_shape_nominal.GetYaxis().SetLabelSize(labelSize)
            h_shape_nominal.GetYaxis().SetTitleSize(titleSize)
            h_shape_nominal.GetYaxis().SetTitleOffset(titleOffsetYaxis)
            h_shape_nominal.GetYaxis().SetNdivisions(5, 5, 0, True)
            
            # lower plot
            h_ratio_up.GetXaxis().SetLabelSize(labelSize)
            h_ratio_up.GetXaxis().SetTitleSize(titleSize)
            h_ratio_up.GetXaxis().SetTitleOffset(titleOffsetXaxis)
            h_ratio_up.GetXaxis().SetNdivisions(5, 5, 0, True)
            h_ratio_up.GetYaxis().SetLabelSize(labelSize)
            h_ratio_up.GetYaxis().SetTitleSize(titleSize)
            h_ratio_up.GetYaxis().SetTitleOffset(titleOffsetYaxis)
            h_ratio_up.GetYaxis().SetNdivisions(3, 5, 0, True)
            
            # draw
            # note: use different variables for different legends, one legend for each plot

            # shapes
            pad = c.cd(1)
            # set ticks on all sides of plot
            pad.SetTickx()
            pad.SetTicky()
            pad.SetLeftMargin(0.15)
            pad.SetRightMargin(0.10)
            pad.SetTopMargin(0.15)
            pad.SetBottomMargin(0.01)
            
            # new legend for each plot
            # legend: TLegend(x1,y1,x2,y2)
            legend_x1 = 0.40
            legend_x2 = 0.90
            legend_y1 = 0.50
            legend_y2 = 0.85
            # legend: TLegend(x1,y1,x2,y2)
            legend1 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend1.SetFillStyle(0)
            legend1.SetBorderSize(0)
            legend1.SetLineWidth(1)
            legend1.SetNColumns(1)
            legend1.SetTextFont(42)
            legend1.AddEntry(h_shape_nominal,  "Shape (nominal)",                "l")
            legend1.AddEntry(h_shape_up,       "Shape ({0} up)".format(key),     "l")
            legend1.AddEntry(h_shape_down,     "Shape ({0} down)".format(key),   "l")
            h_shape_nominal.Draw(draw_option)
            h_shape_up.Draw(draw_option + " same")
            h_shape_down.Draw(draw_option + " same")
            legend1.Draw()
            # ratios
            pad = c.cd(2)
            pad.SetGridy()
            # set ticks on all sides of plot
            pad.SetTickx()
            pad.SetTicky()
            pad.SetLeftMargin(0.15)
            pad.SetRightMargin(0.10)
            pad.SetTopMargin(0.01)
            pad.SetBottomMargin(0.25)
            
            # skip legend in lower plot for now
            # legend: TLegend(x1,y1,x2,y2)
            # legend2 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            # legend2.AddEntry(h_ratio_up,       "({0} up) / nominal".format(key),     "l")
            # legend2.AddEntry(h_ratio_down,     "({0} down) / nominal".format(key),   "l")
            
            h_ratio_up.Draw(draw_option)
            h_ratio_down.Draw(draw_option + " same")
            # legend2.Draw()
            
            # save histograms
            plot_name = "{0}VaryShapes_{1}_{2}_{3}_{4}".format(self.plot_dir, region, variable, key, era)
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



