# systematics.py

from colors import getColorIndex
import json
import os
import numpy as np
import ROOT
import tools
from tools import setupHist, getNormalizedRatio

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


class Systematic:
    def __init__(self, plot_dir, normalization, shape):
        self.plot_dir   = plot_dir
        self.N          = normalization
        self.S          = shape
        self.regions    = ["LowDM", "HighDM"]
        self.particles  = ["Electron", "Muon"]

    def getZRatio(self, root_file, region, selection, era, variable):
        debug = True
        eraTag = "_" + era
        selectionTag = "_" + selection
        
        # histogram names example 
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLDatadata;1  metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLDYstack;1   metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLt#bar{t}stack;1 metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLSingle tstack;1 metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLRarestack;1 metWithLL
        
        h_map_norm = {}
        for particle in self.particles:
            h_map_norm[particle] = { 
                "Data"      : "DataMC_" + particle + "_" + region + "_met" + selectionTag + eraTag + 2 * variable + "Datadata",
                "DY"        : "DataMC_" + particle + "_" + region + "_met" + selectionTag + eraTag + 2 * variable + "DYstack",
                "TTbar"     : "DataMC_" + particle + "_" + region + "_met" + selectionTag + eraTag + 2 * variable + "t#bar{t}stack",
                "SingleT"   : "DataMC_" + particle + "_" + region + "_met" + selectionTag + eraTag + 2 * variable + "Single tstack",
                "Rare"      : "DataMC_" + particle + "_" + region + "_met" + selectionTag + eraTag + 2 * variable + "Rarestack"
            }
        
            if debug:
                print "# Loading {0} histograms".format(particle)
                print str(variable + "/" + h_map_norm[particle]["Data"]     )
                print str(variable + "/" + h_map_norm[particle]["DY"]       )
                print str(variable + "/" + h_map_norm[particle]["TTbar"]    )
                print str(variable + "/" + h_map_norm[particle]["SingleT"]  )
                print str(variable + "/" + h_map_norm[particle]["Rare"]     )



        #WARNING: strings loaded from json file have type 'unicode'
        # ROOT cannot load histograms using unicode input: use type 'str'
        h_Data_Electron     = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Data"]    ) )
        h_DY_Electron       = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["DY"]      ) )
        h_TTbar_Electron    = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["TTbar"]   ) )
        h_SingleT_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["SingleT"] ) )
        h_Rare_Electron     = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Rare"]    ) )
        h_Data_Muon         = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Data"]        ) )
        h_DY_Muon           = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["DY"]          ) )
        h_TTbar_Muon        = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["TTbar"]       ) )
        h_SingleT_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["SingleT"]     ) )
        h_Rare_Muon         = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Rare"]        ) )
        
        # combine all MC in denominator
        h_mc = h_DY_Electron.Clone("h_mc") 
        h_mc.Add(h_TTbar_Electron)
        h_mc.Add(h_SingleT_Electron)
        h_mc.Add(h_Rare_Electron)
        h_mc.Add(h_DY_Muon)
        h_mc.Add(h_TTbar_Muon)
        h_mc.Add(h_SingleT_Muon)
        h_mc.Add(h_Rare_Muon)
        # combine all Data in numerator
        h_Data = h_Data_Electron.Clone("h_Data") 
        h_Data.Add(h_Data_Muon)
        
        # numerator = Data
        h_num = h_Data.Clone("h_num") 
        # denominator = MC
        h_den = h_mc.Clone("h_den")
        # contstant binning for plots
        h_num.Rebin(4)
        h_den.Rebin(4)
        h_ratio_normalized = getNormalizedRatio(h_num, h_den)
        return h_ratio_normalized

    def getPhotonRatio(self, root_file, region, selection, era, variable):
        eraTag = "_" + era
        selectionTag = "_" + selection
        # getSimpleMap(self, region, nameTag, selectionTag, eraTag, variable)
        h_map_shape = self.S.getSimpleMap(region, "_met", selectionTag, eraTag, variable)
        #WARNING: strings loaded from json file have type 'unicode'
        # ROOT cannot load histograms using unicode input: use type 'str'
        h_Data              = root_file.Get( str(variable + "/" + h_map_shape["Data"]            ) )
        h_GJets             = root_file.Get( str(variable + "/" + h_map_shape["GJets"]           ) )
        if self.S.splitQCD:
            h_QCD_Fragmented    = root_file.Get( str(variable + "/" + h_map_shape["QCD_Fragmented"]  ) )
            h_QCD_Fake          = root_file.Get( str(variable + "/" + h_map_shape["QCD_Fake"]        ) )
        else:
            h_QCD               = root_file.Get( str(variable + "/" + h_map_shape["QCD"]             ) )
        h_WJets             = root_file.Get( str(variable + "/" + h_map_shape["WJets"]           ) )
        h_TTG               = root_file.Get( str(variable + "/" + h_map_shape["TTG"]             ) )
        h_TTbar             = root_file.Get( str(variable + "/" + h_map_shape["TTbar"]           ) )
        h_tW                = root_file.Get( str(variable + "/" + h_map_shape["tW"]              ) )
        h_Rare              = root_file.Get( str(variable + "/" + h_map_shape["Rare"]            ) )
        
        # combine all MC in denominator
        h_mc = h_GJets.Clone("h_mc") 
        if self.S.splitQCD:
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
        
        # contstant binning for plots
        h_num.Rebin(4)
        h_den.Rebin(4)
        h_ratio_normalized = getNormalizedRatio(h_num, h_den)
        return h_ratio_normalized

    def makeZvsPhoton(self, file_name, era):
        draw_option = "hist error"
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        f = ROOT.TFile(file_name)
        c = ROOT.TCanvas("c", "c", 800, 800)
        c.Divide(1, 2)
        metPhoton = "metWithPhoton"
        metLepton = "metWithLL"
        selection = "jetpt30"
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 
        for region in self.regions:
            # MET with Z
            # getZRatio(region, selection, era, variable)
            h_ratio_lepton = self.getZRatio(f, region, selection, era, metLepton)
            # MET with photon
            # getPhotonRatio(self, region, selection, era, variable)
            h_ratio_photon = self.getPhotonRatio(f, region, selection, era, metPhoton)
            
            # Double Ratio for MET: Z / photon
            h_ratio_ZoverPhoton = h_ratio_lepton.Clone("h_ratio_ZoverPhoton")
            h_ratio_ZoverPhoton.Divide(h_ratio_photon)    

            title = "Z vs. Photon, {0}, {1}".format(region, era)
            x_title = "MET (GeV)" 
            y_title = "Data / MC"
            y_min = 0.0
            y_max = 2.0
            
            # TODO: fix plotting range problems that are happening when rebin is used
            # rebin 
            #xbins = np.array([250, 350, 450, 550, 650, 1000])
            #n_bins = len(xbins) - 1
            #h_ratio_lepton          = h_ratio_lepton.Rebin(n_bins, "", xbins)
            #h_ratio_photon          = h_ratio_photon.Rebin(n_bins, "", xbins)
            #h_ratio_ZoverPhoton     = h_ratio_ZoverPhoton.Rebin(n_bins, "", xbins)
            
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h_ratio_lepton,       title, x_title, y_title,                "vermillion",      y_min, y_max)
            setupHist(h_ratio_photon,       title, x_title, y_title,                "electric blue",   y_min, y_max)
            setupHist(h_ratio_ZoverPhoton,  title, x_title, "(Z to LL) / Photon",   "black",           y_min, y_max)
        
            
            # histograms
            pad = c.cd(1)
            pad.SetGrid()
           
            # draw
            h_ratio_lepton.Draw(draw_option)
            h_ratio_photon.Draw(draw_option + " same")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_ratio_lepton,     "Z to LL",       "l")
            legend.AddEntry(h_ratio_photon,     "Photon",        "l")
            legend.Draw()
            
            # ratio
            pad = c.cd(2)
            pad.SetGrid()
            
            # draw
            h_ratio_ZoverPhoton.Draw(draw_option)
            
            # save histograms
            plot_name = "{0}ZvsPhoton_{1}_{2}".format(self.plot_dir, region, era)
            c.Update()
            c.SaveAs(plot_name + ".pdf")
            c.SaveAs(plot_name + ".png")



