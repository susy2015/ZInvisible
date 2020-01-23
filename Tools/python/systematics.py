# systematics.py

import copy
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
        # WARNING: for rebinning, xbins need to be a numpy array of floats
        # rebinned
        self.xbins   = np.array([0.0, 250.0, 350.0, 450.0, 550.0, 650.0, 1000.0])
        self.n_bins  = len(self.xbins) - 1
        self.met_min = 250.0
        self.met_max = 1000.0
        self.h_map_syst = {}

    def getZRatio(self, root_file, region, selection, era, variable, rebin):
        debug = False
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
                "Data"      : "DataMC_" + particle + "_" + region + "_met" + selectionTag + 2 * variable + "Datadata",
                "DY"        : "DataMC_" + particle + "_" + region + "_met" + selectionTag + 2 * variable + "DYstack",
                "TTbar"     : "DataMC_" + particle + "_" + region + "_met" + selectionTag + 2 * variable + "t#bar{t}stack",
                "SingleT"   : "DataMC_" + particle + "_" + region + "_met" + selectionTag + 2 * variable + "Single tstack",
                "Rare"      : "DataMC_" + particle + "_" + region + "_met" + selectionTag + 2 * variable + "Rarestack"
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
        
        # WARNING: do not rebin ratios; first rebin, then get ratio
        # numerator = Data
        # denominator = MC
        if rebin:
            # variable binning
            h_num = h_Data.Rebin(self.n_bins, "h_num", self.xbins)
            h_den = h_mc.Rebin(self.n_bins, "h_den", self.xbins)
        else:
            h_num = h_Data.Clone("h_num") 
            h_den = h_mc.Clone("h_den")
            # contstant binning for plots
            h_num.Rebin(2)
            h_den.Rebin(2)
        h_ratio_normalized = getNormalizedRatio(h_num, h_den)
        return h_ratio_normalized

    def getPhotonRatio(self, root_file, region, selection, era, variable, rebin):
        eraTag = "_" + era
        selectionTag = "_" + selection
        # getSimpleMap(self, region, nameTag, dataSelectionTag, mcSelectionTag, eraTag, variable):
        h_map_shape = self.S.getSimpleMap(region, "_met", selectionTag, selectionTag, eraTag, variable)
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

        # WARNING: do not rebin ratios; first rebin, then get ratio
        # numerator = Data
        # denominator = MC
        if rebin:
            # variable binning
            h_num = h_Data.Rebin(self.n_bins, "h_num", self.xbins)
            h_den = h_mc.Rebin(self.n_bins, "h_den", self.xbins)
        else:
            h_num = h_Data.Clone("h_num") 
            h_den = h_mc.Clone("h_den")
            # contstant binning for plots
            h_num.Rebin(2)
            h_den.Rebin(2)
        h_ratio_normalized = getNormalizedRatio(h_num, h_den)
        return h_ratio_normalized

    def makeZvsPhoton(self, file_name, era, rebin):
        doFit = False
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
            h_ratio_lepton = self.getZRatio(f, region, selection, era, metLepton, rebin)
            # MET with photon
            # getPhotonRatio(self, region, selection, era, variable)
            h_ratio_photon = self.getPhotonRatio(f, region, selection, era, metPhoton, rebin)
            
            # Double Ratio for MET: Z / photon
            h_ratio_ZoverPhoton = h_ratio_lepton.Clone("h_ratio_ZoverPhoton")
            h_ratio_ZoverPhoton.Divide(h_ratio_photon)    
                
            # fit the Z over Photon ratio
            if doFit:
                # the first MET bin is not fit
                # we use a linear 2 parameter fit
                nBinsFit = self.n_bins - 1
                nDegFree = nBinsFit - 2
                fit = ROOT.TF1("f1", "pol1", self.met_min, self.met_max)
                h_ratio_ZoverPhoton.Fit(fit, "N", "", self.met_min, self.met_max)
                fit.SetLineColor(getColorIndex("violet"))
                fit.SetLineWidth(5)
                p0      = fit.GetParameter(0)
                p1      = fit.GetParameter(1)
                p0_err  = fit.GetParError(0)
                p1_err  = fit.GetParError(1)
                chisq   = fit.GetChisquare()
                chisq_r = chisq / nDegFree
                mark = ROOT.TLatex()
                mark.SetTextSize(0.05)
            
            
            # histogram info 
            title = "Z vs. Photon, {0}, {1}".format(region, era)
            x_title = "MET (GeV)" 
            y_title = "Data / MC"
            y_min = 0.0
            y_max = 2.0
            
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h_ratio_lepton,       title, x_title, y_title,                "vermillion",      y_min, y_max)
            setupHist(h_ratio_photon,       title, x_title, y_title,                "electric blue",   y_min, y_max)
            setupHist(h_ratio_ZoverPhoton,  title, x_title, "(Z to LL) / Photon",   "black",           y_min, y_max)
            # set x axis range
            h_ratio_lepton.GetXaxis().SetRangeUser(self.met_min, self.met_max)
            h_ratio_photon.GetXaxis().SetRangeUser(self.met_min, self.met_max)
            h_ratio_ZoverPhoton.GetXaxis().SetRangeUser(self.met_min, self.met_max)
            
            # do Run 2 systematic
            if era == "Run2" and rebin:
                h_syst = ROOT.TH1F("h_syst", "h_syst", self.n_bins, self.xbins)
                for i in xrange(1, self.n_bins + 1):
                    # syst = max(stat uncertainty in double ratio, |(double ratio) - 1|)
                    value    = h_ratio_ZoverPhoton.GetBinContent(i)
                    stat_err = h_ratio_ZoverPhoton.GetBinError(i)
                    diff     = abs(value - 1)
                    syst_err = max(stat_err, diff) 
                    h_syst.SetBinContent(i, syst_err)
                    h_syst.SetBinError(i, 0)
                setupHist(h_syst,       title, x_title, "syst.",   "irish green",      y_min, y_max)
                h_syst.GetXaxis().SetRangeUser(self.met_min, self.met_max)
                self.h_map_syst[region] = copy.deepcopy(h_syst)
            
            # pad for histograms
            pad = c.cd(1)
            pad.SetGrid()
            
            # draw
            h_ratio_lepton.Draw(draw_option)
            h_ratio_photon.Draw(draw_option + " same")
            # legend: TLegend(x1,y1,x2,y2)
            legend1 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend1.AddEntry(h_ratio_lepton,         "Z to LL Data/MC",              "l")
            legend1.AddEntry(h_ratio_photon,         "Photon Data/MC",               "l")
            legend1.Draw()
            
            # pad for ratio
            pad = c.cd(2)
            pad.SetGrid()
            
            # draw
            h_ratio_ZoverPhoton.Draw(draw_option)
            if doFit:
                fit.Draw("same")
                # write chisq_r
                # give x, y coordinates (same as plot coordinates)
                #print "Fit: f(x) = (%.5f #pm %.5f) * x + (%.5f #pm %.5f)" % (p1, p1_err, p0, p0_err)
                mark.DrawLatex(300.0, y_max - 0.2, "Fit: f(x) = %.5f + %.5f * x" % (p0, p1))
                mark.DrawLatex(300.0, y_max - 0.4, "#chi_{r}^{2} = %.3f" % chisq_r)
            if era == "Run2" and rebin:
                h_syst.Draw(draw_option + " same")
            
            # legend: TLegend(x1,y1,x2,y2)
            legend2 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend2.AddEntry(h_ratio_ZoverPhoton,    "(Z to LL) / Photon",           "l")
            if doFit:
                legend2.AddEntry(fit,                "Fit to (Z to LL) / Photon",    "l")
            if era == "Run2" and rebin:
                legend2.AddEntry(h_syst,             "syst. unc.",    "l")
            legend2.Draw()
            
            # save histograms
            if rebin:
                plot_name = "{0}ZvsPhoton_{1}_rebinned_{2}".format(self.plot_dir, region, era)
            else:
                plot_name = "{0}ZvsPhoton_{1}_{2}".format(self.plot_dir, region, era)
            c.Update()
            c.SaveAs(plot_name + ".pdf")
            c.SaveAs(plot_name + ".png")



