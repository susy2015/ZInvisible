# systematics.py

import copy
from colors import getColorIndex
import json
import os
import numpy as np
import ROOT
import tools
from tools import setupHist, getRatio, getNormalizedRatio, normalizeHistToHist

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
        self.x_min = 250.0
        self.x_max = 1000.0
        self.h_map_syst = {}

    # get make of Z histograms
    def getZHistoMap(self, region, nameTag, selectionTag, variable):
        # histogram name examples
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLDatadata;1  metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLDYstack;1   metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLt#bar{t}stack;1 metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLSingle tstack;1 metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30_2016metWithLLmetWithLLRarestack;1 metWithLL
        debug = False
        h_map_norm = {}
        for particle in self.particles:
            h_map_norm[particle] = { 
                "Data"      : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Datadata",
                "DY"        : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "DYstack",
                "TTbar"     : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "t#bar{t}stack",
                "SingleT"   : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Single tstack",
                "Rare"      : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Rarestack"
            }
        
            if debug:
                print "# Loading {0} histograms".format(particle)
                print str(variable + "/" + h_map_norm[particle]["Data"]     )
                print str(variable + "/" + h_map_norm[particle]["DY"]       )
                print str(variable + "/" + h_map_norm[particle]["TTbar"]    )
                print str(variable + "/" + h_map_norm[particle]["SingleT"]  )
                print str(variable + "/" + h_map_norm[particle]["Rare"]     )
        return h_map_norm

    def getHistRatio(self, num, den, rebin):
        # WARNING: do not rebin ratios; first rebin, then get ratio
        if rebin:
            # variable binning
            h_num = num.Rebin(self.n_bins, "h_num", self.xbins)
            h_den = den.Rebin(self.n_bins, "h_den", self.xbins)
        else:
            h_num = num.Clone("h_num") 
            h_den = den.Clone("h_den")
        #h_ratio = getNormalizedRatio(h_num, h_den)
        h_ratio = getRatio(h_num, h_den)
        return h_ratio
    
    # WARNING: strings loaded from json file have type 'unicode'
    # ROOT cannot load histograms using unicode input: use type 'str'

    def getZData(self, root_file, variable, h_map_norm):
        h_data_electron = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Data"]    ) )
        h_data_muon     = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Data"]        ) )
        h_data_lepton   = h_data_electron.Clone("h_data_lepton") 
        h_data_lepton.Add(h_data_muon)
        return h_data_lepton

    def getZMC(self, root_file, variable, h_map_norm):
        h_DY_Electron       = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["DY"]      ) )
        h_TTbar_Electron    = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["TTbar"]   ) )
        h_SingleT_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["SingleT"] ) )
        h_Rare_Electron     = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Rare"]    ) )
        h_DY_Muon           = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["DY"]          ) )
        h_TTbar_Muon        = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["TTbar"]       ) )
        h_SingleT_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["SingleT"]     ) )
        h_Rare_Muon         = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Rare"]        ) )
        
        h_mc_Z = h_DY_Electron.Clone("h_mc_Z") 
        h_mc_Z.Add(h_TTbar_Electron)
        h_mc_Z.Add(h_SingleT_Electron)
        h_mc_Z.Add(h_Rare_Electron)
        h_mc_Z.Add(h_DY_Muon)
        h_mc_Z.Add(h_TTbar_Muon)
        h_mc_Z.Add(h_SingleT_Muon)
        h_mc_Z.Add(h_Rare_Muon)
        
        return h_mc_Z

    def getPhotonData(self, root_file, variable, h_map_shape):
        h_data_photon = root_file.Get( str(variable + "/" + h_map_shape["Data"]) )
        return h_data_photon
    
    def getPhotonMC(self, root_file, variable, h_map_shape):
        h_GJets             = root_file.Get( str(variable + "/" + h_map_shape["GJets"]           ) )
        if self.S.splitQCD:
            h_QCD_Direct        = root_file.Get( str(variable + "/" + h_map_shape["QCD_Direct"]         ) )
            h_QCD_Fragmentation = root_file.Get( str(variable + "/" + h_map_shape["QCD_Fragmentation"]  ) )
            h_QCD_NonPrompt     = root_file.Get( str(variable + "/" + h_map_shape["QCD_NonPrompt"]      ) )
            h_QCD_Fake          = root_file.Get( str(variable + "/" + h_map_shape["QCD_Fake"]           ) )
        else:
            h_QCD               = root_file.Get( str(variable + "/" + h_map_shape["QCD"]  ) )
        h_WJets             = root_file.Get( str(variable + "/" + h_map_shape["WJets"]    ) )
        h_TTG               = root_file.Get( str(variable + "/" + h_map_shape["TTG"]      ) )
        h_tW                = root_file.Get( str(variable + "/" + h_map_shape["tW"]       ) )
        h_Rare              = root_file.Get( str(variable + "/" + h_map_shape["Rare"]     ) )
        
        h_mc_photon = h_GJets.Clone("h_mc_photon") 
        if self.S.splitQCD:
            h_mc_photon.Add(h_QCD_Direct)
            h_mc_photon.Add(h_QCD_Fragmentation)
            h_mc_photon.Add(h_QCD_NonPrompt)
            h_mc_photon.Add(h_QCD_Fake)
        else:
            h_mc_photon.Add(h_QCD)
        h_mc_photon.Add(h_WJets)
        h_mc_photon.Add(h_TTG)
        h_mc_photon.Add(h_tW)
        h_mc_photon.Add(h_Rare)

        return h_mc_photon

    # get data Z / Photon ratio
    def getDataRatio(self, root_file, region, selection, name, varLepton, varPhoton, rebin):
        selectionTag    = "_" + selection
        nameTag         = "_" + name
        h_map_norm  = self.getZHistoMap(region, nameTag, selectionTag, varLepton)
        h_map_shape = self.S.getSimpleMap(region, nameTag, selectionTag, selectionTag, varPhoton)
        # numerator: Z data
        h_data_lepton = self.getZData(root_file, varLepton, h_map_norm)
        # denominator: Photon data
        h_data_photon = self.getPhotonData(root_file, varPhoton, h_map_shape)
        # take ratio
        h_ratio = self.getHistRatio(h_data_lepton, h_data_photon, rebin) 
        return h_ratio

    # get MC Z / Photon ratio
    def getMCRatio(self, root_file, region, selection, name, varLepton, varPhoton, rebin):
        selectionTag    = "_" + selection
        nameTag         = "_" + name
        h_map_norm  = self.getZHistoMap(region, nameTag, selectionTag, varLepton)
        h_map_shape = self.S.getSimpleMap(region, nameTag, selectionTag, selectionTag, varPhoton)
        # numerator: Z MC (normalize to Z data)
        h_mc_Z   = self.getZMC(root_file, varLepton, h_map_norm)
        h_data_Z = self.getZData(root_file, varLepton, h_map_norm)
        normalizeHistToHist(h_mc_Z, h_data_Z)
        # denominator: Photon MC (normalize to Photon data)
        h_mc_photon   = self.getPhotonMC(root_file, varPhoton, h_map_shape)
        h_data_photon = self.getPhotonData(root_file, varPhoton, h_map_shape)
        normalizeHistToHist(h_mc_photon, h_data_photon)
        # take ratio
        h_ratio = self.getHistRatio(h_mc_Z, h_mc_photon, rebin) 
        return h_ratio

    def getZRatio(self, root_file, region, selection, name, variable, rebin):
        selectionTag    = "_" + selection
        nameTag         = "_" + name
        h_map_norm = self.getZHistoMap(region, nameTag, selectionTag, variable)
        # numerator: data
        h_data = self.getZData(root_file, variable, h_map_norm)
        # denominator: MC 
        h_mc = self.getZMC(root_file, variable, h_map_norm)
        # normalize MC to data
        normalizeHistToHist(h_mc, h_data)
        # take ratio
        h_ratio = self.getHistRatio(h_data, h_mc, rebin)
        return h_ratio

    def getPhotonRatio(self, root_file, region, selection, name, variable, rebin):
        selectionTag    = "_" + selection
        nameTag         = "_" + name
        # getSimpleMap(self, region, nameTag, dataSelectionTag, mcSelectionTag, variable)
        h_map_shape = self.S.getSimpleMap(region, nameTag, selectionTag, selectionTag, variable)
        # numerator: data
        h_data = self.getPhotonData(root_file, variable, h_map_shape)
        # denominator: MC
        h_mc = self.getPhotonMC(root_file, variable, h_map_shape)
        # normalize MC to data
        normalizeHistToHist(h_mc, h_data)
        # take ratio
        h_ratio = self.getHistRatio(h_data, h_mc, rebin)
        return h_ratio

    # double ratio: Z (data/MC) over Photon (data/MC), or data (Z/Photon) over MC (Z/Photon)
    def makeZvsPhoton(self, file_name, var, varPhoton, varLepton, era, rebin, useForSyst, doDataOverData=False, xbins = np.array([]), n_bins = 0, x_min=0, x_max=0):
        # set variables based on mode
        if doDataOverData:
            fileTag = "DataOverData"
            num_label = "(Z to LL Data)/(Photon Data)"
            den_label = "(Z to LL MC)/(Photon MC)" 
            y_min = 0.0
            y_max = 0.2
        else:
            fileTag = "ZvsPhoton"
            num_label = "Z to LL Data/MC"
            den_label = "Photon Data/MC" 
            y_min = 0.0
            y_max = 2.0
        # redefine xbins and n_bins if provided
        if xbins.any():
            self.xbins  = xbins
            self.n_bins = n_bins
        if x_max:
            self.x_min = x_min
            self.x_max = x_max
        doFit = False
        draw_option = "hist error"
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        f = ROOT.TFile(file_name)
        c = ROOT.TCanvas("c", "c", 800, 800)
        c.Divide(1, 2)
        selection = "jetpt30"
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 
        for region in self.regions:
            if doDataOverData:
                h_ratio_num = self.getDataRatio(f, region, selection, var, varLepton, varPhoton, rebin)
                h_ratio_den = self.getMCRatio(f, region, selection, var, varLepton, varPhoton, rebin)
            else:
                h_ratio_num = self.getZRatio(f, region, selection, var, varLepton, rebin)
                h_ratio_den = self.getPhotonRatio(f, region, selection, var, varPhoton, rebin)
            
            # Double Ratio: Z / photon
            h_ratio_ZoverPhoton = h_ratio_num.Clone("h_ratio_ZoverPhoton")
            h_ratio_ZoverPhoton.Divide(h_ratio_den)    
                
            # fit the Z over Photon ratio
            if doFit:
                # the first MET bin is not fit
                # we use a linear 2 parameter fit
                nBinsFit = self.n_bins - 1
                nDegFree = nBinsFit - 2
                fit = ROOT.TF1("f1", "pol1", self.x_min, self.x_max)
                h_ratio_ZoverPhoton.Fit(fit, "N", "", self.x_min, self.x_max)
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
            x_title = var
            y_title = "Data / MC"
            ratio_y_min = 0.0
            ratio_y_max = 2.0
            
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h_ratio_num,          title, x_title, y_title,                "vermillion",      y_min, y_max)
            setupHist(h_ratio_den,          title, x_title, y_title,                "electric blue",   y_min, y_max)
            setupHist(h_ratio_ZoverPhoton,  title, x_title, "(Z to LL) / Photon",   "black",           ratio_y_min, ratio_y_max)
            # set x axis range
            #print "self.x_min = {0}".format(self.x_min)
            #print "self.x_max = {0}".format(self.x_max)
            h_ratio_num.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            h_ratio_den.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            h_ratio_ZoverPhoton.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            
            # do Run 2 systematic
            if era == "Run2":
                #h_syst = ROOT.TH1F("h_syst", "h_syst", self.n_bins, self.xbins)
                h_syst = h_ratio_ZoverPhoton.Clone("h_syst")
                for i in xrange(1, h_syst.GetNbinsX() + 1):
                    # syst = max(stat uncertainty in double ratio, |(double ratio) - 1|)
                    value    = h_ratio_ZoverPhoton.GetBinContent(i)
                    stat_err = h_ratio_ZoverPhoton.GetBinError(i)
                    # default
                    syst_err = 0
                    if value != 0:
                        diff     = abs(value - 1)
                        syst_err = max(stat_err, diff) 
                    h_syst.SetBinContent(i, syst_err)
                    h_syst.SetBinError(i, 0)
                setupHist(h_syst,       title, x_title, "syst.",   "irish green",      y_min, y_max)
                h_syst.GetXaxis().SetRangeUser(self.x_min, self.x_max)
                if useForSyst:
                    self.h_map_syst[region] = copy.deepcopy(h_syst)
            
            # pad for histograms
            pad = c.cd(1)
            pad.SetGrid()
            
            # draw
            h_ratio_num.Draw(draw_option)
            h_ratio_den.Draw(draw_option + " same")
            # legend: TLegend(x1,y1,x2,y2)
            legend1 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend1.AddEntry(h_ratio_num, num_label, "l")
            legend1.AddEntry(h_ratio_den, den_label, "l")
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
            if era == "Run2":
                h_syst.Draw(draw_option + " same")
            
            # legend: TLegend(x1,y1,x2,y2)
            legend2 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend2.AddEntry(h_ratio_ZoverPhoton,    "(Z to LL) / Photon",           "l")
            if doFit:
                legend2.AddEntry(fit,                "Fit to (Z to LL) / Photon",    "l")
            if era == "Run2":
                legend2.AddEntry(h_syst,             "syst. unc.",    "l")
            legend2.Draw()
            
            # save histograms
            if rebin:
                plot_name = "{0}{1}_{2}_{3}_rebinned_{4}".format(self.plot_dir, fileTag, var, region, era)
            else:
                plot_name = "{0}{1}_{2}_{3}_{4}".format(self.plot_dir, fileTag, var, region, era)
            c.Update()
            c.SaveAs(plot_name + ".pdf")
            c.SaveAs(plot_name + ".png")



