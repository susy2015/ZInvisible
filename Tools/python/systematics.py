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
        self.region_labels = {"LowDM" : "Low #Deltam", "HighDM": "High #Deltam"}
        self.labels = {
                        "met"   : "Modified p_{T}^{miss} [GeV]", 
                        "ht"    : "H_{T} [GeV]",
                        "nj"    : "N_{j}",
                        "nb"    : "N_{b}",
                        "nmt"   : "N_{merged tops}",
                        "nrt"   : "N_{resolved tops}",
                        "nw"    : "N_{W}"
                      }
        self.lumis = {
                        "2016" :  35815.165,
                        "2017" :  41486.136,
                        "2018" :  59699.489,
                        "Run2" : 137000.790,
                        "2016and2017" : 77301.301
                     }

    # get make of Z histograms
    def getZHistoMap(self, region, nameTag, selectionTag, variable):
        # histogram name examples
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLDatadata;1   metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLDYstack;1    metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLDibosonstack;1   metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLRarestack;1  metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLt#bar{t}stack;1  metWithLL
        # KEY: TH1D  DataMC_Electron_LowDM_met_jetpt30metWithLLmetWithLLSingle tstack;1  metWithLL  
        
        debug = False
        h_map_norm = {}
        for particle in self.particles:
            h_map_norm[particle] = { 
                "Data"      : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Datadata",
                "DY"        : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "DYstack",
                "Diboson"   : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Dibosonstack",
                "Rare"      : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Rarestack",
                "TTbar"     : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "t#bar{t}stack",
                "SingleT"   : "DataMC_" + particle + "_" + region + nameTag + selectionTag + 2 * variable + "Single tstack"
            }
        
            if debug:
                print "# Loading {0} histograms".format(particle)
                print str(variable + "/" + h_map_norm[particle]["Data"]     )
                print str(variable + "/" + h_map_norm[particle]["DY"]       )
                print str(variable + "/" + h_map_norm[particle]["Diboson"]  )
                print str(variable + "/" + h_map_norm[particle]["Rare"]     )
                print str(variable + "/" + h_map_norm[particle]["TTbar"]    )
                print str(variable + "/" + h_map_norm[particle]["SingleT"]  )
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
        h_Diboson_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Diboson"] ) )
        h_Rare_Electron     = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Rare"]    ) )
        h_TTbar_Electron    = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["TTbar"]   ) )
        h_SingleT_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["SingleT"] ) )
        h_DY_Muon           = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["DY"]          ) )
        h_Diboson_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Diboson"]     ) )
        h_Rare_Muon         = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Rare"]        ) )
        h_TTbar_Muon        = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["TTbar"]       ) )
        h_SingleT_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["SingleT"]     ) )
        
        h_mc_Z = h_DY_Electron.Clone("h_mc_Z") 
        h_mc_Z.Add(h_Diboson_Electron)
        h_mc_Z.Add(h_Rare_Electron)
        h_mc_Z.Add(h_TTbar_Electron)
        h_mc_Z.Add(h_SingleT_Electron)
        h_mc_Z.Add(h_DY_Muon)
        h_mc_Z.Add(h_Diboson_Muon)
        h_mc_Z.Add(h_Rare_Muon)
        h_mc_Z.Add(h_TTbar_Muon)
        h_mc_Z.Add(h_SingleT_Muon)
        
        return h_mc_Z
    
    def getZMCHists(self, root_file, variable, h_map_norm):
        h_DY_Electron       = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["DY"]      ) )
        h_Diboson_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Diboson"] ) )
        h_Rare_Electron     = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["Rare"]    ) )
        h_TTbar_Electron    = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["TTbar"]   ) )
        h_SingleT_Electron  = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["SingleT"] ) )
        h_DY_Muon           = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["DY"]          ) )
        h_Diboson_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Diboson"]     ) )
        h_Rare_Muon         = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["Rare"]        ) )
        h_TTbar_Muon        = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["TTbar"]       ) )
        h_SingleT_Muon      = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["SingleT"]     ) )
        
        h_DY        = h_DY_Electron.Clone("h_DY") 
        h_Diboson   = h_Diboson_Electron.Clone("h_Diboson") 
        h_Rare      = h_Rare_Electron.Clone("h_Rare") 
        h_TTbar     = h_TTbar_Electron.Clone("h_TTbar") 
        h_SingleT   = h_SingleT_Electron.Clone("h_SingleT") 
        h_DY.Add(h_DY_Muon)
        h_Diboson.Add(h_Diboson_Muon)
        h_Rare.Add(h_Rare_Muon)
        h_TTbar.Add(h_TTbar_Muon)
        h_SingleT.Add(h_SingleT_Muon)
        
        hist_map = {}
        hist_map["DY"]          = h_DY
        hist_map["Diboson"]     = h_Diboson
        hist_map["Rare"]        = h_Rare
        hist_map["TTbar"]       = h_TTbar
        hist_map["SingleT"]     = h_SingleT
        
        return hist_map
    
    def getDY(self, root_file, variable, h_map_norm):
        h_DY_Electron   = root_file.Get( str(variable + "/" + h_map_norm["Electron"]["DY"]      ) )
        h_DY_Muon       = root_file.Get( str(variable + "/" + h_map_norm["Muon"]["DY"]          ) )
        h_DY            = h_DY_Electron.Clone("h_DY") 
        h_DY.Add(h_DY_Muon)
        return h_DY

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
    
    def getPhotonMCHists(self, root_file, variable, h_map_shape):
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
        
        hist_map = {}
        hist_map["GJets"] = h_GJets
        if self.S.splitQCD:
            hist_map["QCD_Direct"]          = h_QCD_Direct
            hist_map["QCD_Fragmentation"]   = h_QCD_Fragmentation
            hist_map["QCD_NonPrompt"]       = h_QCD_NonPrompt
            hist_map["QCD_Fake"]            = h_QCD_Fake
        else:
            hist_map["QCD"] = h_QCD
        hist_map["WJets"] = h_WJets
        hist_map["TTG"]   = h_TTG
        hist_map["tW"]    = h_tW
        hist_map["Rare"]  = h_Rare

        return hist_map
    
    def getGJets(self, root_file, variable, h_map_shape):
        h_GJets = root_file.Get( str(variable + "/" + h_map_shape["GJets"]           ) )
        return h_GJets

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
    def getMCRatio(self, root_file, region, selection, name, varLepton, varPhoton, rebin, val_map):
        selectionTag    = "_" + selection
        nameTag         = "_" + name
        h_map_norm      = self.getZHistoMap(region, nameTag, selectionTag, varLepton)
        h_map_shape     = self.S.getSimpleMap(region, nameTag, selectionTag, selectionTag, varPhoton)
        # WARNING: first calculate purity; normalize histograms at the end
        # numerator: Z MC (normalize to Z data)
        h_DY     = self.getDY(root_file, varLepton, h_map_norm)
        h_mc_Z   = self.getZMC(root_file, varLepton, h_map_norm)
        h_data_Z = self.getZData(root_file, varLepton, h_map_norm)
        # denominator: Photon MC (normalize to Photon data)
        h_GJets       = self.getGJets(root_file, varPhoton, h_map_shape)
        h_mc_photon   = self.getPhotonMC(root_file, varPhoton, h_map_shape)
        h_data_photon = self.getPhotonData(root_file, varPhoton, h_map_shape)
        # add purity values to map
        num_DY          = h_DY.Integral(0, h_DY.GetNbinsX() + 1)
        num_mc_Z        = h_mc_Z.Integral(0, h_mc_Z.GetNbinsX() + 1)
        num_GJets       = h_GJets.Integral(0, h_GJets.GetNbinsX() + 1)
        num_mc_photon   = h_mc_photon.Integral(0, h_mc_photon.GetNbinsX() + 1)
        purity_DY       = num_DY / num_mc_Z
        purity_GJets    = num_GJets / num_mc_photon
        val_map[region]["num_DY"]           = num_DY
        val_map[region]["num_mc_Z"]         = num_mc_Z
        val_map[region]["num_GJets"]        = num_GJets
        val_map[region]["num_mc_photon"]    = num_mc_photon
        val_map[region]["purity_DY"]        = purity_DY
        val_map[region]["purity_GJets"]     = purity_GJets
        # add norm values to map
        norm_Z      = normalizeHistToHist(h_mc_Z, h_data_Z)
        norm_photon = normalizeHistToHist(h_mc_photon, h_data_photon)
        val_map[region]["norm_Z"]      = norm_Z 
        val_map[region]["norm_photon"] = norm_photon
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
    def makeZvsPhoton(self, file_name, var, varPhoton, varLepton, era, x_min, x_max, n_bins=0, xbins=np.array([]), rebin=False, useForSyst=False, doDataOverData=False):
        doFit = False
        draw_option = "hist"
        data_style  = "PZ0"
        drawSyst = (era == "Run2") and (not doDataOverData)
        ROOT.gStyle.SetErrorX(0) # remove horizontal bar for data points
        
        # set lower pad height as percentage
        lowerPadHeight = 0.30
        padHeightRatio = lowerPadHeight / (1.0 - lowerPadHeight)
        
        # set variables based on mode
        if doDataOverData:
            fileTag     = "DataOverData"
            num_label   = "Data"
            den_label   = "Simulation"
            y_title     = "(Z #rightarrow ll)/#gamma"
            ratio_title = "Data/Sim."
            num_color   = "black"
            num_legend_style = "pe"
            # avoid 0 label which is cutoff
            y_min = 0.01
            y_max = 0.20
        else:
            fileTag     = "ZvsPhoton"
            num_label   = "Z #rightarrow ll"
            den_label   = "#gamma" 
            y_title     = "Data/Sim."
            ratio_title = "(Z #rightarrow ll)/#gamma"
            num_color   = "vermillion"
            num_legend_style = "l"
            # avoid 0 label which is cutoff
            y_min = 0.1
            y_max = 2.0
        
        # redefine xbins and n_bins if provided
        if xbins.any():
            self.xbins  = xbins
            self.n_bins = n_bins
        self.x_min = x_min
        self.x_max = x_max
        
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        f = ROOT.TFile(file_name)
        c = ROOT.TCanvas("c", "c", 800, 800)
        c.Divide(1, 2)
        selection    = "jetpt30"
        output_name  = "output/values_{0}_{1}.json".format(var, era)
        val_map = {}
        for region in self.regions:
            val_map[region] = {}
            if doDataOverData:
                h_ratio_num = self.getDataRatio(f, region, selection, var, varLepton, varPhoton, rebin)
                h_ratio_den = self.getMCRatio(f, region, selection, var, varLepton, varPhoton, rebin, val_map)
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
            
            # histogram info 
            # turn off title
            #title = "Z vs. Photon, {0}, {1}".format(region, era)
            title = ""
            ratio_y_min = 0.7
            ratio_y_max = 1.3
            x_title = var
            if var in self.labels:
                x_title = self.labels[var]
            
            # turn off x-axis titles for upper plot
            # setupHist(hist, title, x_title, y_title, color, y_min, y_max, adjust=False, lineWidth=5, turnOffStats=True)
            setupHist(h_ratio_num,          title, "",      y_title,       num_color,         y_min,        y_max,        True,  1)
            setupHist(h_ratio_den,          title, "",      y_title,       "electric blue",   y_min,        y_max,        True,  3)
            setupHist(h_ratio_ZoverPhoton,  title, x_title, ratio_title,   "black",           ratio_y_min,  ratio_y_max,  True,  1)
            h_ratio_num.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            h_ratio_den.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            h_ratio_ZoverPhoton.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            
            # label and title formatting
            labelSize           = 0.14
            titleSize           = 0.14
            titleOffsetXaxis    = 1.20
            titleOffsetYaxis    = 0.60
            
            # upper plot
            h_ratio_den.GetXaxis().SetLabelSize(0) # turn off x-axis labels for upper plot
            h_ratio_den.GetYaxis().SetLabelSize(padHeightRatio * labelSize)
            h_ratio_den.GetYaxis().SetTitleSize(padHeightRatio * titleSize)
            h_ratio_den.GetYaxis().SetTitleOffset(titleOffsetYaxis/padHeightRatio)
            h_ratio_den.GetYaxis().SetNdivisions(5, 5, 0, True)
            
            # lower plot
            h_ratio_ZoverPhoton.GetXaxis().SetLabelSize(labelSize)
            h_ratio_ZoverPhoton.GetXaxis().SetTitleSize(titleSize)
            h_ratio_ZoverPhoton.GetXaxis().SetTitleOffset(titleOffsetXaxis)
            h_ratio_ZoverPhoton.GetXaxis().SetNdivisions(5, 5, 0, True)
            h_ratio_ZoverPhoton.GetYaxis().SetLabelSize(labelSize)
            h_ratio_ZoverPhoton.GetYaxis().SetTitleSize(titleSize)
            h_ratio_ZoverPhoton.GetYaxis().SetTitleOffset(titleOffsetYaxis)
            h_ratio_ZoverPhoton.GetYaxis().SetNdivisions(3, 5, 0, True)
            
            # do Run 2 systematic
            if era == "Run2":
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
                # setup
                setupHist(h_syst,       title, x_title, "syst.",   "irish green",      y_min, y_max)
                h_syst.GetXaxis().SetRangeUser(self.x_min, self.x_max)
                if useForSyst:
                    self.h_map_syst[region] = copy.deepcopy(h_syst)
            
            # manually (bin by bin) create TGraphErrors to show stat. unc. as error bands (rectangles)
            
            nBins      = h_ratio_den.GetNbinsX()
            # simulation (Z/photon) values and errors (upper plot)
            xValsSim   = []
            yValsSim   = []
            xErrorsSim = []
            yErrorsSim = []
            # relative stat. unc. for simulation (lower plot)
            xValsSimRel   = []
            yValsSimRel   = []
            xErrorsSimRel = []
            yErrorsSimRel = []
            for n in range(1, nBins + 1):
                # center at bin center
                x_val_sim = h_ratio_den.GetBinCenter(n)
                y_val_sim = h_ratio_den.GetBinContent(n)
                xValsSim.append(x_val_sim)
                yValsSim.append(y_val_sim)
                xValsSimRel.append(x_val_sim)
                yValsSimRel.append(1.0) # WARNING: needs to be a float (1.0), not int (1) to work
                # use 1/2 of bin width for x_error to draw rectangle errors
                x_error_sim  = h_ratio_den.GetBinWidth(n) / 2.0
                y_error_sim  = h_ratio_den.GetBinError(n)
                y_error_data = h_ratio_num.GetBinError(n)
                # relative uncertainties for ratio plot
                simRelUnc  = 0.0 
                dataRelUnc = 0.0 
                if y_val_sim != 0.0:
                    simRelUnc  = y_error_sim  / y_val_sim
                    dataRelUnc = y_error_data / y_val_sim
                # stat. unc. for simulation 
                xErrorsSim.append(x_error_sim)
                yErrorsSim.append(y_error_sim)
                xErrorsSimRel.append(x_error_sim)
                yErrorsSimRel.append(simRelUnc)
                # relative stat. unc. for data (lower plot)
                h_ratio_ZoverPhoton.SetBinError(n, dataRelUnc)
            
            # WARNING: np.array() of floats is required for TGraphErrors constructor
            # - needs to be an array of floats
            # - can't be a list of floats
            # - can't be an array of ints
            xValsSim        = np.array(xValsSim)
            yValsSim        = np.array(yValsSim)
            xErrorsSim      = np.array(xErrorsSim)
            yErrorsSim      = np.array(yErrorsSim)
            xValsSimRel     = np.array(xValsSimRel)
            yValsSimRel     = np.array(yValsSimRel)
            xErrorsSimRel   = np.array(xErrorsSimRel)
            yErrorsSimRel   = np.array(yErrorsSimRel)
            
            # simulation (Z/photon) values and errors (upper plot)
            sim_stat_unc = ROOT.TGraphErrors(nBins, xValsSim, yValsSim, xErrorsSim, yErrorsSim)
            sim_stat_unc.SetFillColor(getColorIndex("electric blue"))
            sim_stat_unc.SetFillStyle(3013)
            sim_stat_unc.SetLineStyle(0)
            sim_stat_unc.SetLineWidth(0)
            sim_stat_unc.SetMarkerSize(0)
            
            # relative stat. unc. for simulation (lower plot)
            sim_rel_stat_unc = ROOT.TGraphErrors(nBins, xValsSimRel, yValsSimRel, xErrorsSimRel, yErrorsSimRel)
            sim_rel_stat_unc.SetFillColor(getColorIndex("electric blue"))
            sim_rel_stat_unc.SetFillStyle(3013)
            sim_rel_stat_unc.SetLineStyle(0)
            sim_rel_stat_unc.SetLineWidth(0)
            sim_rel_stat_unc.SetMarkerSize(0)
            
            # pad for histograms
            pad = c.cd(1)
            # resize pad
            # SetPad(xlow, ylow, xup, yup)
            pad.SetPad(0, lowerPadHeight, 1, 1)
            # set ticks on all sides of plot
            pad.SetTickx()
            pad.SetTicky()
            pad.SetLeftMargin(0.2)
            pad.SetRightMargin(0.1)
            pad.SetTopMargin(0.1)
            pad.SetBottomMargin(0.01)

            # --- testing ---
            printTests = False

            if printTests:

                for n in range(1, h_ratio_den.GetNbinsX() + 1):
                    print "h_ratio_den_{0}_{1}_{2}: bin: {3}, content: {4}, error: {5}".format(era, var, region, n, h_ratio_den.GetBinContent(n), h_ratio_den.GetBinError(n))
                
                points_x = sim_stat_unc.GetX()
                points_y = sim_stat_unc.GetY()
                print "TGraphErrors_{0}: sim_stat_unc.GetN(): {1}".format(era, sim_stat_unc.GetN())
                for n in range(1, sim_stat_unc.GetN() + 1):
                    print "TGraphErrors_{0}_{1}_{2}_upper_plot: bin: {3}, x: {4}, y: {5}, x_error: {6}, y_error: {7}".format(era, var, region, n, points_x[n-1], points_y[n-1], sim_stat_unc.GetErrorX(n), sim_stat_unc.GetErrorY(n))
                
                points_x = sim_rel_stat_unc.GetX()
                points_y = sim_rel_stat_unc.GetY()
                print "TGraphErrors_{0}: sim_rel_stat_unc.GetN(): {1}".format(era, sim_rel_stat_unc.GetN())
                for n in range(1, sim_rel_stat_unc.GetN() + 1):
                    print "TGraphErrors_{0}_{1}_{2}_lower_plot: bin: {3}, x: {4}, y: {5}, x_error: {6}, y_error: {7}".format(era, var, region, n, points_x[n-1], points_y[n-1], sim_rel_stat_unc.GetErrorX(n), sim_rel_stat_unc.GetErrorY(n))
            
            # --- draw --- 
            h_ratio_den.Draw(draw_option)
            sim_stat_unc.Draw("E2 same") # E2: error rectangles
            if doDataOverData:
                h_ratio_num.SetMarkerStyle(ROOT.kFullDotLarge)
                h_ratio_num.Draw(data_style + " same")
            else:
                h_ratio_num.Draw(draw_option + " same")
            
            # legend: TLegend(x1,y1,x2,y2)
            legend_x1 = 0.60
            legend_x2 = 0.90
            legend_y1 = 0.65
            legend_y2 = 0.85
            legend1 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend1.SetFillStyle(0)
            legend1.SetBorderSize(0)
            legend1.SetLineWidth(1)
            legend1.SetNColumns(1)
            legend1.SetTextFont(42)
            legend1.AddEntry(h_ratio_num, num_label, num_legend_style)
            legend1.AddEntry(h_ratio_den, den_label, "l")
            legend1.AddEntry(sim_stat_unc, "Stat. unc.", "F")
            legend1.Draw()
            
            # parameters for CMS mark, other text, and lumi stamp
            font_scale  = 0.05
            left        = self.x_min
            right       = self.x_max
            x_offset    = 0.06  # x offset from left edge
            y_offset    = 0.90  # y offset from bottom edge
            y_spacing   = 0.075 # vertical spacing between marks
            width       = right - left
            mark_x1     = left + x_offset * width
            mark_x2     = left + (x_offset + 0.15) * width
            mark_y1     = (y_offset - 0 * y_spacing) * y_max
            mark_y2     = (y_offset - 1 * y_spacing) * y_max
            mark_y3     = (y_offset - 2 * y_spacing) * y_max
            lumi_x      = right 
            lumi_y      = 1.02 * y_max
            
            # Draw CMS and other text
            cms_mark = ROOT.TLatex()
            cms_mark.SetTextAlign(11) # left aligned
            cms_mark.SetTextFont(62)
            cms_mark.SetTextSize(1.25 * font_scale)
            cms_mark.DrawLatex(mark_x1, mark_y1, "CMS")
            cms_mark.SetTextFont(52)
            cms_mark.SetTextSize(font_scale)
            cms_mark.DrawLatex(mark_x2, mark_y1, "Supplementary")
            cms_mark.SetTextFont(42)
            cms_mark.SetTextSize(font_scale)
            cms_mark.DrawLatex(mark_x1, mark_y2, "arXiv:2103.01290")
            cms_mark.SetTextFont(42)
            cms_mark.SetTextSize(font_scale)
            cms_mark.DrawLatex(mark_x1, mark_y3, self.region_labels[region])

            # Draw lumi stamp
            lumi = self.lumis[era]
            lumistamp = "{0:.1f} fb^{{-1}} (13 TeV)".format(lumi / 1000.0)
            cms_mark.SetTextAlign(31) # right aligned
            cms_mark.SetTextFont(42)
            cms_mark.SetTextSize(font_scale)
            cms_mark.DrawLatex(lumi_x, lumi_y, lumistamp)
            
            # pad for ratio
            pad = c.cd(2)
            # resize pad
            pad.SetGridy()
            # SetPad(xlow, ylow, xup, yup)
            pad.SetPad(0, 0, 1, lowerPadHeight)
            # set ticks on all sides of plot
            pad.SetTickx()
            pad.SetTicky()
            pad.SetLeftMargin(0.2)
            pad.SetRightMargin(0.1)
            pad.SetTopMargin(0.01)
            pad.SetBottomMargin(0.4)
            
            # --- draw --- 

            # First draw to setup axis, labels, etc.

            if doDataOverData:
                h_ratio_ZoverPhoton.SetMarkerStyle(ROOT.kFullDotLarge)
                h_ratio_ZoverPhoton.Draw(data_style)
            else:
                h_ratio_ZoverPhoton.Draw(draw_option)
            
            sim_rel_stat_unc.Draw("E2 same") # E2: error rectangles

            # Repeat first draw to get data points on top (over stat unc)
            
            if doDataOverData:
                h_ratio_ZoverPhoton.SetMarkerStyle(ROOT.kFullDotLarge)
                h_ratio_ZoverPhoton.Draw(data_style + " same")
            else:
                h_ratio_ZoverPhoton.Draw(draw_option + " same")
            
            if doFit:
                fit.Draw("same")
                # write chisq_r
                # give x, y coordinates (same as plot coordinates)
                #print "Fit: f(x) = (%.5f #pm %.5f) * x + (%.5f #pm %.5f)" % (p1, p1_err, p0, p0_err)
                mark = ROOT.TLatex()
                mark.SetTextSize(0.05)
                mark.DrawLatex(300.0, y_max - 0.2, "Fit: f(x) = %.5f + %.5f * x" % (p0, p1))
                mark.DrawLatex(300.0, y_max - 0.4, "#chi_{r}^{2} = %.3f" % chisq_r)
            
            if drawSyst:
                h_syst.Draw(draw_option + " same")
                # legend: TLegend(x1,y1,x2,y2)
                legend_x1 = 0.65
                legend_x2 = 0.90
                legend_y1 = 0.80
                legend_y2 = 0.95
                legend2 = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                legend2.SetFillStyle(0)
                legend2.SetBorderSize(0)
                legend2.SetLineWidth(1)
                legend2.SetNColumns(1)
                legend2.SetTextFont(42)
                legend2.AddEntry(h_ratio_ZoverPhoton,    ratio_title,           num_legend_style)
                if doFit:
                    legend2.AddEntry(fit,                "Fit to " + ratio_title,    "l")
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

            # ---------------------- #
            # --- study stat unc --- #
            # ---------------------- #
            
            selectionTag    = "_" + selection
            nameTag         = "_" + var
            h_map_norm      = self.getZHistoMap(region, nameTag, selectionTag, varLepton)
            h_map_shape     = self.S.getSimpleMap(region, nameTag, selectionTag, selectionTag, varPhoton)
            h_mc_map_Z      = self.getZMCHists(f, varLepton, h_map_norm)
            h_mc_map_photon = self.getPhotonMCHists(f, varPhoton, h_map_shape)
            mc_list_Z       = ["DY", "Diboson", "Rare", "TTbar", "SingleT"]
            if self.S.splitQCD:
                mc_list_photon = ["GJets", "QCD_Direct", "QCD_Fragmentation", "QCD_NonPrompt", "QCD_Fake", "WJets", "TTG", "tW", "Rare"]
            else:
                mc_list_photon = ["GJets", "QCD", "WJets", "TTG", "tW", "Rare"]
            
            colors = ["cherry red", "orange", "apple green", "cerulean", "bright purple", "fuchsia", "marigold", "lightish blue", "purpley blue", "terracotta"]
            
            c1 = ROOT.TCanvas("c1", "c1", 800, 800)
            

            
            # --- draw --- #
            c1.SetLogy(1) # set log y
            
            # legend: TLegend(x1,y1,x2,y2)
            #legend_x1 = 0.60
            #legend_x2 = 0.90
            #legend_y1 = 0.65
            #legend_y2 = 0.85
            #legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            #legend.SetFillStyle(0)
            #legend.SetBorderSize(0)
            #legend.SetLineWidth(1)
            #legend.SetNColumns(1)
            #legend.SetTextFont(42)
            
            # histograms to show stat. unc.
            h_mc_statunc_map_Z    = {}
            h_mc_relstatunc_map_Z = {}

            # use list to define order
            i = 0
            for histName in mc_list_Z:
                # load hist and create stat unc hist
                hist      = h_mc_map_Z[histName]
                # rebin 
                if rebin: 
                    h_new = hist.Rebin(self.n_bins, "h_new", self.xbins)
                else:
                    h_new = hist.Clone("h_new")
                # create stat unc hist
                h_statunc    = h_new.Clone("h_statunc")
                h_relstatunc = h_new.Clone("h_relstatunc")
                for b in range(1, h_new.GetNbinsX() + 1):
                    # statunc
                    h_statunc.SetBinContent(b, h_new.GetBinError(b))
                    h_statunc.SetBinError(b, 0.0)
                    # relstatunc
                    rel_stat_unc = 0.0
                    if h_new.GetBinContent(b) != 0.0:
                        rel_stat_unc = h_new.GetBinError(b) / h_new.GetBinContent(b)
                    h_relstatunc.SetBinContent(b, rel_stat_unc)
                    h_relstatunc.SetBinError(b, 0.0)
                h_mc_statunc_map_Z[histName]    = h_statunc 
                h_mc_relstatunc_map_Z[histName] = h_relstatunc 
                
                title   = "MC Yields in Z CR: {0} {1} {2}".format(self.labels[var], self.region_labels[region], era)
                x_title = self.labels[var]
                y_title = "Events"
                y_min = 0.01
                y_max = 10.0 ** 5
                setupHist(h_new,   title,  x_title,  y_title,  colors[i],   y_min,   y_max,   True,  3)
                h_new.GetXaxis().SetNdivisions(5, 5, 0, True)
                h_new.GetYaxis().SetNdivisions(5, 5, 0, True)
            
                # adding hists to legend causes seg fault at the moment
                #legend.AddEntry(h_new, histName, "l")
                
                if i == 0:
                    h_new.Draw("hist")
                else:
                    h_new.Draw("hist same")

                i += 1
            
            # adding hists to legend causes seg fault at the moment
            #legend.Draw()

            # save histograms
            if rebin:
                plot_name = "{0}{1}_Yield_Z_{2}_{3}_rebinned_{4}".format(self.plot_dir, fileTag, var, region, era)
            else:
                plot_name = "{0}{1}_Yield_Z_{2}_{3}_{4}".format(self.plot_dir, fileTag, var, region, era)
            
            c1.Update()
            c1.SaveAs(plot_name + ".pdf")
            c1.SaveAs(plot_name + ".png")
            
            
            # --- draw --- #
            c1.SetLogy(0) # unset log y
            
            # use list to define order
            i = 0
            for histName in mc_list_Z:
                hist = h_mc_statunc_map_Z[histName]
                
                title   = "MC Stat. Unc. in Z CR: {0} {1} {2}".format(self.labels[var], self.region_labels[region], era)
                x_title = self.labels[var]
                y_title = "Stat. Unc."
                y_min = 0.01
                if region == "LowDM":
                    y_max = 100.0
                else:
                    y_max = 10.0
                setupHist(hist,   title,  x_title,  y_title,  colors[i],   y_min,   y_max,   True,  3)
                hist.GetXaxis().SetNdivisions(5, 5, 0, True)
                hist.GetYaxis().SetNdivisions(5, 5, 0, True)
            
                if i == 0:
                    hist.Draw("hist")
                else:
                    hist.Draw("hist same")

                i += 1
            
            # save histograms
            if rebin:
                plot_name = "{0}{1}_StatUnc_Z_{2}_{3}_rebinned_{4}".format(self.plot_dir, fileTag, var, region, era)
            else:
                plot_name = "{0}{1}_StatUnc_Z_{2}_{3}_{4}".format(self.plot_dir, fileTag, var, region, era)
            
            c1.Update()
            c1.SaveAs(plot_name + ".pdf")
            c1.SaveAs(plot_name + ".png")
            
            # --- draw --- #
            c1.SetLogy(0) # unset log y
            
            # use list to define order
            i = 0
            for histName in mc_list_Z:
                hist = h_mc_relstatunc_map_Z[histName]
                
                title   = "MC Rel. Stat. Unc. in Z CR: {0} {1} {2}".format(self.labels[var], self.region_labels[region], era)
                x_title = self.labels[var]
                y_title = "Rel. Stat. Unc."
                y_min = 0.01
                y_max = 1.00
                setupHist(hist,   title,  x_title,  y_title,  colors[i],   y_min,   y_max,   True,  3)
                hist.GetXaxis().SetNdivisions(5, 5, 0, True)
                hist.GetYaxis().SetNdivisions(5, 5, 0, True)
            
                if i == 0:
                    hist.Draw("hist")
                else:
                    hist.Draw("hist same")

                i += 1
            
            # save histograms
            if rebin:
                plot_name = "{0}{1}_RelStatUnc_Z_{2}_{3}_rebinned_{4}".format(self.plot_dir, fileTag, var, region, era)
            else:
                plot_name = "{0}{1}_RelStatUnc_Z_{2}_{3}_{4}".format(self.plot_dir, fileTag, var, region, era)
            
            c1.Update()
            c1.SaveAs(plot_name + ".pdf")
            c1.SaveAs(plot_name + ".png")


        # write map to json file 
        if doDataOverData:
            with open (output_name, "w") as j:
                json.dump(val_map, j, sort_keys=True, indent=4, separators=(',', ' : '))


