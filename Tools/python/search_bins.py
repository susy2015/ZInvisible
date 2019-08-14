# search_bins.py

import ROOT
from tools import setupHist, getMultiplicationErrorList

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# search bins and validation bins

# search bins 
class SearchBins:
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.low_dm_bins  = range( 0,  53)
        self.high_dm_bins = range(53, 204)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.binValues = {}
        self.eras = []

# vadliation bins
class ValidationBins:
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.binValues = {}
        self.values = ["norm", "shape", "mc", "pred"]
        # bins 0 to 45; bins 19, 20, and 21 are not included
        self.low_dm_bins_normal     = range( 0, 15)
        self.low_dm_bins_highmet    = range(15, 19)
        self.low_dm_bins            = range( 0, 19)
        self.high_dm_bins           = range(22, 46)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.low_dm_bins_met_250to300  = [4, 5, 8, 9, 10, 11]
        self.low_dm_bins_met_300toINF  = [17, 18]
        self.low_dm_bins_met_250to400  = [0, 1, 2, 3, 6, 7, 12, 13, 14]
        self.low_dm_bins_met_400toINF  = [15, 16]
        self.low_dm_bins_met_250toINF  = [15, 16, 17, 18]
        self.high_dm_bins_met_250to400 = range(22, 46, 2)
        self.high_dm_bins_met_400toINF = range(23, 46, 2)
        # colors
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_black  = "black"
        self.setupBins()

    def setupBins(self):
        self.bins = {
             0 : {"region" : "LowDM",  "selection" : "NBeq0_NSVeq0", "met" : "met_250to400"},
             1 : {"region" : "LowDM",  "selection" : "NBeq0_NSVeq0", "met" : "met_250to400"},
             2 : {"region" : "LowDM",  "selection" : "NBeq0_NSVge1", "met" : "met_250to400"},
             3 : {"region" : "LowDM",  "selection" : "NBeq0_NSVge1", "met" : "met_250to400"},
             4 : {"region" : "LowDM",  "selection" : "NBeq1_NSVeq0", "met" : "met_250to300"},
             5 : {"region" : "LowDM",  "selection" : "NBeq1_NSVeq0", "met" : "met_250to300"},
             6 : {"region" : "LowDM",  "selection" : "NBeq1_NSVeq0", "met" : "met_250to400"},
             7 : {"region" : "LowDM",  "selection" : "NBeq1_NSVeq0", "met" : "met_250to400"},
             8 : {"region" : "LowDM",  "selection" : "NBeq1_NSVge1", "met" : "met_250to300"},
             9 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to300"},
            10 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to300"},
            11 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to300"},
            12 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to400"},
            13 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to400"},
            14 : {"region" : "LowDM",  "selection" : "NBge2",        "met" : "met_250to400"},
            15 : {"region" : "LowDM",  "selection" : "NBeq0_NSVeq0", "met" : "met_250toINF"},
            16 : {"region" : "LowDM",  "selection" : "NBeq0_NSVge1", "met" : "met_250toINF"},
            17 : {"region" : "LowDM",  "selection" : "NBge1_NSVeq0", "met" : "met_250toINF"},
            18 : {"region" : "LowDM",  "selection" : "NBge1_NSVge1", "met" : "met_250toINF"},
            22 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            23 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            24 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            25 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
            26 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            27 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            28 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            29 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
            30 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            31 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            32 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            33 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            34 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            35 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            36 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_250to400"},
            37 : {"region" : "HighDM", "selection" : "NBeq1",        "met" : "met_400toINF"},
            38 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            39 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
            40 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            41 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
            42 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            43 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
            44 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_250to400"},
            45 : {"region" : "HighDM", "selection" : "NBge2",        "met" : "met_400toINF"},
        }

    def getValues(self, file_name, era):
        # new root file to save validation bin histograms
        new_file = "validationBinsZinv_" + era + ".root"
        draw_option = "hist error"
        eraTag = "_" + era
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]        = self.N.norm_map[era]["Combined"][region][selection]["R_Z"]
            self.binValues[era][b]["norm_error"]  = self.N.norm_map[era]["Combined"][region][selection]["R_Z_error"]
            self.binValues[era][b]["shape"]       = self.S.shape_map[era][region][met]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era][region][met + "_error"]

        # Z to NuNu MC histograms
        #TDirectoryFile*: nValidationBinLowDM_jetpt20
        #TH1D:  ZNuNu_nValidationBin_LowDM_jetpt20_2016nValidationBinLowDM_jetpt20nValidationBinLowDM_jetpt20ZJetsToNuNu Validation Bin Low DMdata
        #TDirectoryFile*: nValidationBinLowDMHighMET_jetpt20 
        #ZNuNu_nValidationBin_LowDM_HighMET_jetpt20_2016nValidationBinLowDMHighMET_jetpt20nValidationBinLowDMHighMET_jetpt20ZJetsToNuNu Validation Bin Low DM High METdata
        
        f = ROOT.TFile(file_name)
        h_validation_lowdm          = f.Get("nValidationBinLowDM_jetpt20/ZNuNu_nValidationBin_LowDM_jetpt20_" + era + "nValidationBinLowDM_jetpt20nValidationBinLowDM_jetpt20ZJetsToNuNu Validation Bin Low DMdata")
        h_validation_lowdm_highmet  = f.Get("nValidationBinLowDMHighMET_jetpt20/ZNuNu_nValidationBin_LowDM_HighMET_jetpt20_" + era + "nValidationBinLowDMHighMET_jetpt20nValidationBinLowDMHighMET_jetpt20ZJetsToNuNu Validation Bin Low DM High METdata")
        h_validation_highdm         = f.Get("nValidationBinHighDM_jetpt20/ZNuNu_nValidationBin_HighDM_jetpt20_" + era + "nValidationBinHighDM_jetpt20nValidationBinHighDM_jetpt20ZJetsToNuNu Validation Bin High DMdata")
        
        # Note: bin_i and b are different
        # bin_i is histogram bin number
        # b is validation bin number
        bin_i = 1
        for b in self.low_dm_bins_normal:
            self.binValues[era][b]["mc"]       = h_validation_lowdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_lowdm.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.low_dm_bins_highmet:
            self.binValues[era][b]["mc"]       = h_validation_lowdm_highmet.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_lowdm_highmet.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            self.binValues[era][b]["mc"]       = h_validation_highdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_highdm.GetBinError(bin_i)
            bin_i += 1

        f = ROOT.TFile(new_file, "recreate")
        # define histograms 
        h_mc_lowdm    = ROOT.TH1F("validation_mc_lowdm",    "validation_mc_lowdm",    19,  0, 19) 
        h_mc_highdm   = ROOT.TH1F("validation_mc_highdm",   "validation_mc_highdm",   24, 22, 46) 
        h_pred_lowdm  = ROOT.TH1F("validation_pred_lowdm",  "validation_pred_lowdm",  19,  0, 19) 
        h_pred_highdm = ROOT.TH1F("validation_pred_highdm", "validation_pred_highdm", 24, 22, 46) 

        # setup histograms
        #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
        setupHist(h_mc_lowdm,    "validation_mc_lowdm" + eraTag,    "Validation Bin", "Events", self.color_red,  10.0 ** -2, 10.0 ** 4)
        setupHist(h_mc_highdm,   "validation_mc_highdm" + eraTag,   "Validation Bin", "Events", self.color_red,  10.0 ** -2, 10.0 ** 4)
        setupHist(h_pred_lowdm,  "validation_pred_lowdm" + eraTag,  "Validation Bin", "Events", self.color_blue, 10.0 ** -2, 10.0 ** 4)
        setupHist(h_pred_highdm, "validation_pred_highdm" + eraTag, "Validation Bin", "Events", self.color_blue, 10.0 ** -2, 10.0 ** 4)

        # print values 
        if self.verbose:
            print era
        for b in self.all_bins:
            n       = self.binValues[era][b]["norm"]
            n_error = self.binValues[era][b]["norm_error"]
            s       = self.binValues[era][b]["shape"]
            s_error = self.binValues[era][b]["shape_error"]
            m       = self.binValues[era][b]["mc"]
            m_error = self.binValues[era][b]["mc_error"]
            p       = n * s * m
            x_list = [n, s, m]
            dx_list = [n_error, s_error, m_error]
            p_error = getMultiplicationErrorList(p, x_list, dx_list)
            self.binValues[era][b]["pred"] = p
            self.binValues[era][b]["pred_error"] = p_error
            
            for value in self.values:
                self.binValues[era][b][value + "_tex"] = "${0:.3f} \pm {1:.3f}$".format(self.binValues[era][b][value], self.binValues[era][b][value + "_error"])

            if self.verbose:
                print "bin {0}: N = {1:.3f} +/- {2:.3f} S = {3:.3f} +/- {4:.3f} M = {5:.3f} +/- {6:.3f} P = {7:.3f} +/- {8:.3f}".format(
                            b, n, n_error, s, s_error, m, m_error, p, p_error 
                        )
        # set histogram content and error
        bin_i = 1
        for b in self.low_dm_bins:
            h_mc_lowdm.SetBinContent(bin_i, self.binValues[era][b]["mc"])
            h_mc_lowdm.SetBinError(bin_i, self.binValues[era][b]["mc_error"])
            h_pred_lowdm.SetBinContent(bin_i, self.binValues[era][b]["pred"])
            h_pred_lowdm.SetBinError(bin_i, self.binValues[era][b]["pred_error"])
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            h_mc_highdm.SetBinContent(bin_i, self.binValues[era][b]["mc"])
            h_mc_highdm.SetBinError(bin_i, self.binValues[era][b]["mc_error"])
            h_pred_highdm.SetBinContent(bin_i, self.binValues[era][b]["pred"])
            h_pred_highdm.SetBinError(bin_i, self.binValues[era][b]["pred_error"])
            bin_i += 1

        h_map = {}
        h_map["lowdm"] = {}
        h_map["lowdm"]["mc"]   = h_mc_lowdm
        h_map["lowdm"]["pred"] = h_pred_lowdm
        h_map["highdm"] = {}
        h_map["highdm"]["mc"]   = h_mc_highdm
        h_map["highdm"]["pred"] = h_pred_highdm

        # draw histograms
        c = ROOT.TCanvas("c", "c", 800, 800)
        c.Divide(1, 2)
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.5
        legend_x2 = 0.9 
        legend_y1 = 0.7 
        legend_y2 = 0.9 

        ###################
        # Draw Histograms #
        ###################

        for region in h_map:
            h_mc   = h_map[region]["mc"]
            h_pred = h_map[region]["pred"]
            h_ratio = h_pred.Clone("h_ratio")
            h_ratio.Divide(h_mc)
        
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h_ratio, "h_ratio", "Validation Bin", "Events", self.color_black, 0.5, 1.5)

            # histograms
            c.cd(1)
            ROOT.gPad.SetLogy(1) # set log y
            # ZInv MC and Prediction
            h_mc.Draw(draw_option)
            h_pred.Draw("error same")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h_mc,   "MC",   "l")
            legend.AddEntry(h_pred, "Pred", "l")
            legend.Draw()
           
            # ratios
            c.cd(2)
            h_ratio.Draw(draw_option)
                
            # save histograms
            plot_name = self.plot_dir + "validation_" + region
            c.Update()
            c.SaveAs(plot_name + eraTag + ".pdf")
            c.SaveAs(plot_name + eraTag + ".png")
            # write histograms to file
            h_mc.Write()
            h_pred.Write()
        
        f.Close()
    
    def writeLine(self, line):
        self.output_file.write(line + "\n")

    def makeTexFile(self, output_name):
        # write latex file with table
        with open(output_name, "w+") as f:
            self.output_file = f
            # begin document
            self.writeLine("\documentclass{article}")
            self.writeLine("\usepackage[utf8]{inputenc}")
            self.writeLine("\usepackage{geometry}")
            self.writeLine("\usepackage{longtable}")
            self.writeLine("\geometry{margin=1in}")
            self.writeLine("")
            self.writeLine("\\begin{document}")
            for era in self.eras:
                era_tex = era.replace("_", " ")
                # begin table
                self.writeLine("\\begin{table}[ht]")
                self.writeLine("\\caption{{Z Invisible Predition for Validation Bins ({0})}}".format(era_tex))
                self.writeLine("\\vspace{2mm}")
                self.writeLine("\\centering")
                self.writeLine("\\begin{tabular}{|c|c|c|c|c|}")
                self.writeLine("\hline Bin & $R_Z$ & $S_\gamma$ & MC & Pred. \\\\")
                # write values to table
                for b in self.all_bins:
                    norm  = self.binValues[era][b]["norm_tex"]
                    shape = self.binValues[era][b]["shape_tex"]
                    mc    = self.binValues[era][b]["mc_tex"]
                    pred  = self.binValues[era][b]["pred_tex"]
                    self.writeLine("\hline {0} & {1} & {2} & {3} & {4} \\\\".format(b, norm, shape, mc, pred))
                self.writeLine("\hline")
                self.writeLine("\end{tabular}")
                # end table
                self.writeLine("\end{table}")
            # end document
            self.writeLine("\end{document}")



