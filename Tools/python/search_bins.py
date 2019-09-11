# search_bins.py
import json
import ROOT
from tools import setupHist, getMultiplicationErrorList, removeCuts

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# classes for search bins and validation bins

# Common class is parent class
# SearchBins and ValidationBins inherit from Common
class Common:
    def __init__(self):
        # colors
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_black  = "black"
    
    def writeLine(self, line):
        self.output_file.write(line + "\n")

    def makeTexFile(self, caption, output_name):
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
            self.writeLine("\\footnotesize")
            self.writeLine("\\tabcolsep=0.01cm")
            for era in self.eras:
                era_tex = era.replace("_", " ")
                # begin table
                self.writeLine("\\centering")
                # *n{} syntax with vertical lines for n columns; put last | in expression *5{|}
                self.writeLine("\\begin{longtable}{|*7{p{0.16\\textwidth}|}}")
                self.writeLine("\hline Bin & $R_{Z}$ & $S_{\gamma}$ & $N_{MC}$ & $N_{p}$ & $\\langle w \\rangle$ & $N_{eff}$ \\\\")
                # write values to table
                for b in self.all_bins:
                    norm  = self.binValues[era][b]["norm_tex"]
                    shape = self.binValues[era][b]["shape_tex"]
                    mc    = self.binValues[era][b]["mc_tex"]
                    pred  = self.binValues[era][b]["pred_tex"]
                    avg_w = self.binValues[era][b]["avg_w_tex"]
                    n_eff = self.binValues[era][b]["n_eff_tex"]
                    self.writeLine("\hline {0} & {1} & {2} & {3} & {4} & {5} & {6}\\\\".format(b, norm, shape, mc, pred, avg_w, n_eff))
                self.writeLine("\hline")
                # for longtable, caption must go at the bottom of the table... it is not working at the top
                self.writeLine("\\caption{{{0} ({1})}}".format(caption, era_tex))
                # end table
                self.writeLine("\end{longtable}")
            # end document
            self.writeLine("\end{document}")

    # ---------------------------------------------------------------------- #
    # makeHistos():                                                          #
    #    - make, plot, and save histograms                                   #
    #    - save relevant values to map                                       #
    # ---------------------------------------------------------------------- #
    def makeHistos(self, output_file, x_title, name, era):
        eraTag = "_" + era
        draw_option = "hist error"
        f_out = ROOT.TFile(output_file, "recreate")
        # define histograms 
        h_mc_lowdm    = ROOT.TH1F("mc_lowdm",    "mc_lowdm",    self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1) 
        h_mc_highdm   = ROOT.TH1F("mc_highdm",   "mc_highdm",   self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1) 
        h_pred_lowdm  = ROOT.TH1F("pred_lowdm",  "pred_lowdm",  self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1) 
        h_pred_highdm = ROOT.TH1F("pred_highdm", "pred_highdm", self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1) 

        # setup histograms
        #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
        setupHist(h_mc_lowdm,    "Z to Invisible MC and Prediction " + era, x_title, "Events", self.color_red,  10.0 ** -2, 10.0 ** 4)
        setupHist(h_mc_highdm,   "Z to Invisible MC and Prediction " + era, x_title, "Events", self.color_red,  10.0 ** -2, 10.0 ** 4)
        setupHist(h_pred_lowdm,  "Z to Invisible MC and Prediction " + era, x_title, "Events", self.color_blue, 10.0 ** -2, 10.0 ** 4)
        setupHist(h_pred_highdm, "Z to Invisible MC and Prediction " + era, x_title, "Events", self.color_blue, 10.0 ** -2, 10.0 ** 4)

        # save values in map
        if self.verbose:
            print era
        for b in self.all_bins:
            n       = self.binValues[era][b]["norm"]
            n_error = self.binValues[era][b]["norm_error"]
            s       = self.binValues[era][b]["shape"]
            s_error = self.binValues[era][b]["shape_error"]
            m       = self.binValues[era][b]["mc"]
            m_error = self.binValues[era][b]["mc_error"]
            
            # prediction:                   p     = bin value
            # uncertainty:                  sigma = bin error
            # average weight:               avg_w = sigma^2 / p 
            # effective number of events:   N_eff = p / avg_w
            
            p       = n * s * m
            x_list = [n, s, m]
            dx_list = [n_error, s_error, m_error]
            p_error = getMultiplicationErrorList(p, x_list, dx_list)
            avg_w = (p_error ** 2) / p
            if p == 0:
                print "ERROR: bin {0}, pred = {1}; seting avg weight to 1.0".format(b, p)
                avg_w   = 1.0
            else:
                avg_w   = (p_error ** 2) / p
            n_eff = p / avg_w
            n_eff_final = int(n_eff)
            if n_eff_final == 0:
                print "ERROR: bin {0}, n_eff_final = {1}; leaving avg weight unchanged".format(b, n_eff_final)
                avg_w_final = avg_w
            else:
                avg_w_final = p / n_eff_final

            self.binValues[era][b]["pred"]          = p
            self.binValues[era][b]["pred_error"]    = p_error
            self.binValues[era][b]["avg_w"]         = avg_w_final
            self.binValues[era][b]["n_eff"]         = n_eff_final
            
            for value in self.values:
                self.binValues[era][b][value + "_tex"] = "${0:.3f} \pm {1:.3f}$".format(self.binValues[era][b][value], self.binValues[era][b][value + "_error"])
                
            for value in ["avg_w", "n_eff"]:
                self.binValues[era][b][value + "_tex"] = "${0:.3f}$".format(self.binValues[era][b][value])

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
            setupHist(h_ratio, "Z to Invisible Prediction / MC", x_title, "Pred / MC", self.color_black, 0.5, 1.5)

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
            plot_name = self.plot_dir + name + "_" + region
            c.Update()
            c.SaveAs(plot_name + eraTag + ".pdf")
            c.SaveAs(plot_name + eraTag + ".png")
            # write histograms to file
            h_mc.Write()
            h_pred.Write()
        
        f_out.Close()


# search bins 
class SearchBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.low_dm_start   = 0
        self.low_dm_end     = 52
        self.high_dm_start  = 53
        self.high_dm_end    = 203
        self.low_dm_nbins   = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins  = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins    = list(str(b) for b in range( self.low_dm_start,  self.low_dm_end + 1)) 
        self.high_dm_bins   = list(str(b) for b in range( self.high_dm_start, self.high_dm_end + 1)) 
        self.all_bins       = self.low_dm_bins + self.high_dm_bins
        self.binValues = {}
        self.values = ["norm", "shape", "mc", "pred"]
        with open("search_bins.json", "r") as j:
            self.bins = json.load(j)
    
    def getValues(self, file_name, era):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            # remove cuts from selection for norm and shape
            selection_norm  = removeCuts(selection, "NJ")
            selection_shape = removeCuts(selection, "NSV")
            if self.verbose:
                print "{0}: {1} {2} {3} {4}".format(b, region, selection_norm, selection_shape, met)
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]        = self.N.norm_map[era]["search"]["Combined"][region][selection_norm]["R_Z"]
            self.binValues[era][b]["norm_error"]  = self.N.norm_map[era]["search"]["Combined"][region][selection_norm]["R_Z_error"]
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["search"][region][selection_shape][met]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["search"][region][selection_shape][met + "_error"]
        
        # Z to NuNu MC histograms
        f_in        = ROOT.TFile(file_name, "read")
        h_lowdm     = f_in.Get("nSearchBinLowDM_jetpt20/ZNuNu_nSearchBin_LowDM_jetpt20_" + era + "nSearchBinLowDM_jetpt20nSearchBinLowDM_jetpt20ZJetsToNuNu Search Bin Low DMdata")
        h_highdm    = f_in.Get("nSearchBinHighDM_jetpt20/ZNuNu_nSearchBin_HighDM_jetpt20_" + era + "nSearchBinHighDM_jetpt20nSearchBinHighDM_jetpt20ZJetsToNuNu Search Bin High DMdata")
        # Note: bin_i and b are different
        # bin_i is histogram bin number
        # b is search bin number
        bin_i = 1
        for b in self.low_dm_bins:
            self.binValues[era][b]["mc"]       = h_lowdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_lowdm.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            self.binValues[era][b]["mc"]       = h_highdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_highdm.GetBinError(bin_i)
            bin_i += 1

        # new root file to save search bin histograms
        new_file = "searchBinsZinv_" + era + ".root"
        self.makeHistos(new_file, "Search Bin", "search", era)
        f_in.Close()


# vadliation bins
class ValidationBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.binValues = {}
        self.values = ["norm", "shape", "mc", "pred"]
        # bins 0 to 45; bins 19, 20, and 21 are not included
        # old version
        #self.low_dm_bins_normal     = list(str(b) for b in range( 0, 15)) 
        #self.low_dm_bins_highmet    = list(str(b) for b in range(15, 19))
        #self.low_dm_bins            = list(str(b) for b in range( 0, 19))
        #self.high_dm_bins           = list(str(b) for b in range(22, 46))
        #self.all_bins               = self.low_dm_bins + self.high_dm_bins
        # new version
        self.low_dm_start           = 0
        self.low_dm_end             = 18
        self.high_dm_start          = 22
        self.high_dm_end            = 45
        self.low_dm_nbins           = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins          = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins            = list(str(b) for b in range( self.low_dm_start,  self.low_dm_end + 1)) 
        self.high_dm_bins           = list(str(b) for b in range( self.high_dm_start, self.high_dm_end + 1)) 
        self.low_dm_bins_normal     = list(str(b) for b in range( 0, 15)) 
        self.low_dm_bins_highmet    = list(str(b) for b in range(15, 19))
        self.all_bins               = self.low_dm_bins + self.high_dm_bins
        with open("validation_bins.json", "r") as j:
            self.bins = json.load(j)

    def getValues(self, file_name, era):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            # remove cuts from selection for norm and shape
            selection_norm  = removeCuts(selection, "NJ")
            selection_shape = removeCuts(selection, "NSV")
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]        = self.N.norm_map[era]["validation"]["Combined"][region][selection_norm]["R_Z"]
            self.binValues[era][b]["norm_error"]  = self.N.norm_map[era]["validation"]["Combined"][region][selection_norm]["R_Z_error"]
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["validation"][region][selection_shape][met]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["validation"][region][selection_shape][met + "_error"]

        # Z to NuNu MC histograms
        #TDirectoryFile*: nValidationBinLowDM_jetpt20
        #TH1D:  ZNuNu_nValidationBin_LowDM_jetpt20_2016nValidationBinLowDM_jetpt20nValidationBinLowDM_jetpt20ZJetsToNuNu Validation Bin Low DMdata
        #TDirectoryFile*: nValidationBinLowDMHighMET_jetpt20 
        #ZNuNu_nValidationBin_LowDM_HighMET_jetpt20_2016nValidationBinLowDMHighMET_jetpt20nValidationBinLowDMHighMET_jetpt20ZJetsToNuNu Validation Bin Low DM High METdata
        f_in                = ROOT.TFile(file_name, "read")
        h_lowdm             = f_in.Get("nValidationBinLowDM_jetpt20/ZNuNu_nValidationBin_LowDM_jetpt20_" + era + "nValidationBinLowDM_jetpt20nValidationBinLowDM_jetpt20ZJetsToNuNu Validation Bin Low DMdata")
        h_lowdm_highmet     = f_in.Get("nValidationBinLowDMHighMET_jetpt20/ZNuNu_nValidationBin_LowDM_HighMET_jetpt20_" + era + "nValidationBinLowDMHighMET_jetpt20nValidationBinLowDMHighMET_jetpt20ZJetsToNuNu Validation Bin Low DM High METdata")
        h_highdm            = f_in.Get("nValidationBinHighDM_jetpt20/ZNuNu_nValidationBin_HighDM_jetpt20_" + era + "nValidationBinHighDM_jetpt20nValidationBinHighDM_jetpt20ZJetsToNuNu Validation Bin High DMdata")
        
        # Note: bin_i and b are different
        # bin_i is histogram bin number
        # b is validation bin number
        bin_i = 1
        for b in self.low_dm_bins_normal:
            self.binValues[era][b]["mc"]       = h_lowdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_lowdm.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.low_dm_bins_highmet:
            self.binValues[era][b]["mc"]       = h_lowdm_highmet.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_lowdm_highmet.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            self.binValues[era][b]["mc"]       = h_highdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_highdm.GetBinError(bin_i)
            bin_i += 1

        # new root file to save validation bin histograms
        new_file = "validationBinsZinv_" + era + ".root"
        self.makeHistos(new_file, "Validation Bin", "validation", era)
        f_in.Close()
  



