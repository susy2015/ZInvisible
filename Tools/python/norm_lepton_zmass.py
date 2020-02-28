# norm_lepton_zmass.py

from colors import getColorIndex
import json
import os
import numpy as np
import ROOT
import tools
from tools import setupHist 

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Normalization:
    def __init__(self, plot_dir, verbose):
        self.plot_dir   = plot_dir
        self.verbose    = verbose
        self.ERROR_CODE = tools.ERROR_CODE
        self.norm_map       = {}
        self.norm_map_tex   = {}
        self.rz_syst_map    = {}
        self.histos         = {}
        self.output_file = 0 
        self.root_file   = 0
        self.eras = []
        self.systTag = ""
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "bestRecoZM"
        self.bin_types  = ["validation", "validationMetStudy", "search"]
        self.particles = ["Electron", "Muon"]
        self.channels = self.particles + ["Combined"]
        self.factors = ["R_Z", "R_T"]
        self.regions   = ["LowDM", "HighDM"]
        self.regions_tex = {
                             "LowDM"  : "Low $\Delta m$",
                             "HighDM" : "High $\Delta m$"
                           }
       
        self.bin_maps = {}
        with open("validation_bins_v3.json", "r") as j:
            self.bin_maps["validation"] = tools.stringifyMap(json.load(j))
        with open("validation_bins_metStudy.json", "r") as j:
            self.bin_maps["validationMetStudy"] = tools.stringifyMap(json.load(j))
        with open("search_bins_v4.json", "r") as j:
            self.bin_maps["search"] = tools.stringifyMap(json.load(j))
        
        # colors
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_black  = "black"
        
        # get selections from json file
        self.selections = {}
        self.selections_tex = {}
        for bin_type in self.bin_types:
            # selections per region
            self.selections[bin_type] = tools.getSelections(self.bin_maps[bin_type], bin_type, ["NJ"])
            self.selections_tex[bin_type] = {}
            for region in self.regions:
                self.selections_tex[bin_type][region] = {}
                for selection in self.selections[bin_type][region]:
                    self.selections_tex[bin_type][region][selection] = tools.getTexSelection(selection) 
                    if self.verbose:
                        print "For Normalization: {0} {1} {2}".format(bin_type, selection, self.selections_tex[bin_type][region][selection])
    
    def setupHistoMap(self, era):
        # histogram examples
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDatadata
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMZToLLstack
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMNoZToLLstack 
        # systematics
        # KEY: TH1D    DataMC_Muon_LowDM_Normalization_bestRecoZM_0to400_NBeq0_NSVeq0_pileup_syst_up_jetpt30bestRecoZMbestRecoZMZToLLstack;1   bestRecoZM
        # KEY: TH1D    DataMC_Muon_HighDM_Normalization_bestRecoZM_0to400_NBeq1_pileup_syst_down_jetpt30bestRecoZMbestRecoZMZToLLstack;1   bestRecoZM
        eraTag = "_" + era
        self.histos[era] = {}
        for bin_type in self.bin_types:
            self.histos[era][bin_type] = {}
            for particle in self.particles:
                self.histos[era][bin_type][particle] = {}
                for region in self.regions:
                    # using ZToLL and NoZToLL MC for normalization 
                    self.histos[era][bin_type][particle][region] = {}
                    for selection in self.selections[bin_type][region]: 
                        # apply syst. to MC only
                        dataSelectionTag    = "_" + selection + "_jetpt30"
                        mcSelectionTag      = "_" + selection + self.systTag + "_jetpt30"
                        # using Nb and Nsv selection
                        # testing
                        #print "In setupHistoMap(): " + "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + dataSelectionTag + 2 * self.variable + "Datadata"
                        #print "In setupHistoMap(): " + "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + mcSelectionTag   + 2 * self.variable + "ZToLLstack"
                        #print "In setupHistoMap(): " + "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + mcSelectionTag   + 2 * self.variable + "NoZToLLstack"
                        self.histos[era][bin_type][particle][region][selection] = { 
                            "Data"     : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + dataSelectionTag + 2 * self.variable + "Datadata",
                            "ZToLL"    : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + mcSelectionTag   + 2 * self.variable + "ZToLLstack",
                            "NoZToLL"  : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + mcSelectionTag   + 2 * self.variable + "NoZToLLstack",
                        }

    def calcNorm(self, A, x):
        # --- Solve for b --- #
        # x = A * b           #
        # b = A^(-1) * x      #
        # ------------------- #
        
        # inverse of matrix
        try:
            Ainverse = np.linalg.inv(A)
        except:
            print "ERROR: Cannot calculate inverse because matrix is singular (or possibly other error...)."
            print "A:"
            print A
            return [self.ERROR_CODE, self.ERROR_CODE]
        b = np.dot(Ainverse, x)
        if self.verbose:
            x_calc = np.dot(A, b)
            print "-----------------------------------"
            print "A:"
            print A
            print "Inverse of A:"
            print Ainverse
            print "x:"
            print x
            print "Calculated x (for validation):"
            print x_calc
            print "Difference between x and calculated x:"
            print x - x_calc
            print "b:"
            print b
            print "-----------------------------------"
    
        return b
    
    def calcError(self, A, A_error, x, x_error):
        # inverse of matrix
        try:
            Ainverse = np.linalg.inv(A)
        except:
            print "ERROR: Cannot calculate inverse because matrix is singular (or possibly other error...)."
            print "A:"
            print A
            return [self.ERROR_CODE, self.ERROR_CODE]
        Ainverse_error = tools.getMatrixInverseError(A, A_error)
        # from equation b = A^(-1) * x
        b_error = np.zeros(2)
        b1_error_1 = tools.getMultiplicationError(Ainverse[0][0] * x[0], Ainverse[0][0], Ainverse_error[0][0], x[0], x_error[0])
        b1_error_2 = tools.getMultiplicationError(Ainverse[0][1] * x[1], Ainverse[0][1], Ainverse_error[0][1], x[1], x_error[1])
        b2_error_1 = tools.getMultiplicationError(Ainverse[1][0] * x[0], Ainverse[1][0], Ainverse_error[1][0], x[0], x_error[0])
        b2_error_2 = tools.getMultiplicationError(Ainverse[1][1] * x[1], Ainverse[1][1], Ainverse_error[1][1], x[1], x_error[1])
        b_error[0] = tools.getAdditionError(b1_error_1, b1_error_2)
        b_error[1] = tools.getAdditionError(b2_error_1, b2_error_2)
        if self.verbose:
            print "-----------------------------------"
            print "A_error:"
            print A_error
            print "x_error:"
            print x_error
            print "b_error:"
            print b_error
            print "-----------------------------------"
        return b_error 
    
    
    def getNormAndError(self, file_name, era, systTag = ""):
        self.systTag = systTag
        self.eras.append(era)
        
        # define maps for every channel (not just particles)
        self.norm_map[era] = {}
        self.norm_map_tex[era] = {}
        for bin_type in self.bin_types:
            self.norm_map[era][bin_type] = {}
            self.norm_map_tex[era][bin_type] = {}
            for channel in self.channels:
                self.norm_map[era][bin_type][channel] = {}
                self.norm_map_tex[era][bin_type][channel] = {}
                for region in self.regions:
                    self.norm_map[era][bin_type][channel][region] = {}
                    self.norm_map_tex[era][bin_type][channel][region] = {}
                    for selection in self.selections[bin_type][region]:
                        self.norm_map[era][bin_type][channel][region][selection] = {}
                        self.norm_map_tex[era][bin_type][channel][region][selection] = {}
        
        if self.verbose:
            print "------------------------------------------------------------------------------"
            print "Era: {0}".format(era)
            print "File: {0}".format(file_name)
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        self.root_file = ROOT.TFile(file_name)
        # setup histogram map
        self.setupHistoMap(era)
        
        # calculate normalizations
        for bin_type in self.bin_types:
            for particle in self.particles:
                for region in self.regions:
                    for selection in self.selections[bin_type][region]:
                        self.runCalculation(era, bin_type, particle, region, selection)
                
        # weighted average: combine channels
        for bin_type in self.bin_types:
            for region in self.regions:
                for selection in self.selections[bin_type][region]:
                    for factor in self.factors:
                        self.calcCombined(factor, era, bin_type, region, selection)

        if self.verbose:
            print "------------------------------------------------------------------------------"
   

    def runCalculation(self, era, bin_type, particle, region, selection): 
        #WARNING: strings loaded from json file have type 'unicode'
        # ROOT cannot load histograms using unicode input: use type 'str'
        h_Data    = self.root_file.Get( str(self.variable + "/" + self.histos[era][bin_type][particle][region][selection]["Data"]    ) )
        h_ZToLL   = self.root_file.Get( str(self.variable + "/" + self.histos[era][bin_type][particle][region][selection]["ZToLL"]   ) )
        h_NoZToLL = self.root_file.Get( str(self.variable + "/" + self.histos[era][bin_type][particle][region][selection]["NoZToLL"] ) )
    
        if not h_Data:
            print "ERROR: Unable to load Data histogram {0}".format(self.variable + "/" + self.histos[era][bin_type][particle][region][selection]["Data"])

        #############################
        # Calculating Normalization #
        #############################
        minMass = 50.0
        Zmass = 91.0
        window = 10.0
        minZmass = Zmass - window
        maxZmass = Zmass + window
        bin_1 = h_Data.FindBin(minMass)
        bin_2 = h_Data.FindBin(minZmass) - 1
        bin_3 = bin_2 + 1
        bin_4 = h_Data.FindBin(maxZmass) - 1
        bin_5 = bin_4 + 1
        bin_6 = h_Data.GetNbinsX() + 1
        bins = [bin_1, bin_2, bin_3, bin_4, bin_5, bin_6]
        if self.verbose:
            for i, bin_i in enumerate(bins):
                lower_edge = h_Data.GetXaxis().GetBinLowEdge(bin_i)
                upper_edge = h_Data.GetXaxis().GetBinUpEdge(bin_i)
                print "bin_{0}: {1} [{2}, {3}]".format(i+1, bin_i, lower_edge, upper_edge)
        
        # pass error to IntegralAndError as Double_t & error 
        # use ROOT.Double for pass-by-ref of doubles
        # Double_t TH1::IntegralAndError  (Int_t binx1,
        #                                  Int_t binx2,
        #                                  Double_t & error,
        #                                  Option_t * option = "" 
        #                                 ) const
        
        # MC matrix (A)
        # a11: ZToLL on Z mass peak
        # a12: NoZToLL on Z mass peak
        # a21: ZToLL off Z mass peak
        # a22: NoZToLL off Z mass peak
        # Data (x)
        # x1: data on Z mass peak
        # x2: data off Z mass peak
        # Normalization (b)
        # b1: R_Z (ZToLL normalization)
        # b2: R_T (NoZToLL normalization)
        
        # Integration limis
        # on Z mass peak:  bin_3 to bin_4
        # off Z mass peak: bin_1 to bin_2 and bin_5 to bin_6 
    
        # MC matrix
        a11_error = ROOT.Double()
        a12_error = ROOT.Double()
        a21_error_1 = ROOT.Double()
        a21_error_2 = ROOT.Double()
        a22_error_1 = ROOT.Double()
        a22_error_2 = ROOT.Double()
        a11 = h_ZToLL.IntegralAndError(bin_3, bin_4, a11_error) 
        a12 = h_NoZToLL.IntegralAndError(bin_3, bin_4, a12_error) 
        a21 = h_ZToLL.IntegralAndError(bin_1, bin_2, a21_error_1) + h_ZToLL.IntegralAndError(bin_5, bin_6, a21_error_2)
        a22 = h_NoZToLL.IntegralAndError(bin_1, bin_2, a22_error_1) + h_NoZToLL.IntegralAndError(bin_5, bin_6, a22_error_2)
        A = np.array([[a11, a12], [a21, a22]])
        a21_error = tools.getAdditionError(a21_error_1, a21_error_2) 
        a22_error = tools.getAdditionError(a22_error_1, a22_error_2) 
        A_error = np.array([[a11_error, a12_error], [a21_error, a22_error]])
        # Data
        x1_error = ROOT.Double()
        x2_error_1 = ROOT.Double()
        x2_error_2 = ROOT.Double()
        x1 = h_Data.IntegralAndError(bin_3, bin_4, x1_error)
        x2 = h_Data.IntegralAndError(bin_1, bin_2, x2_error_1) + h_Data.IntegralAndError(bin_5, bin_6, x2_error_2)
        x = [x1, x2]
        x2_error = tools.getAdditionError(x2_error_1, x2_error_2)
        x_error = [x1_error, x2_error]
        # calculate normalization and error
        norm  = self.calcNorm(A, x)
        error = self.calcError(A, A_error, x, x_error)
        b1 = norm[0]
        b2 = norm[1]
        b1_error = error[0]
        b2_error = error[1]
        
        values = {}
        values["R_Z"] = b1
        values["R_T"] = b2
        values["R_Z_error"] = b1_error
        values["R_T_error"] = b2_error
        
        # loop over factors
        for factor in self.factors:
            value       = values[factor]
            value_error = values[factor + "_error"]
            # can be used to print if needed
            factor_print = "{0} = {1:.3f} +/- {2:.3f}".format(factor, value, value_error)
            factor_tex   = "${0:.3f} \pm {1:.3f}$".format(value, value_error)
            self.norm_map[era][bin_type][particle][region][selection][factor] = value
            self.norm_map[era][bin_type][particle][region][selection][factor + "_error"] = value_error
            self.norm_map_tex[era][bin_type][particle][region][selection][factor] = factor_tex

    def calcCombined(self, factor, era, bin_type, region, selection):
        # take the weighted average over particles
        weighted_average_numerator = 0.0
        weighted_average_denominator = 0.0
        errors_numerator = {}
        error_numerator = 0.0
        for particle in self.particles:
            value       = self.norm_map[era][bin_type][particle][region][selection][factor]
            value_error = self.norm_map[era][bin_type][particle][region][selection][factor + "_error"]
            weight = 1.0 / (value_error ** 2)
            weighted_average_numerator   += weight * value
            weighted_average_denominator += weight
            errors_numerator[particle] = tools.getConstantMultiplicationError(weight, value_error) 
        for i,particle in enumerate(self.particles):
            if i == 0:
                error_numerator = errors_numerator[particle]
            else:
                error_numerator = tools.getAdditionError(error_numerator, errors_numerator[particle])
        value       = weighted_average_numerator / weighted_average_denominator
        value_error = tools.getConstantMultiplicationError(1.0 / weighted_average_denominator, error_numerator)
        # can be used to print if needed
        factor_print = "{0} = {1:.3f} +/- {2:.3f}".format(factor, value, value_error)
        factor_tex   = "${0:.3f} \pm {1:.3f}$".format(value, value_error)
        self.norm_map[era][bin_type]["Combined"][region][selection][factor] = value
        self.norm_map[era][bin_type]["Combined"][region][selection][factor + "_error"] = value_error
        self.norm_map_tex[era][bin_type]["Combined"][region][selection][factor] = factor_tex

    def writeLine(self, line):
        self.output_file.write(line + "\n")

    # tex file with Rz values in all eras
    def makeTexFile(self, bin_type, output_name):
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
            
            # begin long table
            self.writeLine("\\centering")
            self.writeLine("\\begin{longtable}{|*6{p{0.16\\textwidth}|}}")
            self.writeLine("\hline Region & Selection & Era & CR & $R_Z$ & $R_T$ \\\\")
            # write values to table
            for region in self.regions:
                for selection in self.selections[bin_type][region]:
                    for era in self.eras:
                        for channel in self.channels:
                            R_Z = self.norm_map_tex[era][bin_type][channel][region][selection]["R_Z"] 
                            R_T = self.norm_map_tex[era][bin_type][channel][region][selection]["R_T"] 
                            region_tex     = self.regions_tex[region]
                            selections_tex = self.selections_tex[bin_type][region][selection]
                            era_tex = era.replace("_", " ")
                            self.writeLine("\hline {0} & {1} & {2} & {3} & {4} & {5} \\\\".format(region_tex, selections_tex, era_tex, channel, R_Z, R_T))
            self.writeLine("\hline")
            # end table
            self.writeLine("\\caption{Normalization}")
            self.writeLine("\end{longtable}")
            # end document
            self.writeLine("\end{document}")

    # Run 2 Rz values in table for analysis note 
    def makeTable(self, output_name, makeDoc=False):
        era      = "Run2"
        bin_type = "search"
        channelsForTable    = self.channels
        header              = "$\Nb$ & $\Nsv$ & $\Rz^{ee}$ & $\Rz^{\\mu\\mu}$ & $\\langle \Rz \\rangle$ & $\\langle \Rz \\rangle$ \\\\"
        header             += "\n& & & & & (w/ full unc.)\\\\"
        caption  = "Summary of the different regions used to derive the \Rz and $R_T$ factors."
        caption += "\nThe \Rz factors from the di-electron and di-muon control regions for the full Run 2 dataset are shown, as well as the weighted average $\\langle \Rz \\rangle$, all with statistical uncertainties."
        caption += "\nAn additional systematic uncertainty is obtained to account for differences in \Rz for different eras as shown in Figs.~\\ref{fig:norm_eras_lowdm}--\\ref{fig:norm_eras_highdm}, and the full uncertainty is listed in the last column."
        caption += "\nThe \Rz value obtained with $\Nb\geq2$ is used for search bins with $\Nb=2$, $\Nb\geq2$, $\Nb\geq3$."
        with open(output_name, "w+") as f:
            self.output_file = f
            
            if makeDoc:
                # begin document
                self.writeLine("\\documentclass{article}")
                self.writeLine("\\usepackage[utf8]{inputenc}")
                self.writeLine("\\usepackage{geometry}")
                self.writeLine("\\usepackage{longtable}")
                self.writeLine("\\usepackage{xspace}")
                self.writeLine("\\usepackage{amsmath}")
                self.writeLine("\\usepackage{graphicx}")
                self.writeLine("\\usepackage{cancel}")
                self.writeLine("\\usepackage{amsmath}")
                self.writeLine("\\geometry{margin=1in}")
                self.writeLine("\\input{VariableNames.tex}")
                self.writeLine("\\begin{document}")
            
            # begin table
            self.writeLine("\\begin{table}")
            self.writeLine("\\begin{center}")
            self.writeLine("\\caption{")
            self.writeLine(caption)
            self.writeLine("}")
            self.writeLine("\\label{tab:RZregions}")
            self.writeLine("\\begin{tabular}{%s}" % ( "c" * (3 + len(channelsForTable)) ) )
            self.writeLine(header)
            self.writeLine("\\hline")
            self.writeLine("\\multicolumn{2}{c}{low \dm normalization regions} \\\\")
            self.writeLine("\\hline")
            self.writeLine("0       & 0        & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["LowDM"]["NBeq0_NSVeq0"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["LowDM"]["NBeq0_NSVeq0"]["R_Z_total_unc"]) )
            self.writeLine("0       & $\geq$1  & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["LowDM"]["NBeq0_NSVge1"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["LowDM"]["NBeq0_NSVge1"]["R_Z_total_unc"]) )
            self.writeLine("1       & 0        & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["LowDM"]["NBeq1_NSVeq0"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["LowDM"]["NBeq1_NSVeq0"]["R_Z_total_unc"]) )
            self.writeLine("1       & $\geq$1  & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["LowDM"]["NBeq1_NSVge1"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["LowDM"]["NBeq1_NSVge1"]["R_Z_total_unc"]) )
            self.writeLine("$\geq$2 & --       & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["LowDM"]["NBge2"]["R_Z"]        for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["LowDM"]["NBge2"]["R_Z_total_unc"]       ) )
            self.writeLine("\\hline")
            self.writeLine("\\multicolumn{2}{c}{high \dm normalization regions} \\\\")
            self.writeLine("\\hline")
            # Nb = 2 is no longer used
            self.writeLine("1       & -- & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["HighDM"]["NBeq1"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["HighDM"]["NBeq1"]["R_Z_total_unc"])  ) 
            self.writeLine("$\geq$2 & -- & %s & %s \\\\" % (" & ".join(self.norm_map_tex[era][bin_type][channel]["HighDM"]["NBge2"]["R_Z"] for channel in channelsForTable), self.norm_map_tex[era][bin_type]["Combined"]["HighDM"]["NBge2"]["R_Z_total_unc"])  )
            self.writeLine("\\hline")
            # end table
            self.writeLine("\\end{tabular}")
            self.writeLine("\\end{center}")
            self.writeLine("\\end{table}")
            
            if makeDoc:
                # end document
                self.writeLine("\\end{document}")


    # make a plot of nomralization (y-axis) vs. era (x-axis) for different selections
    def makeComparison(self, bin_type):
        doFit = False
        draw_option = "hist error"
        total_era = "Run2"
        # treat Run2 era differentley
        eras = self.eras
        if self.eras[-1] == total_era:
            eras  = self.eras[:-1]
        nBins = len(eras)
        
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
        
        self.rz_syst_map[bin_type] = {}
        for region in self.regions:
            self.rz_syst_map[bin_type][region] = {}
            for selection in self.selections[bin_type][region]:
                region_tex              = self.regions_tex[region]
                selections_tex          = self.selections_tex[bin_type][region][selection]
                region_root_tex         = region_tex.replace("\\", "#")
                selections_root_tex     = selections_tex.replace("\\", "#")
                region_root_tex         = region_root_tex.replace("$", "")
                selections_root_tex     = selections_root_tex.replace("$", "")
                #print "{0} : {1}".format(region_tex, region_root_tex)
                #print "{0} : {1}".format(selections_tex, selections_root_tex)
                h_Electron       = ROOT.TH1F("h_Electron",        "h_Electron",        nBins, 0, nBins)
                h_Muon           = ROOT.TH1F("h_Muon",            "h_Muon",            nBins, 0, nBins)
                h_Combined       = ROOT.TH1F("h_Combined",        "h_Combined",        nBins, 0, nBins)
                h_Combined_Run2  = ROOT.TH1F("h_Combined_Run2",   "h_Combined_Run2",   nBins, 0, nBins)
                title = "Norm. for {0} bins, {1}, {2}".format(bin_type, region_root_tex, selections_root_tex)
                x_title = "Era" 
                y_title = "Norm. #left(R_{Z}#right)"
                y_min = -1.0
                y_max = 5.0
                #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                setupHist(h_Electron,       title, x_title, y_title, self.color_red,    y_min, y_max)
                setupHist(h_Muon,           title, x_title, y_title, self.color_blue,   y_min, y_max)
                setupHist(h_Combined,       title, x_title, y_title, self.color_black,  y_min, y_max)
                setupHist(h_Combined_Run2,  title, x_title, y_title, self.color_green,  y_min, y_max)
                Run2_norm     = self.norm_map["Run2"][bin_type]["Combined"][region][selection]["R_Z"]
                Run2_stat_err = self.norm_map["Run2"][bin_type]["Combined"][region][selection]["R_Z_error"]
                final_Run2_total_err  = -1
                for i in xrange(nBins):
                    era = eras[i]
                    # set bin contents and error
                    h_Electron.SetBinContent(           i + 1,      self.norm_map[era][bin_type]["Electron"][region][selection]["R_Z"])
                    h_Electron.SetBinError(             i + 1,      self.norm_map[era][bin_type]["Electron"][region][selection]["R_Z_error"])
                    h_Muon.SetBinContent(               i + 1,      self.norm_map[era][bin_type]["Muon"][region][selection]["R_Z"])
                    h_Muon.SetBinError(                 i + 1,      self.norm_map[era][bin_type]["Muon"][region][selection]["R_Z_error"])
                    h_Combined.SetBinContent(           i + 1,      self.norm_map[era][bin_type]["Combined"][region][selection]["R_Z"])
                    h_Combined.SetBinError(             i + 1,      self.norm_map[era][bin_type]["Combined"][region][selection]["R_Z_error"])
                    h_Combined_Run2.SetBinContent(      i + 1,      Run2_norm)
                    h_Combined_Run2.SetBinError(        i + 1,      self.norm_map["Run2"][bin_type]["Combined"][region][selection]["R_Z_error"])
                    # set bin labels
                    h_Electron.GetXaxis().SetBinLabel(      i + 1, era)
                    h_Muon.GetXaxis().SetBinLabel(          i + 1, era)
                    h_Combined.GetXaxis().SetBinLabel(      i + 1, era)
                    h_Combined_Run2.GetXaxis().SetBinLabel( i + 1, era)
                
                

                # loop over a range of Run 2 error
                # set Run 2 error
                # calculate chi sq
                # find Run 2 error such that reduced chi sq is 1.0
                
                # for Chi2Test()
                # "WW" = MC MC comparison (weighted-weighted)
                # "CHI2" = returns chi2 instead of p-value
                # nDegFree = nBins - 1 for chisq_r
                # nDegFree = nBins for chisq_Run2_r 
                
                # do a fit
                fit_value   = -1.0
                fit_error   = -1.0
                chisq_fit   = -1.0
                chisq_fit_r = -1.0
                if doFit:
                    f_Combined = ROOT.TF1("f1", "pol0", 0, 5)
                    h_Combined.Fit(f_Combined,  "", "", 0, 5)
                    f_Combined.SetLineColor(getColorIndex("violet"))
                    f_Combined.SetLineWidth(5)
                    fit_value = f_Combined.GetParameter(0)
                    fit_error = f_Combined.GetParError(0)
                    chisq_fit      = f_Combined.GetChisquare()
                    chisq_fit_r    = chisq_fit / (nBins - 1)
                
                chisq_Run2 = h_Combined_Run2.Chi2Test(h_Combined, "WW CHI2")
                chisq_Run2_r = chisq_Run2 / (nBins - 1)
                scale        = np.sqrt(chisq_Run2_r)
                if scale >= 1.0:
                    final_Run2_total_err = scale * Run2_stat_err
                else:
                    final_Run2_total_err = Run2_stat_err
                # set total Run 2 errror
                for i in xrange(nBins):
                    h_Combined_Run2.SetBinError(        i + 1,      final_Run2_total_err)
                mark = ROOT.TLatex()
                mark.SetTextSize(0.03)
                
                # title font size
                h_Electron.SetTitleSize(0.1)
                h_Muon.SetTitleSize(0.1)
                h_Combined.SetTitleSize(0.1)
                h_Combined_Run2.SetTitleSize(0.1)
                
                # draw
                h_Combined.Draw(draw_option)        
                h_Electron.Draw(draw_option + " same")        
                h_Muon.Draw(draw_option + " same")        
                if doFit:
                    f_Combined.Draw("same")
                h_Combined_Run2.Draw("same")

                # legend: TLegend(x1,y1,x2,y2)
                legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                legend.AddEntry(h_Electron,         "Electron",                  "l")
                legend.AddEntry(h_Muon,             "Muon",                      "l")
                legend.AddEntry(h_Combined,         "Combined e, #mu",           "l")
                legend.AddEntry(h_Combined_Run2,    "Run2",                      "l")
                if doFit:
                    legend.AddEntry(f_Combined,         "Fit to Combined e, #mu",    "l")
                legend.Draw()

                # write chisq
                # - give x, y coordinates (same as plot coordinates)
                # - y list is for positioning the text
                y_list = np.arange(y_max, 0.0, -0.3)
                if doFit: 
                    mark.DrawLatex(0.2, y_list[1], "Fit: f(x) = %.3f #pm %.3f"                % (fit_value, fit_error))
                    mark.DrawLatex(0.2, y_list[2], "Fit #chi_{r}^{2} = %.3f"                  % chisq_fit_r)
                    mark.DrawLatex(0.2, y_list[3], "Run 2 #chi_{r}^{2} = %.3f"                % chisq_Run2_r)
                    mark.DrawLatex(0.2, y_list[4], "S = %.3f"                                 % scale)
                    mark.DrawLatex(0.2, y_list[5], "R_{Z} #pm #sigma_{stat} = %.3f #pm %.3f"  % (Run2_norm, Run2_stat_err))
                    mark.DrawLatex(0.2, y_list[6], "Run 2 #sigma_{total} = %.3f"              % final_Run2_total_err)
                else:
                    mark.DrawLatex(0.2, y_list[1], "Run 2 #chi^{2} = %.3f"                    % chisq_Run2)
                    mark.DrawLatex(0.2, y_list[2], "Run 2 #chi_{r}^{2} = %.3f"                % chisq_Run2_r)
                    mark.DrawLatex(0.2, y_list[3], "S = %.3f"                                 % scale)
                    mark.DrawLatex(0.2, y_list[4], "#sigma_{stat} = %.3f"  %                  (Run2_stat_err))
                    mark.DrawLatex(0.2, y_list[5], "R_{Z} #pm #sigma_{total} = %.3f #pm %.3f"  % (Run2_norm, final_Run2_total_err))
                
                # save final Rz syst.
                # WARNING: this needs to be saved separately for validaiton and search bins, otherwise it will be overwritten
                self.rz_syst_map[bin_type][region][selection] = final_Run2_total_err
                
                # save final Rz syst. to norm_map
                # can be used to print if needed
                factor      = "R_Z_total_unc"
                value       = self.norm_map[total_era][bin_type]["Combined"][region][selection]["R_Z"] 
                value_error = final_Run2_total_err
                factor_print = "{0} = {1:.3f} +/- {2:.3f}".format(factor, value, value_error)
                factor_tex   = "${0:.3f} \pm {1:.3f}$".format(value, value_error)
                self.norm_map[total_era][bin_type]["Combined"][region][selection][factor] = value
                self.norm_map[total_era][bin_type]["Combined"][region][selection][factor + "_error"] = value_error
                self.norm_map_tex[total_era][bin_type]["Combined"][region][selection][factor] = factor_tex
                print "{0}".format(factor_print)

                # save histograms
                plot_name = "{0}Normalization_{1}_{2}_{3}".format(self.plot_dir, bin_type, region, selection)
                c.Update()
                c.SaveAs(plot_name + ".pdf")
                c.SaveAs(plot_name + ".png")
               
                # delete histograms to avoid memory leak
                del h_Electron
                del h_Muon
                del h_Combined
                del h_Combined_Run2



def main():
    json_file = "runs/run_2019-07-24.json"
    eras = ["2016", "2017", "2018_AB", "2018_CD"]
    #eras = ["2016", "2017", "2018_PreHEM", "2018_PostHEM"]
    latex_dir = "latex_files"
    # add "/" to directory if not present
    if latex_dir[-1] != "/":
        latex_dir += "/"
    verbose = False
    N = Normalization(verbose)
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        for era in eras:
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            N.getNormAndError(result_file, era)
    
        N.makeTexFile("validation", latex_dir + "validationBins_normalization_Zmass.tex")
        N.makeTexFile("search",     latex_dir + "searchBins_normalization_Zmass.tex")


if __name__ == "__main__":
    main()






