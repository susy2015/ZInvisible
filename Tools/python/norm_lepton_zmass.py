# norm_lepton_zmass.py

import json
import os
import numpy as np
import ROOT
import tools

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Normalization:
    def __init__(self, verbose):
        self.verbose = verbose
        self.ERROR_CODE = tools.ERROR_CODE
        self.norm_map = {}
        self.norm_map_tex = {}
        self.histos = {}
        self.output_file = 0 
        self.root_file   = 0
        self.eras = []
        
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "bestRecoZM"
        self.bin_types  = ["validation", "search"]
        self.particles = ["Electron", "Muon"]
        self.channels = self.particles + ["Combined"]
        self.factors = ["R_Z", "R_T"]
        self.regions   = ["LowDM", "HighDM"]
        self.regions_tex = {
                             "LowDM"  : "Low $\Delta m$",
                             "HighDM" : "High $\Delta m$"
                           }
       
        
        self.bin_maps = {}
        with open("validation_bins.json", "r") as j:
            self.bin_maps["validation"] = json.load(j)
        with open("search_bins.json", "r") as j:
            self.bin_maps["search"] = json.load(j)
        
        # get selections from json file
        self.selections = {}
        self.selections_tex = {}
        for bin_type in self.bin_types:
            # selections per region
            self.selections[bin_type] = tools.getSelections(self.bin_maps[bin_type], bin_type, "NJ")
            self.selections_tex[bin_type] = {}
            for region in self.regions:
                self.selections_tex[bin_type][region] = {}
                for selection in self.selections[bin_type][region]:
                    self.selections_tex[bin_type][region][selection] = tools.getTexSelection(selection) 
                    print "For Normalization: {0} {1} {2}".format(bin_type, selection, self.selections_tex[bin_type][region][selection])
         

        #TODO: delete once json file method works
        #self.selections  = {
        #                      "LowDM"  : ["NBeq0_NSVeq0", "NBeq0_NSVge1", "NBeq1_NSVeq0", "NBeq1_NSVge1", "NBge1_NSVeq0", "NBge1_NSVge1", "NBge2"],
        #                      "HighDM" : ["NBeq1", "NBge2"]
        #                   }
        #self.selections_tex = {
        #                         "NBeq0_NSVeq0" : "$N_{b} = 0, N_{sv} = 0$",
        #                         "NBeq0_NSVge1" : "$N_{b} = 0, N_{sv} \geq 1$",
        #                         "NBeq1_NSVeq0" : "$N_{b} = 1, N_{sv} = 0$",
        #                         "NBeq1_NSVge1" : "$N_{b} = 1, N_{sv} \geq 1$",
        #                         "NBge1_NSVeq0" : "$N_{b} \geq 1, N_{sv} = 0$",
        #                         "NBge1_NSVge1" : "$N_{b} \geq 1, N_{sv} \geq 1$",
        #                         "NBeq1"        : "$N_{b} = 1$",
        #                         "NBge2"        : "$N_{b} \geq 2$"
        #                      }
    
    def setupHistoMap(self, era):
        # histogram examples
        # DataMC_Electron_LowDM_bestRecoZM_50to250_NBeq0_NSVeq0_jetpt20_2016bestRecoZMbestRecoZMDatadata
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDatadata
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMZToLLstack
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMNoZToLLstack 
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
                        selectionTag = "_" + selection + "_jetpt20"
                        # using Nb and Nsv selection
                        self.histos[era][bin_type][particle][region][selection] = { 
                            "Data"     : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + selectionTag + eraTag + 2 * self.variable + "Datadata",
                            "ZToLL"    : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + selectionTag + eraTag + 2 * self.variable + "ZToLLstack",
                            "NoZToLL"  : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400" + selectionTag + eraTag + 2 * self.variable + "NoZToLLstack",
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
    
    
    def getNormAndError(self, file_name, era):
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
        print "File: {0}".format(file_name)
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
            
            # begin table
            #self.writeLine("\\begin{table}[ht]")
            #self.writeLine("\\caption{Normalization}")
            #self.writeLine("\\vspace{2mm}")
            #self.writeLine("\\centering")
            
            # begin long table
            self.writeLine("\\centering")
            #self.writeLine("\\begin{tabular}{|c|c|c|c|c|c|}")
            self.writeLine("\\begin{longtable}{| p{0.16\\textwidth} | p{0.16\\textwidth} | p{0.16\\textwidth} | p{0.16\\textwidth} | p{0.16\\textwidth} | p{0.16\\textwidth} |}")
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
            #self.writeLine("\end{tabular}")
            #self.writeLine("\end{table}")
            self.writeLine("\\caption{Normalization}")
            self.writeLine("\end{longtable}")
            # end document
            self.writeLine("\end{document}")

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






