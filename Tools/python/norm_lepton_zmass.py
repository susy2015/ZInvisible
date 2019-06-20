# norm_lepton_zmass.py

import os
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Normalization:
    def __init__(self, useAllMC, verbose):
        self.useAllMC = useAllMC
        self.verbose = verbose
        self.ERROR_CODE = -999
        self.norm_map = {}
        self.norm_map_tex = {}
        self.histos = {}
        self.output_file = 0 
        self.eras = []
        self.particles = ["Electron", "Muon"]
        self.channels = self.particles + ["Combined"]
        self.regions   = ["LowDM", "HighDM"]
        self.regions_tex = {
                             "LowDM"  : "Low $\Delta m$",
                             "HighDM" : "High $\Delta m$"
                           }
        self.factors = ["R_Z", "R_T"]
        # variable is also TDirectoryFile that holds histograms 
        self.variable = "bestRecoZM"
    
    def setupHistoMap(self, year):
        # histograms
        # examples (DY and ttbar only)
        # DataMC_Electron_LowDM_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDatadata
        # DataMC_Electron_LowDM_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDYstack
        # DataMC_Electron_LowDM_bestRecoZM_0to400_2016bestRecoZMbestRecoZMt#bar{t}stack 
        # examples (all MC)
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDatadata
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMZToLLstack
        # DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_2016bestRecoZMbestRecoZMNoZToLLstack 
        self.histos[year] = {}
        for particle in self.particles:
            self.histos[year][particle] = {}
            for region in self.regions:
                if self.useAllMC:
                    # using ZToLL and NoZToLL MC for normalization 
                    self.histos[year][particle][region] = { 
                        "Data"     : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMDatadata",
                        "ZToLL"    : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMZToLLstack",
                        "NoZToLL"  : "DataMC_" + particle + "_" + region + "_Normalization_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMNoZToLLstack",
                    }
                else:
                    # using only DY and ttbar for normalization 
                    self.histos[year][particle][region] = { 
                        "Data"     : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMDatadata",
                        "ZToLL"    : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMDYstack",
                        "NoZToLL"  : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMt#bar{t}stack",
                    }

    def setUseAllMC(self, value):
        self.useAllMC = value

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
        Ainverse_error = self.getMatrixInverseError(A, A_error)
        # from equation b = A^(-1) * x
        b_error = np.zeros(2)
        b1_error_1 = self.getMultiplicationError(Ainverse[0][0] * x[0], Ainverse[0][0], Ainverse_error[0][0], x[0], x_error[0])
        b1_error_2 = self.getMultiplicationError(Ainverse[0][1] * x[1], Ainverse[0][1], Ainverse_error[0][1], x[1], x_error[1])
        b2_error_1 = self.getMultiplicationError(Ainverse[1][0] * x[0], Ainverse[1][0], Ainverse_error[1][0], x[0], x_error[0])
        b2_error_2 = self.getMultiplicationError(Ainverse[1][1] * x[1], Ainverse[1][1], Ainverse_error[1][1], x[1], x_error[1])
        b_error[0] = self.getAdditionError(b1_error_1, b1_error_2)
        b_error[1] = self.getAdditionError(b2_error_1, b2_error_2)
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
    
    def getAdditionError(self, dx, dy):
        # q  = x + y
        # q  = x - y
        # dq = sqrt( dx^2 + dy^2 )
        return abs(np.sqrt( dx**2 + dy**2 ))
    
    def getConstantMultiplicationError(self, a, dx):
        # a is a constant
        # q  = a * x
        # dq = |a| * dx
        return abs(a * dx)
    
    def getMultiplicationError(self, q, x, dx, y, dy):
        # q = x * y 
        # q = x / y 
        # dq = q * sqrt( (dx/x)^2 + (dy/y)^2) )
        return abs(q * np.sqrt( (dx/x)**2 + (dy/y)**2 ))
    
    def getMultiplicationErrorList(self, q, x_list, dx_list):
        # q = x * y * ... 
        # q = x / y / ...
        # dq = q * sqrt( (dx/x)^2 + (dy/y)^2) + ... )
        if len(x_list) != len(dx_list):
            print "ERROR in getMultiplicationErrorList(): x_list and dx_list do not have the same length"
            return self.ERROR_CODE
        s = 0.0
        for i in xrange(len(x_list)):
            s += (dx_list[i] / x_list[i]) ** 2
        return abs(q * np.sqrt(s))
    
    def getMatrixDeterminantError(self, A, A_error):
        # get error for determinant
        # |A| = a11 * a22 - a12 * a21
        error_1 = self.getMultiplicationError(A[0][0] * A[1][1], A[0][0], A_error[0][0], A[1][1], A_error[1][1])
        error_2 = self.getMultiplicationError(A[0][1] * A[1][0], A[0][1], A_error[0][1], A[1][0], A_error[1][0])
        det_error = self.getAdditionError(error_1, error_2)
        return det_error
    
    def getMatrixInverseError(self, A, A_error):
        # from equation for A^(-1)
        det = np.linalg.det(A)
        det_error = self.getMatrixDeterminantError(A, A_error)
        Ainverse_error = np.zeros((2,2))
        Ainverse_error[0][0] = self.getMultiplicationError( A[1][1] / det,  A[1][1], A_error[1][1], det, det_error)
        Ainverse_error[0][1] = self.getMultiplicationError(-A[0][1] / det, -A[0][1], A_error[0][1], det, det_error)
        Ainverse_error[1][0] = self.getMultiplicationError(-A[1][0] / det, -A[1][0], A_error[1][0], det, det_error)
        Ainverse_error[1][1] = self.getMultiplicationError( A[0][0] / det,  A[0][0], A_error[0][0], det, det_error)
        return Ainverse_error
    
    def getNormAndError(self, file_name, era):
        # currently the histograms are named by year (2018) and not era (2018_AB)
        # we should probalby change the histograms to use era (2018_AB)
        year = era[0:4]
        self.eras.append(era)
        
        # define maps for every channel (not just particles)
        self.norm_map[era] = {}
        self.norm_map_tex[era] = {}
        for channel in self.channels:
            self.norm_map[era][channel] = {}
            self.norm_map_tex[era][channel] = {}
            for region in self.regions:
                self.norm_map[era][channel][region] = {}
                self.norm_map_tex[era][channel][region] = {}
        
        
        if self.verbose:
            print "------------------------------------------------------------------------------"
            print "Era: {0}".format(era)
            print "Year: {0}".format(year)
            print "File: {0}".format(file_name)
        # check that the file exists
        if not os.path.isfile(file_name): 
            print "The file {0} does not exist".format(file_name)
            return
        
        f = ROOT.TFile(file_name)
        # setup histogram map
        self.setupHistoMap(year)
        
        for particle in self.particles:
            if self.verbose:
                print particle
            for region in self.regions:
                if self.verbose:
                    print region
                h_Data    = f.Get(self.variable + "/" + self.histos[year][particle][region]["Data"])
                h_ZToLL   = f.Get(self.variable + "/" + self.histos[year][particle][region]["ZToLL"])
                h_NoZToLL = f.Get(self.variable + "/" + self.histos[year][particle][region]["NoZToLL"])
    
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
                a21_error = self.getAdditionError(a21_error_1, a21_error_2) 
                a22_error = self.getAdditionError(a22_error_1, a22_error_2) 
                A_error = np.array([[a11_error, a12_error], [a21_error, a22_error]])
                # Data
                x1_error = ROOT.Double()
                x2_error_1 = ROOT.Double()
                x2_error_2 = ROOT.Double()
                x1 = h_Data.IntegralAndError(bin_3, bin_4, x1_error)
                x2 = h_Data.IntegralAndError(bin_1, bin_2, x2_error_1) + h_Data.IntegralAndError(bin_5, bin_6, x2_error_2)
                x = [x1, x2]
                x2_error = self.getAdditionError(x2_error_1, x2_error_2)
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
                    factor_print = "{0} = {1:.3f} +/- {2:.3f}".format(factor, value, value_error)
                    factor_tex   = "${0:.3f} \pm {1:.3f}$".format(value, value_error)
                    if self.verbose:
                        print factor_print 
                    self.norm_map[era][particle][region][factor] = value
                    self.norm_map[era][particle][region][factor + "_error"] = value_error
                    self.norm_map_tex[era][particle][region][factor] = factor_tex
                

        # weighted average: combine channels
        if self.verbose:
            print "Combined"
        for region in self.regions:
            if self.verbose:
                print region
            for factor in self.factors:
                # take the weighted average over particles
                weighted_average_numerator = 0.0
                weighted_average_denominator = 0.0
                errors_numerator = {}
                error_numerator = 0.0
                for particle in self.particles:
                    value       = self.norm_map[era][particle][region][factor]
                    value_error = self.norm_map[era][particle][region][factor + "_error"]
                    weight = 1.0 / (value_error ** 2)
                    weighted_average_numerator   += weight * value
                    weighted_average_denominator += weight
                    errors_numerator[particle] = self.getConstantMultiplicationError(weight, value_error) 
                for i,particle in enumerate(self.particles):
                    if i == 0:
                        error_numerator = errors_numerator[particle]
                    else:
                        error_numerator = self.getAdditionError(error_numerator, errors_numerator[particle])
                value       = weighted_average_numerator / weighted_average_denominator
                value_error = self.getConstantMultiplicationError(1.0 / weighted_average_denominator, error_numerator)
                factor_print = "{0} = {1:.3f} +/- {2:.3f}".format(factor, value, value_error)
                factor_tex   = "${0:.3f} \pm {1:.3f}$".format(value, value_error)
                if self.verbose:
                    print factor_print 
                self.norm_map[era]["Combined"][region][factor] = value
                self.norm_map[era]["Combined"][region][factor + "_error"] = value_error
                self.norm_map_tex[era]["Combined"][region][factor] = factor_tex

        if self.verbose:
            print "------------------------------------------------------------------------------"
    
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
            self.writeLine("Normalization")
            # begin table
            self.writeLine("\\begin{table}[h]")
            self.writeLine("\\begin{tabular}{|c|c|c|c|c|}")
            self.writeLine("\hline Selection & Era & CR & $R_Z$ & $R_T$ \\\\")
            # write values to table
            for region in self.regions:
                for era in self.eras:
                    for channel in self.channels:
                        R_Z = self.norm_map_tex[era][channel][region]["R_Z"] 
                        R_T = self.norm_map_tex[era][channel][region]["R_T"] 
                        region_tex = self.regions_tex[region]
                        era_tex = era.replace("_", " ")
                        self.writeLine("\hline {0} & {1} & {2} & {3} & {4} \\\\".format(region_tex, era_tex, channel, R_Z, R_T))
            self.writeLine("\hline")
            self.writeLine("\end{tabular}")
            # end table
            self.writeLine("\end{table}")
            # end document
            self.writeLine("\end{document}")

def main():
    verbose = False
    # set per version
    version = 3
    useAllMC_map = {1:False, 2:False, 3:True}
    useAllMC = useAllMC_map[version]
    N = Normalization(useAllMC, verbose)
    if version == 1:
        # May 5, 2019 Results
        N.getNormAndError("condor/DataMC_2016_submission_2019-05-05_21-57-41/result.root", "2016")
        N.getNormAndError("condor/DataMC_2017_submission_2019-05-05_22-28-09/result.root", "2017")
        N.getNormAndError("condor/DataMC_2018_submission_2019-05-05_22-44-26/result.root", "2018")
    elif version == 2:
        # May 9, 2019 Results
        N.getNormAndError("condor/DataMC_2016_submission_2019-05-09_17-19-42/result.root", "2016")
        N.getNormAndError("condor/DataMC_2017_submission_2019-05-09_17-16-54/result.root", "2017")
        N.getNormAndError("condor/DataMC_2018_submission_2019-05-09_17-15-04/result.root", "2018")
    elif version == 3:
        # May 16, 2019 Results
        N.getNormAndError("condor/DataMC_2016_submission_2019-05-16_10-06-59/result.root",    "2016")
        N.getNormAndError("condor/DataMC_2017_submission_2019-05-16_10-09-29/result.root",    "2017")
        N.getNormAndError("condor/DataMC_2018_AB_submission_2019-05-16_10-10-30/result.root", "2018_AB")
        N.getNormAndError("condor/DataMC_2018_CD_submission_2019-05-16_10-12-04/result.root", "2018_CD")
    N.makeTexFile("normalization.tex")


if __name__ == "__main__":
    main()


