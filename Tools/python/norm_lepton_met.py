# norm_lepton_met.py

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
        self.variable = "metWithLL"
    
    def setupHistoMap(self, year):
       # DataMC_Electron_LowDM_Normalization_met_2016metWithLLmetWithLLZToLLstack
       # DataMC_Electron_LowDM_Normalization_met_2016metWithLLmetWithLLNoZToLLstack
       # DataMC_Electron_HighDM_Normalization_met_2016metWithLLmetWithLLZToLLstack
       # DataMC_Electron_HighDM_Normalization_met_2016metWithLLmetWithLLNoZToLLstack
       # DataMC_Muon_LowDM_Normalization_met_2016metWithLLmetWithLLZToLLstack
       # DataMC_Muon_LowDM_Normalization_met_2016metWithLLmetWithLLNoZToLLstack
       # DataMC_Muon_HighDM_Normalization_met_2016metWithLLmetWithLLZToLLstack
       # DataMC_Muon_HighDM_Normalization_met_2016metWithLLmetWithLLNoZToLLstack

        self.histos[year] = {}
        for particle in self.particles:
            self.histos[year][particle] = {}
            for region in self.regions:
                if self.useAllMC:
                    # using ZToLL and NoZToLL MC for normalization 
                    self.histos[year][particle][region] = { 
                        "Data"     : "DataMC_" + particle + "_" + region + "_Normalization_met_" + year + "metWithLLmetWithLLDatadata",
                        "ZToLL"    : "DataMC_" + particle + "_" + region + "_Normalization_met_" + year + "metWithLLmetWithLLZToLLstack",
                        "NoZToLL"  : "DataMC_" + particle + "_" + region + "_Normalization_met_" + year + "metWithLLmetWithLLNoZToLLstack",
                    }
                else:
                    # using only DY and ttbar for normalization 
                    self.histos[year][particle][region] = { 
                        "Data"     : "DataMC_" + particle + "_" + region + "_met_" + year + "metWithLLmetWithLLDatadata",
                        "ZToLL"    : "DataMC_" + particle + "_" + region + "_met_" + year + "metWithLLmetWithLLDYstack",
                        "NoZToLL"  : "DataMC_" + particle + "_" + region + "_met_" + year + "metWithLLmetWithLLt#bar{t}stack",
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

		print self.variable + "/" + self.histos[year][particle][region]["Data"]

                h_ZToLL   = f.Get(self.variable + "/" + self.histos[year][particle][region]["ZToLL"])
                h_NoZToLL = f.Get(self.variable + "/" + self.histos[year][particle][region]["NoZToLL"])

    		#
		#
		#
		#
                #############################
                # Calculating Normalization #
                #############################
		#
		#
		#
		#
		
		bin_1 = 0
		bin_2 = h_ZToLL.FindBin(1000)
		I_Z = h_ZToLL.Integral(bin_1, bin_2)
		I_B = h_NoZToLL.Integral(bin_1, bin_2)
		I_D = h_Data.Integral(bin_1, bin_2)
		
		print I_Z, "  ", I_B, "  ", I_D
		
                norm  = (I_D-I_B)/I_Z
		               
		self.norm_map_tex[era][particle][region]["R_Z"] = norm
		

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
            self.writeLine("\\begin{tabular}{|c|c|c|c|}")
            self.writeLine("\hline Selection & Era & CR & $R_Z$ \\\\")
            # write values to table
            for region in self.regions:
                for era in self.eras:
                    for particle in self.particles:
                        R_Z = self.norm_map_tex[era][particle][region]["R_Z"]
                        region_tex = self.regions_tex[region]
                        era_tex = era.replace("_", " ")
                        self.writeLine("\hline {0} & {1} & {2} & {3} \\\\".format(region_tex, era_tex, particle, R_Z))
            self.writeLine("\hline")
            self.writeLine("\end{tabular}")
            # end table
            self.writeLine("\end{table}")
            # end document
            self.writeLine("\end{document}")

def main():
    verbose = True
    useAllMC = True
    N = Normalization(useAllMC, verbose)
    N.getNormAndError("/uscms/home/arosado/nobackup/YOURWORKINGAREA/CMSSW_10_2_9/src/ZInvisible/Tools/condor/2018/result.root", "2018")
    N.makeTexFile("normalization.tex")


if __name__ == "__main__":
    main()


