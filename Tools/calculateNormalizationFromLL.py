# calculateNormalizationFromLL.py
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Normalization:
    def __init__(self):
        self.norm_map = {}
        self.output_file = 0 
        self.years = []
        self.particles = ["Electron", "Muon"]
        self.regions   = ["LowDM", "HighDM"]
        self.regions_tex = {
                             "LowDM"  : "Low $\Delta m$",
                             "HighDM" : "High $\Delta m$"
                           }

    def calcNorm(self, A, x, verbose=False):
        # --- Solve for b --- #
        # x = A * b           #
        # b = A^(-1) * x      #
        # ------------------- #
        
        # inverse of matrix
        try:
            Ainverse = np.linalg.inv(A)
        except:
            print "ERROR: Cannot calculate inverse because matrix is singular (or possibly other error...)."
            return [-999, -999]
        b = np.dot(Ainverse, x)
        if verbose:
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
    
    def calcError(self, A, A_error, x, x_error, verbose):
        # inverse of matrix
        try:
            Ainverse = np.linalg.inv(A)
        except:
            print "ERROR: Cannot calculate inverse because matrix is singular (or possibly other error...)."
            return [-999, -999]
        Ainverse_error = self.getMatrixInverseError(A, A_error)
        # from equation b = A^(-1) * x
        b_error = np.zeros(2)
        b1_error_1 = self.getMultiplicationError(Ainverse[0][0] * x[0], Ainverse[0][0], Ainverse_error[0][0], x[0], x_error[0])
        b1_error_2 = self.getMultiplicationError(Ainverse[0][1] * x[1], Ainverse[0][1], Ainverse_error[0][1], x[1], x_error[1])
        b2_error_1 = self.getMultiplicationError(Ainverse[1][0] * x[0], Ainverse[1][0], Ainverse_error[1][0], x[0], x_error[0])
        b2_error_2 = self.getMultiplicationError(Ainverse[1][1] * x[1], Ainverse[1][1], Ainverse_error[1][1], x[1], x_error[1])
        b_error[0] = self.getAdditionError(b1_error_1, b1_error_2)
        b_error[1] = self.getAdditionError(b2_error_1, b2_error_2)
        if verbose:
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
    
    def getMultiplicationError(self, q, x, dx, y, dy):
        # q = x * y 
        # q = x / y 
        # dq = q * sqrt( (dx/x)^2 + (dy/y)^2) )
        return abs(q * np.sqrt( (dx/x)**2 + (dy/y)**2 ))
    
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
    
    def getNormAndError(self, file_name, year, verbose):
        self.years.append(year)
        self.norm_map[year] = {}
        print "------------------------------------------------------------------------------"
        print "Year: {0}".format(year)
        print "File: {0}".format(file_name)
        f = ROOT.TFile(file_name)
        # variable is also TDirectoryFile that holds histograms 
        variable = "bestRecoZM"
        # histograms
        # example: DataMC_Electron_LowDM_bestRecoZM_0to400_2016bestRecoZMbestRecoZMDatadata
        histos = {}
        for particle in self.particles:
            histos[particle] = {}
            for region in self.regions:
                histos[particle][region] = { 
                    "Data"     : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMDatadata",
                    "DY"       : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMDYstack",
                    "TTbar"    : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMt#bar{t}stack",
                    "Single t" : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMSingle tstack",
                    "Rare"     : "DataMC_" + particle + "_" + region + "_bestRecoZM_0to400_" + year + "bestRecoZMbestRecoZMRarestack"
                }
        
        for particle in self.particles:
            self.norm_map[year][particle] = {}
            print particle
            for region in self.regions:
                self.norm_map[year][particle][region] = {}
                print region
                h_Data  = f.Get(variable + "/" + histos[particle][region]["Data"])
                h_DY    = f.Get(variable + "/" + histos[particle][region]["DY"])
                h_TTbar = f.Get(variable + "/" + histos[particle][region]["TTbar"])
    
                #############################
                # Calculating Normalization #
                #############################
                minMass = 50.0
                Zmass = 91.0
                minZmass = Zmass - 10.0
                maxZmass = Zmass + 10.0
                bin_1 = h_Data.FindBin(minMass)
                bin_2 = h_Data.FindBin(minZmass)
                bin_3 = h_Data.FindBin(maxZmass)
                bin_4 = h_Data.GetNbinsX() + 1
                if verbose:
                    for i, bin_i in enumerate([bin_1, bin_2, bin_3, bin_4]):
                        print "bin_{0}: {1}".format(i, bin_i)
                
                # pass error to IntegralAndError as Double_t & error 
                # use ROOT.Double for pass-by-ref of doubles
                # Double_t TH1::IntegralAndError  (Int_t binx1,
                #                                  Int_t binx2,
                #                                  Double_t & error,
                #                                  Option_t * option = "" 
                #                                 ) const
                
                # MC matrix (A)
                # a11: DY on Z mass peak
                # a12: TTbar on Z mass peak
                # a21: DY off Z mass peak
                # a22: TTbar off Z mass peak
                # Data (x)
                # x1: data on Z mass peak
                # x2: data off Z mass peak
                # Normalization (b)
                # b1: R_Z (DY normalization)
                # b2: R_T (TTbar normalization)
    
                # MC matrix
                a11_error = ROOT.Double()
                a12_error = ROOT.Double()
                a21_error_1 = ROOT.Double()
                a21_error_2 = ROOT.Double()
                a22_error_1 = ROOT.Double()
                a22_error_2 = ROOT.Double()
                a11 = h_DY.IntegralAndError(bin_2, bin_3, a11_error) 
                a12 = h_TTbar.IntegralAndError(bin_2, bin_3, a12_error) 
                a21 = h_DY.IntegralAndError(bin_1, bin_2, a21_error_1) + h_DY.IntegralAndError(bin_3, bin_4, a21_error_2)
                a22 = h_TTbar.IntegralAndError(bin_1, bin_2, a22_error_1) + h_TTbar.IntegralAndError(bin_3, bin_4, a22_error_2)
                A = np.array([[a11, a12], [a21, a22]])
                a21_error = self.getAdditionError(a21_error_1, a21_error_2) 
                a22_error = self.getAdditionError(a22_error_1, a22_error_2) 
                A_error = np.array([[a11_error, a12_error], [a21_error, a22_error]])
                # Data
                x1_error = ROOT.Double()
                x2_error_1 = ROOT.Double()
                x2_error_2 = ROOT.Double()
                x1 = h_Data.IntegralAndError(bin_2, bin_3, x1_error)
                x2 = h_Data.IntegralAndError(bin_1, bin_2, x2_error_1) + h_Data.IntegralAndError(bin_3, bin_4, x2_error_2)
                x = [x1, x2]
                x2_error = self.getAdditionError(x2_error_1, x2_error_2)
                x_error = [x1_error, x2_error]
                # calculate normalization and error
                norm  = self.calcNorm(A, x, verbose)
                error = self.calcError(A, A_error, x, x_error, verbose)
                b1 = norm[0]
                b2 = norm[1]
                b1_error = error[0]
                b2_error = error[1]
                
                R_Z_print = "R_Z = {0:.4f} +/- {1:.4f}".format(b1, b1_error)
                R_T_print = "R_T = {0:.4f} +/- {1:.4f}".format(b2, b2_error)
                R_Z_tex = "${0:.4f} \pm {1:.4f}$".format(b1, b1_error)
                R_T_tex = "${0:.4f} \pm {1:.4f}$".format(b2, b2_error)
                print R_Z_print 
                print R_T_print
                self.norm_map[year][particle][region]["R_Z"] = R_Z_tex 
                self.norm_map[year][particle][region]["R_T"] = R_T_tex 
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
            self.writeLine("\hline Year & CR & Selection & $R_Z$ & $R_T$ \\\\")
            # write values to table
            for year in self.years:
                for particle in self.particles:
                    for region in self.regions:
                        R_Z = self.norm_map[year][particle][region]["R_Z"] 
                        R_T = self.norm_map[year][particle][region]["R_T"] 
                        region_tex = self.regions_tex[region]
                        self.writeLine("\hline {0} & {1} & {2} & {3} & {4} \\\\".format(year, particle, region_tex, R_Z, R_T))
            self.writeLine("\hline")
            self.writeLine("\end{tabular}")
            # end table
            self.writeLine("\end{table}")
            # end document
            self.writeLine("\end{document}")

def main():
    verbose = False
    N = Normalization()
    # May 5, 2019 Results
    #N.getNormAndError("condor/DataMC_2016_submission_2019-05-05_21-57-41/result.root", "2016", verbose)
    #N.getNormAndError("condor/DataMC_2017_submission_2019-05-05_22-28-09/result.root", "2017", verbose)
    #N.getNormAndError("condor/DataMC_2018_submission_2019-05-05_22-44-26/result.root", "2018", verbose)
    # May 9, 2019 Results
    N.getNormAndError("condor/DataMC_2016_submission_2019-05-09_17-19-42/result.root", "2016", verbose)
    N.getNormAndError("condor/DataMC_2017_submission_2019-05-09_17-16-54/result.root", "2017", verbose)
    N.getNormAndError("condor/DataMC_2018_submission_2019-05-09_17-15-04/result.root", "2018", verbose)
    N.makeTexFile("normalization.tex")


if __name__ == "__main__":
    main()


