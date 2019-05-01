# calculateNormalizationFromLL.py
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def calcNorm(a11, a12, a21, a22, x1, x2, verbose=False):
    # --- Solve for b --- #
    # x = A * b           #
    # b = A^(-1) * x      #
    # ------------------- #
    A = np.array([[a11, a12], [a21, a22]])
    x = [x1, x2]
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
        print "MC matrix:"
        print A
        print "Inverse of MC matrix:"
        print Ainverse
        print "b:"
        print b
        print "x:"
        print x
        print "Calculated x (for validation):"
        print x_calc
        print "Difference between x and calculated x:"
        print x - x_calc
        print "-----------------------------------"
        


    return b


def main():
    verbose = False
    f = ROOT.TFile("condor/histos_ElectronControlRegionSelection_2016_01_May_2019_2/result.root")
    # histograms
    # Low DM
    # variable is also TDirectoryFile that holds histograms 
    variable = "bestRecoZM"
    histos = { 
        "LowDM": {
            "Data"     : "DataMC_Electron_LowDM_bestRecoZM_2016bestRecoZMbestRecoZMDatadata",
            "DY"       : "DataMC_Electron_LowDM_bestRecoZM_2016bestRecoZMbestRecoZMDYstack",
            "TTbar"    : "DataMC_Electron_LowDM_bestRecoZM_2016bestRecoZMbestRecoZMt#bar{t}stack",
            "Single t" : "DataMC_Electron_LowDM_bestRecoZM_2016bestRecoZMbestRecoZMSingle tstack",
            "Rare"     : "DataMC_Electron_LowDM_bestRecoZM_2016bestRecoZMbestRecoZMRarestack"
        },
        "HighDM" : {
            "Data"     : "DataMC_Electron_HighDM_bestRecoZM_2016bestRecoZMbestRecoZMDatadata",
            "DY"       : "DataMC_Electron_HighDM_bestRecoZM_2016bestRecoZMbestRecoZMDYstack",
            "TTbar"    : "DataMC_Electron_HighDM_bestRecoZM_2016bestRecoZMbestRecoZMt#bar{t}stack",
            "Single t" : "DataMC_Electron_HighDM_bestRecoZM_2016bestRecoZMbestRecoZMSingle tstack",
            "Rare"     : "DataMC_Electron_HighDM_bestRecoZM_2016bestRecoZMbestRecoZMRarestack"
        }
    } 
    for key in histos:
        print key
        h_Data  = f.Get(variable + "/" + histos[key]["Data"])
        h_DY    = f.Get(variable + "/" + histos[key]["DY"])
        h_TTbar = f.Get(variable + "/" + histos[key]["TTbar"])

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
        a11 = h_DY.Integral(bin_2, bin_3) 
        a12 = h_TTbar.Integral(bin_2, bin_3) 
        a21 = h_DY.Integral(bin_1, bin_2) + h_DY.Integral(bin_3, bin_4)
        a22 = h_TTbar.Integral(bin_1, bin_2) + h_TTbar.Integral(bin_3, bin_4)
        # Data
        x1  = h_Data.Integral(bin_2, bin_3)
        x2  = h_Data.Integral(bin_1, bin_2) + h_Data.Integral(bin_3, bin_4)
        # normalization
        norm = calcNorm(a11, a12, a21, a22, x1, x2, verbose)
        b1 = norm[0]
        b2 = norm[1]
        print ("R_Z = {0}".format(b1))
        print ("R_T = {0}".format(b2))

if __name__ == "__main__":
    main()






