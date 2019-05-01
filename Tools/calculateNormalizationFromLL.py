# calculateNormalizationFromLL.py
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True


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
    a11 = 2.0
    a12 = 1.0
    a21 = 1.0
    a22 = 2.0
    # Data
    x1  = 1.1
    x2  = 1.0
    # normalization
    norm = calcNorm(a11, a12, a21, a22, x1, x2, False)
    b1 = norm[0]
    b2 = norm[1]
    print ("R_Z = {0}".format(b1))
    print ("R_T = {0}".format(b2))

if __name__ == "__main__":
    main()






