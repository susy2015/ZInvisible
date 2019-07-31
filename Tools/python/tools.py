# tools.py

from colors import getColorIndex
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

ERROR_CODE = -999

def setupHist(hist, title, x_title, y_title, color, y_min, y_max):
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    y_axis.SetRangeUser(y_min, y_max)
    hist.SetTitle(title)
    x_axis.SetTitle(x_title)
    y_axis.SetTitle(y_title)
    hist.SetStats(ROOT.kFALSE)
    hist.SetLineColor(getColorIndex(color))
    hist.SetLineWidth(3)

def getAdditionError(dx, dy):
    # q  = x + y
    # q  = x - y
    # dq = sqrt( dx^2 + dy^2 )
    return abs(np.sqrt( dx**2 + dy**2 ))

def getConstantMultiplicationError(a, dx):
    # a is a constant
    # q  = a * x
    # dq = |a| * dx
    return abs(a * dx)

def getMultiplicationError(q, x, dx, y, dy):
    # q = x * y 
    # q = x / y 
    # dq = q * sqrt( (dx/x)^2 + (dy/y)^2) )
    if x == 0.0 or y == 0.0:
        print "ERROR in getMultiplicationError(): Cannot divide by zero."
        return ERROR_CODE
    return abs(q * np.sqrt( (dx/x)**2 + (dy/y)**2 ))

def getMultiplicationErrorList(q, x_list, dx_list):
    # q = x * y * ... 
    # q = x / y / ...
    # dq = q * sqrt( (dx/x)^2 + (dy/y)^2) + ... )
    if len(x_list) != len(dx_list):
        print "ERROR in getMultiplicationErrorList(): x_list and dx_list do not have the same length."
        return ERROR_CODE
    s = 0.0
    for i in xrange(len(x_list)):
        if x_list[i] == 0.0:
            print "ERROR in getMultiplicationErrorList(): Cannot divide by zero."
            return ERROR_CODE
        s += (dx_list[i] / x_list[i]) ** 2
    return abs(q * np.sqrt(s))

def getMatrixDeterminantError(A, A_error):
    # get error for determinant
    # |A| = a11 * a22 - a12 * a21
    error_1 = getMultiplicationError(A[0][0] * A[1][1], A[0][0], A_error[0][0], A[1][1], A_error[1][1])
    error_2 = getMultiplicationError(A[0][1] * A[1][0], A[0][1], A_error[0][1], A[1][0], A_error[1][0])
    det_error = getAdditionError(error_1, error_2)
    return det_error

def getMatrixInverseError(A, A_error):
    # from equation for A^(-1)
    det = np.linalg.det(A)
    det_error = getMatrixDeterminantError(A, A_error)
    Ainverse_error = np.zeros((2,2))
    Ainverse_error[0][0] = getMultiplicationError( A[1][1] / det,  A[1][1], A_error[1][1], det, det_error)
    Ainverse_error[0][1] = getMultiplicationError(-A[0][1] / det, -A[0][1], A_error[0][1], det, det_error)
    Ainverse_error[1][0] = getMultiplicationError(-A[1][0] / det, -A[1][0], A_error[1][0], det, det_error)
    Ainverse_error[1][1] = getMultiplicationError( A[0][0] / det,  A[0][0], A_error[0][0], det, det_error)
    return Ainverse_error


