# tools.py

from colors import getColorIndex
import numpy as np
import ROOT
import re
# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS Poisson Errors
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StatisticsCommittee
# https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars

# poisson error for 0 
#ERROR_ZERO = 1.841022
# small error for 0
ERROR_ZERO = 0.001
# general error code
ERROR_CODE = -999
# infinity
INF = float("inf")

# Important: be careful with order of replacing characters if keys match values
tex_map = {
            "LowDM"     : "\\mathrm{Low}~\\Delta m",
            "HighDM"    : "\\mathrm{High}~\\Delta m",
            "met"       : "\\cancel{E}_{T}",
            "ht"        : "H_{T}",
            "NB"        : "N_{b}",
            "NSV"       : "N_{sv}",
            "NJ"        : "N_{j}"
          }

# invert map
def invert(map_):
    inv_map = {v: k for k, v in map_.iteritems()} 
    return inv_map

# when loading json files, strings are loaded as unicode
# unicode can cause problems when using ROOT or Latex
# return map with string values
def stringifyMap(rawMap):
    returnMap = rawMap
    for key1 in returnMap:
        for key2 in returnMap[key1]:
            returnMap[key1][key2] = str(returnMap[key1][key2])
    return returnMap

# https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
# https://docs.python.org/3/whatsnew/3.5.html#pep-485-a-function-for-testing-approximate-equality
# isclose(): used to compare floats up to some relative tolerance and absolute tolerance
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

# get bin error
# return given error for bins with value 0.0
def getBinError(value, value_error, error_for_zero):
    if value == 0 and value_error == 0:
        return error_for_zero
    else:
        return value_error

# used for sorting by first element
def takeFirst(elem):
    return elem[0]

# function to get met binning from map based on selection
# one selection can have mutliple MET binnings
# return list of numpy arrays 
# return dictionary that has
# - list of list of bin names
# - list of numpy arrarys of bin edges
def getMETBinEdges(binMap, selection):
    verbose = False
    max_met = 1000.0
    previous_met_1 = -1.0 * INF
    previous_met_2 = -1.0 * INF
    temp_names = []
    names = []
    met_bins = []
    temp_array = []
    result = []
    for b in binMap:
        cut = removeCuts(binMap[b]["selection"], ["NSV", "MET"])
        if verbose:
            print "b: {0}, selection: {1}, cut: {2}".format(b, selection, cut)
        # get array of bins matching selection used for shape factor
        if selection == cut:
            if temp_array:
                # if temp_array as elements, don't append if previous element is the same
                if binMap[b]["met"] != temp_array[-1][1]:
                    temp_array.append((int(b), binMap[b]["met"]))
            else:
                # append to temp_array if it is empty
                temp_array.append((int(b), binMap[b]["met"]))
    # sort by bin number; assume MET bins increase with bin number (unless new binning starts)
    temp_array.sort(key = takeFirst) 
    for i, elem in enumerate(temp_array):
        # use regex to get values
        # example: met_450to550
        # use first value for MET array (e.g. 450)
        # use second value to check for INF (e.g. 550)
        name = elem[1]
        m = re.search("met_(.*)to(.*)", name)
        value_1 = float(m.group(1))
        value_2 = float(m.group(2))
        # keep adding MET bins as they increase
        if value_1 > previous_met_1:
            temp_names.append(name)
            met_bins.append(value_1)
        # new MET binning starts
        # be careful to append if we have moved to a new binning or if this is the last bin
        # also append previous value 2 if it is not INF
        isNewBinning = bool(value_1 <= previous_met_1)
        isLastBin    = bool(i == len(temp_array) - 1)
        
        if isNewBinning or isLastBin:
            if isNewBinning and (previous_met_2 != INF):
                met_bins.append(previous_met_2)
            if isLastBin and (value_2 != INF):
                met_bins.append(value_2)
            met_bins.append(max_met)
            met_array = np.array(met_bins)
            # avoid duplicate entries: only append if this binning is not already in names 
            if temp_names not in names: 
                names.append(temp_names)
                result.append(met_array)
            # reset lists
            if isNewBinning:
                temp_names = [name]
                met_bins = [value_1]

        previous_met_1 = value_1
        previous_met_2 = value_2
    
    d = {"names": names, "xbins": result}
    if verbose:
        print "{0}: {1}".format(names, result)
    return d

# get selections from json file
def getSelections(bin_map, bin_type, remove_cut_list=[], verbose=False): 
    lowdm_cuts  = []
    highdm_cuts = []
    for b in bin_map:
        r = bin_map[b]["region"]
        cut = bin_map[b]["selection"]
        # only remove cut if remove_cut_list is given
        if remove_cut_list:
            cut = removeCuts(cut, remove_cut_list)
        if r == "LowDM":
            lowdm_cuts.append(cut)
        elif r == "HighDM":
            highdm_cuts.append(cut)

    # remove duplicates from list
    lowdm_cuts  = list(dict.fromkeys(lowdm_cuts))
    highdm_cuts = list(dict.fromkeys(highdm_cuts))
    if verbose:
        print "{0} Low DM cuts: {1}".format(bin_type,  lowdm_cuts)
        print "{0} High DM cuts: {1}".format(bin_type, highdm_cuts)
    result = {}
    result["LowDM"]  = lowdm_cuts
    result["HighDM"] = highdm_cuts
    return result


# special selection with multiple cuts (say met, ht...)
# e.g. met_250to400, met_400toINF
def getTexMultiCut(selection):
    # be careful with order of search and replace
    result = selection.split("_") 
    var  = result[0]
    cuts = result[1]
    cuts = cuts.split("to")
    for key in tex_map:
        var = var.replace(key, tex_map[key])
    if cuts[1] == "INF":
        expression = "${0} \geq {1}$".format(var, cuts[0])
    else:
        expression = "${0} \leq {1} < {2}$".format(cuts[0], var, cuts[1])
    return expression

# get latex expression of selection
def getTexSelection(selection):
    # be careful with order of search and replace
    result = selection.split("_") 
    for i in xrange(len(result)):
        r = result[i]
        for key in tex_map:
            r = r.replace(key, tex_map[key])
        # some keywords contain these strings... be careful
        if "delta" not in r and "Delta" not in r:
            if "eq" in r:
                r = r.replace("eq", " = ")
            elif "ge" in r:
                r = r.replace("ge", " \geq ")
            elif "le" in r:
                r = r.replace("le", " \leq ")
            elif "gt" in r:
                r = r.replace("gt", " > ")
            elif "lt" in r:
                r = r.replace("lt", " < ")
        result[i] = r
    
    expression = "$" + ", ".join(result) + "$"
    return expression

# remove cuts if pattern in cuts
def removeCuts(cutString, pattern_list, delim = "_"):
    newString = cutString
    for pattern in pattern_list:
        cutList = newString.split(delim)
        #print "cutString: {0}, newString: {1}, cutList: {2}".format(cutString, newString, cutList)
        newList = list(c for c in cutList if pattern not in c)
        newString = delim.join(newList)
    return newString

# get ratio
def getNormalizedRatio(h_num, h_den):
    # number of events for normalization
    nNum  = h_num.Integral(0, h_num.GetNbinsX() + 1)
    nDen  = h_den.Integral(0, h_den.GetNbinsX() + 1)
    ratio = float(nNum) / float(nDen)
    
    h_den_normalized = h_den.Clone("h_den_normalized")
    h_den_normalized.Scale(ratio)
    
    # normalized ratio
    h_ratio_normalized = h_num.Clone("h_ratio_normalized")
    h_ratio_normalized.Divide(h_den_normalized)
    return h_ratio_normalized

# setup histogram
def setupHist(hist, title, x_title, y_title, color, y_min, y_max, adjust=False):
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    y_axis.SetRangeUser(y_min, y_max)
    hist.SetTitle(title)
    # adjust font sizes only for some plots
    if adjust:
        x_axis.SetTitleSize(0.06)
        y_axis.SetTitleSize(0.06)
        x_axis.SetTitleOffset(0.75)
        y_axis.SetTitleOffset(0.75)
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

def getAdditionErrorList(dx_list):
    # q  = x + y + ...
    # q  = x - y - ...
    # dq = sqrt( dx^2 + dy^2 + ... )
    dx2_list = [dx**2 for dx in dx_list]
    return abs(np.sqrt(sum(dx2_list)))

def getConstantMultiplicationError(a, dx):
    # a is a constant
    # q  = a * x
    # dq = |a| * dx
    return abs(a) * dx

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


