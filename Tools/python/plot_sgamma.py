# plot_sgamma.py

import ROOT
import json
import numpy as np
from tools import plot

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def run(era):
    verbose = 2
    print "---------- Running {0} ----------".format(era)
    # for datacard, the allowed Sgamma range is [0.01, 5]
    minSgamma = 0.01
    maxSgamma = 5.0
    
    fileSB = "results/SearchBinResults.json"
    fileYields = "datacard_inputs/zinv_yields_" + era + ".json"
    fileBinMap = "dc_BkgPred_BinMaps_master.json"
    with open(fileSB, "r") as f:
        sbResults = json.load(f)
    with open(fileYields, "r") as f:
        yieldResults = json.load(f)
    with open(fileBinMap, "r") as f:
        binMap = json.load(f)

    # histograms
    nbins_1 = 50
    nbins_2 = 60
    limits_1 = [-5, 30]
    limits_2 = [-1, 5]
    # multiple sets of histograms with different limits
    h_sgamma_searchbins_1       = ROOT.TH1F("h_sgamma_searchbins_1",        "h_sgamma_searchbins_1",        nbins_1, limits_1[0], limits_1[1])
    h_sgamma_crunits_1          = ROOT.TH1F("h_sgamma_crunits_1",           "h_sgamma_crunits_1",           nbins_1, limits_1[0], limits_1[1])
    h_sgamma_searchbins_2       = ROOT.TH1F("h_sgamma_searchbins_2",        "h_sgamma_searchbins_2",        nbins_2, limits_2[0], limits_2[1])
    h_sgamma_crunits_2          = ROOT.TH1F("h_sgamma_crunits_2",           "h_sgamma_crunits_2",           nbins_2, limits_2[0], limits_2[1])

    # --------------------------------------------- #
    # get bin and sgamma values for each search bin #
    # --------------------------------------------- #
    
    sgammaForSearchBins = []
    searchBinAndSgamma = list((int(b), sbResults[era][b]["shape"]) for b in sbResults[era])
    searchBinAndSgamma.sort(key = lambda x: x[0])
    for x in searchBinAndSgamma:
        b       = x[0]
        sgamma  = x[1]
        sgammaForSearchBins.append(sgamma)
        h_sgamma_searchbins_1.Fill(sgamma)
        h_sgamma_searchbins_2.Fill(sgamma)
        if verbose > 1 and (sgamma < minSgamma or sgamma > maxSgamma):
            print "search bin {0}, sgamma = {1}".format(b, sgamma)
   
    # ----------------------------------------------------- #
    # get bin and sgamma values for each control region bin #
    # ----------------------------------------------------- #
    
    total_phocr_data    = 0
    total_phocr_gjets   = 0
    total_phocr_back    = 0
    sgammaForCRUnits = []
    CRBinNames = {}
    for binName in binMap["unitCRNum"]["phocr"]:
        b =  int(binMap["unitCRNum"]["phocr"][binName])
        CRBinNames[b] = binName
    
    # keys for yieldResults["yieldsMap"]:
    # znunu
    # phocr_back
    # phocr_data
    # phocr_gjets
    
    for b in xrange(len(CRBinNames)):
        binName = CRBinNames[b]
        
        phocr_data  = yieldResults["yieldsMap"]["phocr_data"][binName][0]
        phocr_gjets = yieldResults["yieldsMap"]["phocr_gjets"][binName][0]
        phocr_back  = yieldResults["yieldsMap"]["phocr_back"][binName][0]
        
        # add to totals
        total_phocr_data    += phocr_data
        total_phocr_gjets   += phocr_gjets
        total_phocr_back    += phocr_back
        
        sgamma = -999
        den = phocr_gjets + phocr_back
        if den > 0.0:
            sgamma = phocr_data / den
        else:
            print "WARNING: CR bin {0}, denominator = {1}".format(b, den)

        sgammaForCRUnits.append(sgamma)
        h_sgamma_crunits_1.Fill(sgamma)
        h_sgamma_crunits_2.Fill(sgamma)
        
        if verbose > 1 and (sgamma < minSgamma or sgamma > maxSgamma):
            print "CR bin {0}, sgamma = {1}; phocr_data = {2}, phocr_gjets = {3}, phocr_back = {4}".format(b, sgamma, phocr_data, phocr_gjets, phocr_back)
        
        if verbose > 2 and phocr_data <= 2.0:
            print "CR bin {0}: phocr_data = {1}".format(b, phocr_data)
        
        if verbose > 2 and phocr_gjets <= 2.0:
            print "CR bin {0}: phocr_gjets = {1}".format(b, phocr_gjets)

    # total normalization
    total_norm = total_phocr_data / (total_phocr_gjets + total_phocr_back)

    histograms_1 = [h_sgamma_searchbins_1, h_sgamma_crunits_1]
    histograms_2 = [h_sgamma_searchbins_2, h_sgamma_crunits_2]
    labels = ["sgamma_searchbins", "sgamma_crunits"]
    
    # --- plot --- #
    # plot(histograms, labels, name, title, x_title, y_title, x_min, x_max, y_min, y_max, era, plot_dir, showStats=False, normalize=False, setLog=False)
    plot(histograms_1, labels, "sgamma_binning1", "Sgamma for " + era, "Sgamma", "Events", limits_1[0], limits_1[1], 0.0, 200.0, era, "more_plots", showStats=True) 
    plot(histograms_2, labels, "sgamma_binning2", "Sgamma for " + era, "Sgamma", "Events", limits_2[0], limits_2[1], 0.0, 60.0,  era, "more_plots", showStats=True) 
   
    if verbose > 0:
        print "Total data  = {0}".format(total_phocr_data)
        print "Total gjets = {0}".format(total_phocr_gjets)
        print "Total other background = {0}".format(total_phocr_back)
        print "Total normalization: data / (gjets + other back) = {0}".format(total_norm)
        print "Sgamma in search bins: mean = {0:.2f}, std_dev = {1:.2f}, min = {2:.2f}, max = {3:.2f}".format(np.mean(sgammaForSearchBins), np.std(sgammaForSearchBins), np.amin(sgammaForSearchBins), np.amax(sgammaForSearchBins))
        print "Sgamma in control bins: mean = {0:.2f}, std_dev = {1:.2f}, min = {2:.2f}, max = {3:.2f}".format(np.mean(sgammaForCRUnits), np.std(sgammaForCRUnits), np.amin(sgammaForCRUnits), np.amax(sgammaForCRUnits))

    del h_sgamma_searchbins_1
    del h_sgamma_searchbins_2
    del h_sgamma_crunits_1
    del h_sgamma_crunits_2

def main():
    eras = ["2016", "2017", "2018", "Run2"]
    for era in eras:
        run(era)

if __name__ == "__main__":
    main()

