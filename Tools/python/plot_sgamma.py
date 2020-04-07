# plot_sgamma.py

import ROOT
import json


# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def run(era):
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

    # get bin and sgamma values for each search bin
    sgammaSearchBins = list((int(b), sbResults[era][b]["shape"]) for b in sbResults[era])
    sgammaSearchBins.sort(key = lambda x: x[0])
    for x in sgammaSearchBins:
        b       = x[0]
        sgamma  = x[1]
        if sgamma < minSgamma or sgamma > maxSgamma:  
            print "search bin {0}, sgamma = {1}".format(b, sgamma)
    # get bin and sgamma values for each control region bin
    sgammaCRUnits = []
    CRBinNames = {}
    for binName in binMap["unitCRNum"]["phocr"]:
        b =  int(binMap["unitCRNum"]["phocr"][binName])
        CRBinNames[b] = binName
    
    #for b in xrange(len(CRBinNames)):
    #    print "{0} : {1}".format(b, CRBinNames[b])
    #for key in yieldResults["yieldsMap"]:
    #    print key
    
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
        sgamma = -999
        den = phocr_gjets + phocr_back
        if den > 0.0:
            sgamma = phocr_data / den
        else:
            print "WARNING: CR bin {0}, denominator = {1}".format(b, den)
        
        if sgamma < minSgamma or sgamma > maxSgamma:  
            print "CR bin {0}, sgamma = {1}".format(b, sgamma)
        
        #if phocr_data < 2.0:
        #    print "CR bin {0}: phocr_data = {1}".format(b, phocr_data)







def main():
    eras = ["2016", "2017", "2018", "Run2"]
    for era in eras:
        run(era)

if __name__ == "__main__":
    main()



