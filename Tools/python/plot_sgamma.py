# plot_sgamma.py

import ROOT
import json


# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def run(era):
    print "Running {0}".format(era)
    minSgamma = 0.2
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

    sgammaSearchBins = list((int(b), sbResults[era][b]["shape"]) for b in sbResults[era])
    sgammaSearchBins.sort(key = lambda x: x[0])
    for x in sgammaSearchBins:
        b       = x[0]
        sgamma  = x[1]
        if sgamma < minSgamma or sgamma > maxSgamma:  
            print "bin {0}, sgamma = {1}".format(b, sgamma)

def main():
    eras = ["2016", "2017", "2018", "Run2"]
    for era in eras:
        run(era)

if __name__ == "__main__":
    main()

