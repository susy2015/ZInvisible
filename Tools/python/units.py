# units.py
# here we go.....
# use those CR and SR units!

import json

def show(jsonFile, key, title):
    with open(jsonFile, "r") as inputFile:
        unitMap = json.load(inputFile)
        # invert map
        binMap = {v: k for k, v in unitMap[key].iteritems()}
        nBins = len(binMap)
        bins = list(str(x) for x in xrange(nBins))
        print "{0} : {1} bins".format(title, nBins)
        for b in bins: 
            print "{0} : {1}".format(b, binMap[b])

def run():
    jsonFile = "dc_BkgPred_BinMaps_master.json"
    show(jsonFile, "binNum",    "Search Bins")
    show(jsonFile, "unitCRNum", "Control Region Units")
    show(jsonFile, "unitSRNum", "Search Region Units")

if __name__ == "__main__":
    run()

