# units.py
# here we go.....
# use those CR and SR units!

import json
from search_bins import SearchBins, ValidationBins, CRUnitBins, SRUnitBins
from tools import invert


def show(jsonFile, key, title):
    with open(jsonFile, "r") as inputFile:
        unitMap = json.load(inputFile)
        # invert map
        binMap = invert(unitMap[key])
        nBins = len(binMap)
        bins = list(str(x) for x in xrange(nBins))
        print "{0} ({1}) : {2} bins".format(title, key, nBins)
        for b in bins: 
            print "{0} : {1}".format(b, binMap[b])

# save final results in json file
def saveResults(inFile, outFile, CRunits, SRunits, eras):
    map_out                                 = {}
    map_out["yieldsMap"]                    = {}
    map_out["yieldsMap"]["znunu"]           = {}
    map_out["yieldsMap"]["phocr_data"]      = {}
    map_out["yieldsMap"]["phocr_gjets"]     = {}
    map_out["yieldsMap"]["phocr_back"]      = {}
    with open (inFile, "r") as f_in:
        map_in = json.load(f_in)
        # invert maps
        SearchBin_map   = invert(map_in["binNum"])
        CRunit_map      = invert(map_in["unitCRNum"])
        SRunit_map      = invert(map_in["unitSRNum"])

    # final yields for json file
    # znunu:        Z nunu MC in SR unit bins
    # phocr_data:   Photon Data in photon CR 
    # phocr_gjets:  GJets MC in photon CR
    # phocr_back:   Other MC (QCD, TTGJets, etc) in photon CR
    
    # combine results over all eras
    # try only 2016 now for development
    # try Run 2 now
    era = "Run2"
    
    print "Search Region Units"
    for b in SRunits.binValues[era]:
        name        = SRunit_map[b]
        mc          = SRunits.binValues[era][b]["mc"]
        mc_error    = SRunits.binValues[era][b]["mc_error"]

        # save value and error
        map_out["yieldsMap"]["znunu"][name] = [mc, mc_error]
        
        # print for testing
        print "{0}: mc = [{1}, {2}]".format(name, mc, mc_error)

    print "Control Region Units"
    for b in CRunits.binValues[era]:
        name                = CRunit_map[b]
        name                = name.replace("lepcr", "phocr")
        data                = CRunits.binValues[era][b]["data"]
        data_error          = CRunits.binValues[era][b]["data_error"]
        mc_gjets            = CRunits.binValues[era][b]["mc_gjets"]
        mc_gjets_error      = CRunits.binValues[era][b]["mc_gjets_error"]
        mc_back             = CRunits.binValues[era][b]["mc_back"]
        mc_back_error       = CRunits.binValues[era][b]["mc_back_error"]

        # save value and error in map
        map_out["yieldsMap"]["phocr_data"][name]    = [data, data_error]
        map_out["yieldsMap"]["phocr_gjets"][name]   = [mc_gjets, mc_gjets_error]
        map_out["yieldsMap"]["phocr_back"][name]    = [mc_back, mc_back_error]
        
        # print for testing
        print "{0}: data = [{1}, {2}], mc_gjets = [{3}, {4}], mc_back = [{5}, {6}]".format(name, data, data_error, mc_gjets, mc_gjets_error, mc_back, mc_back_error)


    with open (outFile, "w") as f_out:
        json.dump(map_out, f_out, sort_keys=True, indent=4, separators=(',', ' : '))

def run():
    jsonFile = "dc_BkgPred_BinMaps_master.json"
    show(jsonFile, "binNum",    "Search Bins")
    show(jsonFile, "unitCRNum", "Control Region Units")
    show(jsonFile, "unitSRNum", "Search Region Units")

if __name__ == "__main__":
    run()

