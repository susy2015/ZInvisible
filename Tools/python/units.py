# units.py
# here we go.....
# use those CR and SR units!

import json
from search_bins import SearchBins, ValidationBins, CRUnitBins, SRUnitBins
from tools import invert, getConstantMultiplicationError


def show(binMap, title):
    # invert map
    invMap = invert(binMap)
    nBins = len(invMap)
    bins = list(str(x) for x in xrange(nBins))
    print "{0}: {1} bins".format(title, nBins)
    for b in bins: 
        print "{0} : {1}".format(b, invMap[b])

# save final results in json file
def saveResults(inFile, outFile, CRunits, SRunits, SB, era):
    verbose = False
    map_out                                 = {}
    map_out["yieldsMap"]                    = {}
    map_out["yieldsMap"]["znunu"]           = {}
    map_out["yieldsMap"]["phocr_data"]      = {}
    map_out["yieldsMap"]["phocr_gjets"]     = {}
    map_out["yieldsMap"]["phocr_back"]      = {}
    unitMap                                 = {}
    normForSRunits                          = {}
    normForCRunits                          = {}
    with open (inFile, "r") as f_in:
        map_in = json.load(f_in)
        # invert maps
        SearchBin_map   = invert(map_in["binNum"])
        CRunit_map      = invert(map_in["unitCRNum"]["phocr"])
        SRunit_map      = invert(map_in["unitSRNum"])

    with open("units.json", "r") as input_file:
        unitMap = json.load(input_file)
    
    # final yields for json file
    # znunu:        (Z nunu MC in SR unit bins) * Rz
    # phocr_data:   Photon Data in photon CR 
    # phocr_gjets:  GJets MC in photon CR normalized to Data
    # phocr_back:   Other MC (QCD, TTGJets, etc) in photon CR normalized to Data

    # ------------------------------------------------------ #
    # get Rz values in search bins and apply to SR unit bins
    # - Rz norm in search bins: SB.binValues[era][b]["norm"] 
    # get Photon Data/MC norm. and apply to CR unit bins
    # - photon norm in search bins: SB.binValues[era][b]["photon_data_mc_norm"] 
    # ------------------------------------------------------ #
    for b in SB.binValues[era]:
        norm             = SB.binValues[era][b]["norm"]
        photonDataMCNorm = SB.binValues[era][b]["photon_data_mc_norm"] 
        sr_units = unitMap["unitBinMapSR"][b]
        cr_units = unitMap["unitBinMapCR_phocr"][b]
        for sr in sr_units:
            normForSRunits[sr] = {}
            normForSRunits[sr]["norm"] = norm
        for cr in cr_units:
            normForCRunits[cr] = {}
            normForCRunits[cr]["norm"] = photonDataMCNorm

    
    if verbose:
        print " --- Search Region Units --- "
    for b in SRunits.binValues[era]:
        name        = SRunit_map[b]
        mc          = SRunits.binValues[era][b]["mc"]
        mc_error    = SRunits.binValues[era][b]["mc_error"]
        # normazliation for SR unit bins
        norm = -999
        try:
            norm = normForSRunits[b]["norm"]
        except:
            print "ERROR: No normalization value found for SR unit bin {0}".format(b)
        # do not propagate norm to mc_error as total Rz unc. (stat and syst) are included in Rz syst.
        mc_with_norm = norm * mc
        # save value and error
        map_out["yieldsMap"]["znunu"][name] = [mc_with_norm, mc_error]
        
        if verbose:
            print "{0}: mc = {1}, mc_with_norm = {2}, mc_error = {3}".format(name, mc, mc_with_norm, mc_error)

    if verbose:
        print " --- Control Region Units --- "
    for b in CRunits.binValues[era]:
        name                = CRunit_map[b]
        data                = CRunits.binValues[era][b]["data"]
        data_error          = CRunits.binValues[era][b]["data_error"]
        mc_gjets            = CRunits.binValues[era][b]["mc_gjets"]
        mc_gjets_error      = CRunits.binValues[era][b]["mc_gjets_error"]
        mc_back             = CRunits.binValues[era][b]["mc_back"]
        mc_back_error       = CRunits.binValues[era][b]["mc_back_error"]
        # normazliation for CR unit bins
        norm = -999
        try:
            norm = normForCRunits[b]["norm"]
        except:
            print "ERROR: No normalization value found for CR unit bin {0}".format(b)
        # propagate norm to photon MC stat. error
        mc_gjets_with_norm          = norm * mc_gjets
        mc_back_with_norm           = norm * mc_back
        mc_gjets_error_with_norm    = getConstantMultiplicationError(norm, mc_gjets_error)
        mc_back_error_with_norm     = getConstantMultiplicationError(norm, mc_back_error)
        # save value and error in map
        map_out["yieldsMap"]["phocr_data"][name]    = [data,                    data_error]
        map_out["yieldsMap"]["phocr_gjets"][name]   = [mc_gjets_with_norm,      mc_gjets_error_with_norm]
        map_out["yieldsMap"]["phocr_back"][name]    = [mc_back_with_norm,       mc_back_error_with_norm]
        
        if verbose:
            print "{0}: data = [{1}, {2}], mc_gjets = [{3}, {4}, {5}], mc_back = [{6}, {7}, {8}]".format(name, data, data_error, mc_gjets, mc_gjets_with_norm, mc_gjets_error, mc_back, mc_back_with_norm, mc_back_error)


    with open (outFile, "w") as f_out:
        json.dump(map_out, f_out, sort_keys=True, indent=4, separators=(',', ' : '))

def run():
    jsonFile = "dc_BkgPred_BinMaps_master.json"
    with open(jsonFile, "r") as inputFile:
        binMap = json.load(inputFile)
        show(binMap["binNum"],              "Search Bins")
        show(binMap["unitCRNum"]["phocr"],  "Control Region Units")
        show(binMap["unitSRNum"],           "Search Region Units")

if __name__ == "__main__":
    run()

