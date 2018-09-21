# calculateEfficiency.py
# Caleb J. Smith
# September 19, 2018

# calculate efficiencies and purities

import ROOT
import argparse
from decimal import Decimal

# get results given input file and list of histograms
def getResults(tag, inputFile, histograms):
    f = ROOT.TFile(inputFile)
    if not f:
        print "ERROR: The file {0} did not load".format(inputFile)
        return
    nEventsPerHist = []
    for histogramMap in histograms:
        #print histogramMap
        d_name, h_name = histogramMap.items()[0]
        h = f.Get(d_name + "/" + h_name)
        #print "directory: {0} ; h: {1}".format(d_name, h)
        if not h:
            print "ERROR: The histogram {0} from directory {1} did not load".format(h_name, d_name)
            return
        #nEventsPerHist.append({d_name : h.GetEntries()})
        nEventsPerHist.append({d_name : h.Integral()})
    #print nEventsPerHist
    return calcEff(tag, nEventsPerHist)

# calculate and print ratios
def calcEff(tag, variables):
    results = {}
    i = 0
    total_efficiency = 1.0
    #print variables
    while i < len(variables) - 1: 
        #print variables[i]
        # get key (name) and value (value) 
        a_name, a_value = variables[i].items()[0]
        b_name, b_value = variables[i+1].items()[0]
        ratio = float(b_value) / float(a_value)
        ratio_percent = 100.0 * ratio
        total_efficiency *= ratio
        names = "{0}/{1}".format(b_name, a_name)
        values = "{0:.3E}/{1:.3E}".format(Decimal(b_value), Decimal(a_value))
        results[names] = ratio
        print("{0:50} = {1:20} = {2:10.4f}: {3:10.2f} %".format(names, values, ratio, ratio_percent))
        i += 1
    total_efficiency_percent = 100.0 * total_efficiency
    print "Total {0} Efficiency = {1:.4f}: {2:.2f} %".format(tag, total_efficiency, total_efficiency_percent)
    return results

def calcFinal(variables):
    output = ""
    total = 1.0
    for i, variable in enumerate(variables):
        output += "{0} ".format(variable)
        if i < len(variables) - 1:
            output += "* "
        total *= variables[variable]
    total_percent = 100.0 * total
    output += "= {0:.4f}: {1:.2f} %".format(total, total_percent)
    print output
    return total


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_file", "-f", default="result.root", help="input root file")
    options = parser.parse_args()
    
    # gen
    histogramMapGen = []
    histogramMapGen += [{"gammaLVecGen"                : "MC_PhotonGenAccPt_singlegammaLVecGenptgammaLVecGen(pt)GJets #gamma Gensingle"}]
    histogramMapGen += [{"gammaLVecGenEta"             : "MC_PhotonGenAccPt_singlegammaLVecGenEtaptgammaLVecGenEta(pt)GJets #gamma GenEtasingle"}]
    histogramMapGen += [{"gammaLVecGenEtaPt"           : "MC_PhotonGenMatchPt_singlegammaLVecGenEtaPtptgammaLVecGenEtaPt(pt)GJets #gamma GenEtaPtsingle"}]
    histogramMapGen += [{"gammaLVecGenEtaPtMatched"    : "MC_PhotonGenMatchPt_singlegammaLVecGenEtaPtMatchedptgammaLVecGenEtaPtMatched(pt)GJets #gamma GenEtaPtMatchedsingle"}]
    
    # reco
    histogramMapReco = []
    histogramMapReco += [{"gammaLVecReco"                : "MC_PhotonRecoAccPt_singlegammaLVecRecoptgammaLVecReco(pt)GJets #gamma Recosingle"}]
    histogramMapReco += [{"gammaLVecRecoEta"             : "MC_PhotonRecoAccPt_singlegammaLVecRecoEtaptgammaLVecRecoEta(pt)GJets #gamma RecoEtasingle"}]
    histogramMapReco += [{"gammaLVecRecoEtaPt"           : "MC_PhotonRecoIsoPt_singlegammaLVecRecoEtaPtptgammaLVecRecoEtaPt(pt)GJets #gamma RecoEtaPtsingle"}]
    histogramMapReco += [{"gammaLVecRecoIso"             : "MC_PhotonRecoIsoPt_singlegammaLVecRecoIsoptgammaLVecRecoIso(pt)GJets #gamma RecoIsosingle"}]
    histogramMapReco += [{"gammaLVecRecoEtaPtMatched"    : "MC_PhotonRecoMatchPt_singlegammaLVecRecoEtaPtMatchedptgammaLVecRecoEtaPtMatched(pt)GJets #gamma RecoEtaPtMatchedsingle"}]
   

    resultsGen  = getResults("Gen",  options.input_file, histogramMapGen)
    resultsReco = getResults("Reco", options.input_file, histogramMapReco)

    # final
    finalGen = {}
    finalGen["genAcc"]           = resultsGen["gammaLVecGenEta/gammaLVecGen"]
    finalGen["genMatchedToReco"] = resultsGen["gammaLVecGenEtaPtMatched/gammaLVecGenEtaPt"]
    genEff = calcFinal(finalGen)

    finalReco = {}
    finalReco["recoAcc"]          = resultsReco["gammaLVecRecoEta/gammaLVecReco"]
    finalReco["recoIso"]          = resultsReco["gammaLVecRecoIso/gammaLVecRecoEtaPt"]
    finalReco["recoMatchedToGen"] = resultsReco["gammaLVecRecoEtaPtMatched/gammaLVecRecoIso"]
    recoEff = calcFinal(finalReco)
    

    # demo
    #calcEff([ 
    #        {"reco":            1 * 10 ** 6}, 
    #        {"eta_cut":         1 * 10 ** 5},
    #        {"eta_pt_cut":      1 * 10 ** 4},
    #        {"gen_matched":     1 * 10 ** 3},
    #        {"loose_isolation": 1 * 10 ** 2}
    #       ])





