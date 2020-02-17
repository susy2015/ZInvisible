# run_modules.py

import argparse
import json
import os
from cutflow_plot import makeCutflows 
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from systematics import Systematic
from search_bins import SearchBins, ValidationBins, ValidationBinsMETStudy, CRUnitBins, SRUnitBins
from make_table import Table
from data_card import makeDataCard
from units import saveResults


def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--json_file",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options     = parser.parse_args()
    json_file   = options.json_file
    verbose     = options.verbose

    doUnits     = True
    doCutflows  = False
    doPhotons   = False
    draw        = False

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    #eras = ["2016"]
    #eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM", "Run2"]
    eras = ["2016", "2017", "2018", "Run2"]
    dirList = []
    plot_dir                = "more_plots"
    latex_dir               = "latex_files"
    results_dir             = "datacard_inputs"
    dirList.append(plot_dir)
    dirList.append(latex_dir)
    dirList.append(results_dir)
    # add "/" to directory if not present
    if plot_dir[-1]  != "/":                plot_dir                += "/"
    if latex_dir[-1] != "/":                latex_dir               += "/"
    if results_dir[-1] != "/":              results_dir             += "/"
    
    for d in dirList:
        # make directory if it does not exist
        if not os.path.exists(d):
            os.makedirs(d)
    # normalization
    N = Normalization(plot_dir, verbose)
    # shape
    S = Shape(plot_dir, draw, doUnits, verbose)
    # validation bins
    VB    = ValidationBins(            N, S, eras, plot_dir, verbose, draw=True, saveRootFile=True )
    VB_MS = ValidationBinsMETStudy(    N, S, eras, plot_dir, verbose, draw=True, saveRootFile=True )
    # search bins
    SB = SearchBins(        N, S, eras, plot_dir, verbose, draw=True, saveRootFile=True )
    if doUnits:
        # control region unit bins  
        CRunits = CRUnitBins(N, S, eras, plot_dir, verbose) 
        # search region unit bins  
        SRunits = SRUnitBins(N, S, eras, plot_dir, verbose) 
    # systematics
    Syst = Systematic(plot_dir, N, S)
    
    with open(json_file, "r") as input_file:
        runMap = json.load(input_file)
        # loop over eras
        for era in eras:
            print "|---------------------------------------------------- Era: {0} ----------------------------------------------------|".format(era)
            runDir = runMap[era]
            result_file = "condor/" + runDir + "/result.root"
            if doCutflows:
                makeCutflows(result_file, era, plot_dir, doPhotons)
            N.getNormAndError(result_file, era)
            S.getShape(result_file, era)
            VB.getValues(result_file, era)
            VB_MS.getValues(result_file, era)
            SB.getValues(result_file, era)
            if doUnits:
                CRunits.getValues(result_file, era)
                SRunits.getValues(result_file, era)
            # WARNING: var, varPhoton, and varLepton are used to load histograms and must match histogram names
            # makeZvsPhoton(self, file_name, var, varPhoton, varLepton, era, rebin, useForSyst)
            Syst.makeZvsPhoton(result_file, "met",  "metWithPhoton",                            "metWithLL",                                era, True, True)
            Syst.makeZvsPhoton(result_file, "ht",   "HT_drPhotonCleaned_jetpt30",               "HT_drLeptonCleaned_jetpt30",               era, False, False)
            Syst.makeZvsPhoton(result_file, "nj",   "nJets_drPhotonCleaned_jetpt30",            "nJets_drLeptonCleaned_jetpt30",            era, False, False)
            Syst.makeZvsPhoton(result_file, "nb",   "nBottoms_drPhotonCleaned_jetpt30",         "nBottoms_drLeptonCleaned_jetpt30",         era, False, False)
            Syst.makeZvsPhoton(result_file, "nmt",  "nMergedTops_drPhotonCleaned_jetpt30",      "nMergedTops_drLeptonCleaned_jetpt30",      era, False, False)
            Syst.makeZvsPhoton(result_file, "nw",   "nWs_drPhotonCleaned_jetpt30",              "nWs_drLeptonCleaned_jetpt30",              era, False, False)
            Syst.makeZvsPhoton(result_file, "nrt",  "nResolvedTops_drPhotonCleaned_jetpt30",    "nResolvedTops_drLeptonCleaned_jetpt30",    era, False, False)

    # Normalization: makeTable(self, output_name, makeDoc=False)
    N.makeTable(latex_dir + "zinv_rz_doc.tex",   True)
    N.makeTable(latex_dir + "zinv_rz_table.tex", False)
    N.makeTexFile("validation", latex_dir + "validationBins_normalization_Zmass.tex")
    N.makeTexFile("search",     latex_dir + "searchBins_normalization_Zmass.tex")
    N.makeComparison("validation")
    N.makeComparison("validationMetStudy")
    N.makeComparison("search")
    S.makeComparison("validation")
    S.makeComparison("validationMetStudy")
    S.makeComparison("search")
    VB.makeTexFile("Z Invisible Per Era Prediction for Validation Bins", latex_dir + "zinv_per_era_prediction_validation_bins.tex")
    SB.makeTexFile("Z Invisible Per Era Prediction for Search Bins",     latex_dir + "zinv_per_era_prediction_search_bins.tex")
    
    # total Run2 prediction
    total_era = "Run2"
    VB.makeTexFile("Z Invisible Total Prediction for Validation Bins", latex_dir + "zinv_total_prediction_validation_bins.tex", total_era)
    SB.makeTexFile("Z Invisible Total Prediction for Search Bins",     latex_dir + "zinv_total_prediction_search_bins.tex",     total_era)
    T = Table()
    # fancy table only supported in search bins right now
    # makeYieldTable(self, BinObject, total_era, output="pred_sr.tex", makeDoc=False, size=0.6)
    T.makeYieldTable(SB, total_era, latex_dir + "zinv_pred_sr_doc.tex",   True,  0.55)
    T.makeYieldTable(SB, total_era, latex_dir + "zinv_pred_sr_table.tex", False, 0.60)
    
    # Get systematics in proper bins: Rz and "Z to LL vs. Photon" systematics
    # must be done after N.makeComparison()
    # must be done after Syst.makeZvsPhoton()
    VB.getRzSyst(               N.rz_syst_map,      "validation",           "RzSyst_ValidationBins.root")
    VB_MS.getRzSyst(            N.rz_syst_map,      "validationMetStudy",   "RzSyst_ValidationBinsMETStudy.root")
    SB.getRzSyst(               N.rz_syst_map,      "search",               "RzSyst_SearchBins.root")
    VB.getZvsPhotonSyst(        Syst.h_map_syst,                            "ZvsPhotonSyst_ValidationBins.root")
    VB_MS.getZvsPhotonSyst(     Syst.h_map_syst,                            "ZvsPhotonSyst_ValidationBinsMETStudy.root")
    SB.getZvsPhotonSyst(        Syst.h_map_syst,                            "ZvsPhotonSyst_SearchBins.root")
    CRunits.getZvsPhotonSyst(   Syst.h_map_syst,                            "ZvsPhotonSyst_CRUnitBins.root")

    if doUnits:
        # saveResults(inFile, outFile, CRunits, SRunits, SB, era)
        saveResults("dc_BkgPred_BinMaps_master.json", results_dir + "zinv_yields_" + total_era + ".json", CRunits, SRunits, SB, total_era)


if __name__ == "__main__":
    main()


