# run_modules.py

import numpy as np
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

    doRun2      = True
    doUnits     = True
    doCutflows  = False
    doPhotons   = False
    draw        = False

    if not os.path.exists(json_file):
        print "The json file \"{0}\" containing runs does not exist.".format(json_file)
        return
    
    eras = ["2016", "2017", "2018", "2016and2017", "Run2"]
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
        CRunits = CRUnitBins(N, S, eras, plot_dir, verbose, saveRootFile=True) 
        # search region unit bins  
        SRunits = SRUnitBins(N, S, eras, plot_dir, verbose, saveRootFile=True) 
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
            if doUnits:
                CRunits.getValues(result_file, era)
                SRunits.getValues(result_file, era)
            VB.getValues(result_file, era)
            VB_MS.getValues(result_file, era)
            SB.getValues(result_file, era, CRunits=CRunits)
            xbins_met   = np.array([0.0, 250.0, 350.0, 450.0, 550.0, 650.0, 1000.0])
            n_bins_met  = len(xbins_met) - 1
            xbins_ht    = np.array([0.0, 300.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0])
            n_bins_ht   = len(xbins_ht) - 1
            # WARNING: var, varPhoton, and varLepton are used to load histograms and must match histogram names
            # makeZvsPhoton(self, file_name, var, varPhoton, varLepton, era, x_min, x_max, n_bins=0, xbins=np.array([]), rebin=False, useForSyst=False, doDataOverData=False)
            Syst.makeZvsPhoton(result_file, "met",  "metWithPhoton",                            "metWithLL",                                era, 250.0, 1000.0, n_bins_met, xbins_met, True,  True,  False )
            Syst.makeZvsPhoton(result_file, "met",  "metWithPhoton",                            "metWithLL",                                era, 250.0, 1000.0, n_bins_met, xbins_met, True,  False, True  )
            Syst.makeZvsPhoton(result_file, "ht",   "HT_drPhotonCleaned_jetpt30",               "HT_drLeptonCleaned_jetpt30",               era, 300.0, 2000.0, n_bins_ht,  xbins_ht,  True,  False, False )
            Syst.makeZvsPhoton(result_file, "ht",   "HT_drPhotonCleaned_jetpt30",               "HT_drLeptonCleaned_jetpt30",               era, 300.0, 2000.0, n_bins_ht,  xbins_ht,  True,  False, True  )
            Syst.makeZvsPhoton(result_file, "nj",   "nJets_drPhotonCleaned_jetpt30",            "nJets_drLeptonCleaned_jetpt30",            era, 0, 11)
            Syst.makeZvsPhoton(result_file, "nb",   "nBottoms_drPhotonCleaned_jetpt30",         "nBottoms_drLeptonCleaned_jetpt30",         era, 0,  6)
            Syst.makeZvsPhoton(result_file, "nmt",  "nMergedTops_drPhotonCleaned_jetpt30",      "nMergedTops_drLeptonCleaned_jetpt30",      era, 0,  6)
            Syst.makeZvsPhoton(result_file, "nw",   "nWs_drPhotonCleaned_jetpt30",              "nWs_drLeptonCleaned_jetpt30",              era, 0,  6)
            Syst.makeZvsPhoton(result_file, "nrt",  "nResolvedTops_drPhotonCleaned_jetpt30",    "nResolvedTops_drLeptonCleaned_jetpt30",    era, 0,  6)
            # studyShapes(file_name, region, era, variable, nameTag, varName, xbins, n_bins)
            S.studyShapes(result_file, "LowDM",  era, "metWithPhoton",              "_met_split_jetpt30", "MET", xbins_met, n_bins_met)
            S.studyShapes(result_file, "HighDM", era, "metWithPhoton",              "_met_split_jetpt30", "MET", xbins_met, n_bins_met)
            S.studyShapes(result_file, "LowDM",  era, "HT_drPhotonCleaned_jetpt30", "_ht_split_jetpt30",  "HT",  xbins_ht,  n_bins_ht)
            S.studyShapes(result_file, "HighDM", era, "HT_drPhotonCleaned_jetpt30", "_ht_split_jetpt30",  "HT",  xbins_ht,  n_bins_ht)

    
    # Normalization: makeTable(self, output_name, makeDoc=False)
    if doRun2:
        N.makeComparison("validation")
        N.makeComparison("validationMetStudy")
        N.makeComparison("search")
        # N.makeTable to be after N.makeComparison
        N.makeTable(latex_dir + "zinv_rz_doc.tex",   True)
        N.makeTable(latex_dir + "zinv_rz_table.tex", False)
        N.makeTexFile("validation", latex_dir + "validationBins_normalization_Zmass.tex")
        N.makeTexFile("search",     latex_dir + "searchBins_normalization_Zmass.tex")
        S.makeComparison("validation")
        S.makeComparison("validationMetStudy")
        S.makeComparison("search")
        S.makeTable(latex_dir + "zinv_q_doc.tex",   True)
        S.makeTable(latex_dir + "zinv_q_table.tex", False)
        VB.makeTexFile("Z Invisible Per Era Prediction for Validation Bins", latex_dir + "zinv_per_era_prediction_validation_bins.tex")
        SB.makeTexFile("Z Invisible Per Era Prediction for Search Bins",     latex_dir + "zinv_per_era_prediction_search_bins.tex")
    
        # total era is combiniation of all eras
        total_era = "Run2"
        # total Run2 prediction
        VB.makeTexFile("Z Invisible Total Prediction for Validation Bins", latex_dir + "zinv_total_prediction_validation_bins.tex", total_era)
        SB.makeTexFile("Z Invisible Total Prediction for Search Bins",     latex_dir + "zinv_total_prediction_search_bins.tex",     total_era)
    
    T = Table()
    
    # Get systematics in proper bins: Rz and "Z to LL vs. Photon" systematics
    # must be done after N.makeComparison()
    # must be done after Syst.makeZvsPhoton()
    if doRun2:
        VB.getRzSyst(               N.rz_syst_map,      "validation",           "RzSyst_ValidationBins.root")
        VB_MS.getRzSyst(            N.rz_syst_map,      "validationMetStudy",   "RzSyst_ValidationBinsMETStudy.root")
        SB.getRzSyst(               N.rz_syst_map,      "search",               "RzSyst_SearchBins.root")
        VB.getZvsPhotonSyst(        Syst.h_map_syst,                            "ZvsPhotonSyst_ValidationBins.root")
        VB_MS.getZvsPhotonSyst(     Syst.h_map_syst,                            "ZvsPhotonSyst_ValidationBinsMETStudy.root")
        SB.getZvsPhotonSyst(        Syst.h_map_syst,                            "ZvsPhotonSyst_SearchBins.root")
        CRunits.getZvsPhotonSyst(   Syst.h_map_syst,                            "ZvsPhotonSyst_CRUnitBins.root")
            
    # save maps to json files
    # makeJson(self, map_, file_)
    VB.makeJson(VB.binValues,           "results/ValidationBinResults.json")
    SB.makeJson(SB.binValues,           "results/SearchBinResults.json")

    if doUnits:
        # save maps to json files
        # makeJson(self, map_, file_)
        CRunits.makeJson(CRunits.binValues, "results/CRUnitsResults.json")
        SRunits.makeJson(SRunits.binValues, "results/SRUnitsResults.json")
        for era in eras:
            # save yields for data card inputs
            # saveResults(inFile, outFile, CRunits, SRunits, SB, era)
            saveResults("dc_BkgPred_BinMaps_master.json", results_dir + "zinv_yields_" + era + ".json", CRunits, SRunits, SB, era)
            # fancy table only supported in search bins right now
            # makeYieldTable(self, BinObject, total_era, output="pred_sr.tex", makeDoc=False, size=1.0)
            T.makeYieldTable(SB, era, latex_dir + "zinv_pred_sr_" + era + "_doc.tex",   True )
            T.makeYieldTable(SB, era, latex_dir + "zinv_pred_sr_" + era + "_table.tex", False)


if __name__ == "__main__":
    main()


