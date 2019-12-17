# Run_modules.py (modified for systematics)
import copy
import numpy as np
import json
import argparse
import os
import ROOT
import run_systematics
from Norm_lepton_zmass import Normalization
from Shape_photon_met import Shape
from Search_bins import  ValidationBins, SearchBins, CRUnitBins
from tools import invert

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


#TODO: combined all systematics to get total syst up/down
# first combined syst in validation bins and give to Hui
# loop over validation bins
# loop over systematics
# syst_up_i       = (p_up - p)   / p
# syst_down_i     = (p - p_down) / p
# syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
# syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 

# use histogram which stores systematic errors
def writeToConfFromSyst(outFile, searchBinMap, syst, h, region, offset):
    zinv = "znunu"
    # sb_i = bin_i - 1 + offset
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = searchBinMap[str(sb_i)]
        # systematic error is stored in bin content
        # treat error as symmetric
        error   = h.GetBinContent(i)
        r_up   = 1 + error 
        r_down = 1 - error 

        outFile.write("{0}  {1}_Up    {2}  {3}\n".format( sb_name, syst, zinv, r_up   ) )
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, syst, zinv, r_down ) )

# use prediction, syst up/down histograms
def writeToConfFromPred(outFile, searchBinMap, syst, h, h_up, h_down, region, offset):
    zinv = "znunu"
    # sb_i = bin_i - 1 + offset
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = searchBinMap[str(sb_i)]
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        #s_up   = (p_up - p)   / p
        #s_down = (p - p_down) / p
        # If there is a zero prediction for nominal, up, and down can you set the systematic to 1.
        # The error is a deviation from 1.
        r_up   = 1
        r_down = 1
        # do not apply SB syst in high dm
        # here syst is already systForConf (systForConf = systMap[syst]["name"])
        if not (region == "highdm" and syst == "ivfunc"):
            if p != 0:
                r_up   = p_up   / p
                r_down = p_down / p
            else:
                print "WARNING: pred = 0 for search bin {0}".format(sb_i)
        outFile.write("{0}  {1}_Up    {2}  {3}\n".format( sb_name, syst, zinv, r_up   ) )
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, syst, zinv, r_down ) )


def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--runs_json",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--syst_json",    "-s", default="",                             help="json file containing systematics")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options         = parser.parse_args()
    runs_json       = options.runs_json
    syst_json       = options.syst_json
    verbose         = options.verbose
    units_json      = "dc_BkgPred_BinMaps_master.json"
    rz_syst_files   = {}
    rz_syst_files["validation"]   = "RzSyst_ValidationBins.root"
    rz_syst_files["search"]       = "RzSyst_SearchBins.root"
    rz_syst_files["controlUnit"]  = "RzSyst_SearchBins.root" # "RzSyst_ControlUnitBins.root" ?????
    ZvPhoton_syst_files = {}
    ZvPhoton_syst_files["validation"]   = "ZvsPhotonSyst_ValidationBins.root"
    ZvPhoton_syst_files["search"]       = "ZvsPhotonSyst_SearchBins.root" 
    ZvPhoton_syst_files["controlUnit"]  = "ZvsPhotonSyst_SearchBins.root" #"ZvsPhotonSyst_ControlUnitBins.root" ????
   
    doUnits      = True
    draw         = False
    runMap       = {}
    systMap      = {}
    unitMap      = {}
    systHistoMap = {}

    # list required files to check if they exist
    file_list = [] 
    file_list.append(runs_json)
    file_list.append(syst_json)
    file_list.append(units_json)
    for key in rz_syst_files:
        file_list.append(rz_syst_files[key])
    for key in ZvPhoton_syst_files:
        file_list.append(ZvPhoton_syst_files[key])

    # check if files exist
    for file_name in file_list:
        if not os.path.exists(file_name):
            print "ERROR: The required file \"{0}\" does not exist.".format(file_name)
            return
    
    # load json files
    with open(runs_json, "r") as input_file:
        runMap  = json.load(input_file)
    with open(syst_json, "r") as input_file:
        systMap = json.load(input_file) 
    with open(units_json, "r") as input_file:
        unitMap = json.load(input_file) 
    
    # map search bin numbers to string names
    searchBinMap = invert(unitMap["binNum"])

    era          = "Run2"
    runDir = runMap[era]
    result_file  = "condor/" + runDir + "/result.root"
    conf_file    = "results/zinv_syst_" + era + ".conf"
    out_dir      = "caleb_syst_plots/" 
    syst_dir     = "caleb_syst_plots/"
    variable     = "mc"
    regions      = ["lowdm", "highdm"]
    directions   = ["up", "", "down"]
    bintypes     = ["validation", "search", "controlUnit"]
    systematics  = ["jes","eff_restoptag_photon","eff_sb_photon","btag","eff_sb","eff_toptag","eff_toptag_photon","eff_wtag","met_trig","pileup","photon_trig","prefire","photon_sf","eff_wtag_photon"] # all systematics available from Znunu_nValidationBin and nCRUnits
    # missing syst: prefire weight (2016, 2017), pdf_weight, MET resolution (uncluster), lepton veto SF, ISR_Weight

    # histo_tmp[region][direction]
    histo_tmp  = {region:dict.fromkeys(directions) for region in regions}

    # histo[bintype][region][direction]
    histo      = {bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} 

    # syst_histo[systemaitc][bintype][region][direction]
    syst_histo = {syst:{bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} for syst in systematics} 

    #-------------------------------------------------------
    # Class instanceses summoning 
    #-------------------------------------------------------
    
    N   =  Normalization(out_dir, verbose)
    S   =  Shape(out_dir, draw, doUnits, verbose)
    VB  =  ValidationBins(N, S, era, out_dir, verbose)
    SB  =  SearchBins(N, S, era, out_dir, verbose)
    CRU =  CRUnitBins(N, S, era, out_dir, verbose)
    
    #-------------------------------------------------------
    # Normal predictions (no systematics) 
    #-------------------------------------------------------
    
    N.getNormAndError(result_file, "", era)
    S.getShape(result_file, "", era)
    VB.getValues(result_file, "", era)
    SB.getValues(result_file, "", era)
    CRU.getValues(result_file, "", era)
    
    histo["validation"]["lowdm"][""]     =  VB.histograms[era]["lowdm"][variable].Clone()
    histo["validation"]["highdm"][""]    =  VB.histograms[era]["highdm"][variable].Clone()
    histo["search"]["lowdm"][""]         =  SB.histograms[era]["lowdm"][variable].Clone()
    histo["search"]["highdm"][""]        =  SB.histograms[era]["highdm"][variable].Clone()
    print("debbuging ")
    print(CRU.histograms.keys())
    histo["controlUnit"]["highdm"][""]   =  CRU.histograms[era]["highdm"]["mc_gjets"].Clone()
    histo["controlUnit"]["lowdm"][""]    =  CRU.histograms[era]["lowdm"]["mc_gjets"].Clone()
    
    #-------------------------------------------------------
    # Calculate normalization and shape factors
    #-------------------------------------------------------
    with open(conf_file, "w") as outFile:
           
        # --- begin loop over systematics
        for syst in systematics:
        
            # --- syst up --- #
            N.getNormAndError(result_file, syst + "_syst_up", era)
            S.getShape(result_file, syst + "_syst_up", era)
            VB.getValues(result_file, syst + "_syst_up", era)
            SB.getValues(result_file, syst + "_syst_up", era)
            CRU.getValues(result_file, syst + "_syst_up", era)
            
            histo["validation"]["highdm"]["up"]   =  VB.histograms[era]["highdm"][variable].Clone()
            histo["validation"]["lowdm"]["up"]    =  VB.histograms[era]["lowdm"][variable].Clone()
            histo["search"]["highdm"]["up"]       =  SB.histograms[era]["highdm"][variable].Clone()
            histo["search"]["lowdm"]["up"]        =  SB.histograms[era]["lowdm"][variable].Clone()
            histo["controlUnit"]["highdm"]["up"]  =  CRU.histograms[era]["highdm"]["mc_gjets"].Clone()
            histo["controlUnit"]["lowdm"]["up"]   =  CRU.histograms[era]["lowdm"]["mc_gjets"].Clone()
            syst_histo[syst]["validation"]["lowdm"]["up"]    = histo["validation"]["lowdm"]["up"]
            syst_histo[syst]["validation"]["highdm"]["up"]   = histo["validation"]["highdm"]["up"]
            syst_histo[syst]["search"]["lowdm"]["up"]        = histo["search"]["lowdm"]["up"]
            syst_histo[syst]["search"]["highdm"]["up"]       = histo["search"]["highdm"]["up"]
            syst_histo[syst]["controlUnit"]["lowdm"]["up"]   = histo["controlUnit"]["lowdm"]["up"]
            syst_histo[syst]["controlUnit"]["highdm"]["up"]  = histo["controlUnit"]["highdm"]["up"]
            
            # --- syst down --- #
            N.getNormAndError(result_file, syst + "_syst_down", era)
            S.getShape(result_file, syst + "_syst_down", era)
            VB.getValues(result_file, syst + "_syst_down", era)
            SB.getValues(result_file, syst + "_syst_down", era)
            CRU.getValues(result_file, syst + "_syst_down", era)
            
            histo["validation"]["highdm"]["down"]   =  VB.histograms[era]["highdm"][variable].Clone()
            histo["validation"]["lowdm"]["down"]    =  VB.histograms[era]["lowdm"][variable].Clone()
            histo["search"]["highdm"]["down"]       =  SB.histograms[era]["highdm"][variable].Clone()
            histo["search"]["lowdm"]["down"]        =  SB.histograms[era]["lowdm"][variable].Clone()
            histo["controlUnit"]["highdm"]["down"]  =  CRU.histograms[era]["highdm"]["mc_gjets"].Clone()
            histo["controlUnit"]["lowdm"]["down"]   =  CRU.histograms[era]["lowdm"]["mc_gjets"].Clone()
            syst_histo[syst]["validation"]["lowdm"]["down"]    =  histo["validation"]["lowdm"]["down"]
            syst_histo[syst]["validation"]["highdm"]["down"]   =  histo["validation"]["highdm"]["down"]
            syst_histo[syst]["search"]["lowdm"]["down"]        =  histo["search"]["lowdm"]["down"]
            syst_histo[syst]["search"]["highdm"]["down"]       =  histo["search"]["highdm"]["down"]
            syst_histo[syst]["controlUnit"]["lowdm"]["down"]   =  histo["controlUnit"]["lowdm"]["down"]
            syst_histo[syst]["controlUnit"]["highdm"]["down"]  =  histo["controlUnit"]["highdm"]["down"]
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------
            # only write search bin systematics to conf
            # TODO: write CR unit systematics to conf
            #writeToConfFromPred(outFile, searchBinMap, syst, h, h_up, h_down, region, offset)
            systForConf = systMap[syst]["name"]  
            writeToConfFromPred(outFile, searchBinMap, systForConf, histo["search"]["lowdm"][""],  histo["search"]["lowdm"]["up"],  histo["search"]["lowdm"]["down"],  "lowdm",  SB.low_dm_start)
            writeToConfFromPred(outFile, searchBinMap, systForConf, histo["search"]["highdm"][""], histo["search"]["highdm"]["up"], histo["search"]["highdm"]["down"], "highdm", SB.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------
            
            for bintype in bintypes:
                for region in regions:
                    # run_systematics.plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir)
                    run_systematics.plot(histo[bintype][region][""], histo[bintype][region]["up"], histo[bintype][region]["down"], syst, bintype, region, era, syst_dir)
                
        # --- end loop over systematics
            
        #-------------------------------------------------------
        # Write to conf (for syst. stored in root files)
        #-------------------------------------------------------
        
        # systematics from syst. histos stored in root files 
        # be careful about using bin_type = validation, search as needed
        # also, you may need to copy histos if they are deleted when file is closed
        
        #systHistoMap[bintype][region][syst]
        # put these histos in a map for later use
        for bintype in bintypes:
            systHistoMap[bintype] = {}
            systHistoMap[bintype]["lowdm"] = {}
            systHistoMap[bintype]["highdm"] = {}
            # --- Rz syst --- #
            f_in = ROOT.TFile(rz_syst_files[bintype], "read")
            # histogram names
            # rz_syst_low_dm
            # rz_syst_high_dm
            systHistoMap[bintype]["lowdm"]["znunu_rzunc"]  = copy.deepcopy(f_in.Get("rz_syst_low_dm"))
            systHistoMap[bintype]["highdm"]["znunu_rzunc"] = copy.deepcopy(f_in.Get("rz_syst_high_dm"))
            # --- Z vs Photon syst --- #
            f_in = ROOT.TFile(ZvPhoton_syst_files[bintype], "read")
            # histogram names
            # ZvsPhoton_syst_low_dm
            # ZvsPhoton_syst_high_dm
            systHistoMap[bintype]["lowdm"]["znunu_zgammdiff"]   = copy.deepcopy(f_in.Get("ZvsPhoton_syst_low_dm"))
            systHistoMap[bintype]["highdm"]["znunu_zgammdiff"]  = copy.deepcopy(f_in.Get("ZvsPhoton_syst_high_dm"))

        # --- Rz syst --- #
        # writeToConfFromSyst(outFile, searchBinMap, syst, h, region, offset)
        systForConf = systMap["znunu_rzunc"]["name"]  
        writeToConfFromSyst(outFile, searchBinMap, systForConf, systHistoMap["search"]["lowdm"]["znunu_rzunc"],  "lowdm",  SB.low_dm_start)
        writeToConfFromSyst(outFile, searchBinMap, systForConf, systHistoMap["search"]["highdm"]["znunu_rzunc"], "highdm", SB.high_dm_start)
        
        # --- Z vs Photon syst --- #
        # writeToConfFromSyst(outFile, searchBinMap, syst, h, region, offset)
        systForConf = systMap["znunu_zgammdiff"]["name"]  
        writeToConfFromSyst(outFile, searchBinMap, systForConf, systHistoMap["search"]["lowdm"]["znunu_zgammdiff"],  "lowdm",  SB.low_dm_start)
        writeToConfFromSyst(outFile, searchBinMap, systForConf, systHistoMap["search"]["highdm"]["znunu_zgammdiff"], "highdm", SB.high_dm_start)

    #-------------------------------------------------------
    # Calculate total systematic up/down
    #-------------------------------------------------------

    #TODO: combined all systematics to get total syst up/down
    # first combined syst in validation bins and give to Hui
    # loop over validation bins
    # loop over systematics
    # syst_up_i       = (p_up - p)   / p
    # syst_down_i     = (p - p_down) / p
    # syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
    # syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 
    
    # --- validation bins ---
    f_out = ROOT.TFile("validationBinsZinv_syst_" + era + ".root", "recreate")
    h_syst_up_lowdm     = ROOT.TH1F("syst_up_lowdm",    "syst_up_lowdm",    VB.low_dm_nbins,  VB.low_dm_start,  VB.low_dm_end  + 1)
    h_syst_up_highdm    = ROOT.TH1F("syst_up_highdm",   "syst_up_highdm",   VB.high_dm_nbins, VB.high_dm_start, VB.high_dm_end + 1)
    h_syst_down_lowdm   = ROOT.TH1F("syst_down_lowdm",  "syst_down_lowdm",  VB.low_dm_nbins,  VB.low_dm_start,  VB.low_dm_end  + 1)
    h_syst_down_highdm  = ROOT.TH1F("syst_down_highdm", "syst_down_highdm", VB.high_dm_nbins, VB.high_dm_start, VB.high_dm_end + 1)
    
    # TODO: remove old code; use grid instead
    # Draw nice ratio at 1
    #    # horizontal line at 1
    #    ratio_1 = histo[region]["down"].Clone()
    #    ratio_1.SetLineColor(ROOT.kBlack)
    #    for k in range(0,45):
    #        ratio_1.SetBinContent(k, 1)

    #    ratio_up.GetYaxis().SetRangeUser(0,1.5)
    #
    #    ratio_up.Draw("hist")
    #    ratio_down.Draw("hist same")
    #    ratio_1.Draw("hist same")
    validationBinMap = {}
    validationBinMap["lowdm"]   = VB.low_dm_bins
    validationBinMap["highdm"]  = VB.high_dm_bins
    # use copy.deepcopy() to avoid modifying original
    # histo_tmp[region][direction]
    validationHistoMap = copy.deepcopy(histo_tmp)
    validationHistoMap["lowdm"]["up"]       = h_syst_up_lowdm
    validationHistoMap["lowdm"]["down"]     = h_syst_down_lowdm
    validationHistoMap["highdm"]["up"]      = h_syst_up_highdm
    validationHistoMap["highdm"]["down"]    = h_syst_down_highdm
    
    h_pred_lowdm    = histo["validation"]["lowdm"][""]
    h_pred_highdm   = histo["validation"]["highdm"][""]

    # validation bins
    # bins are list of strings starting at 0
    # loop over regions (lowdm and highdm)
    print "# validation bin systematics"
    debug = False
    for region in regions:
        print "--- {0} ---".format(region)
        # DEBUG
        if debug:
            #systHistoMap[bintype][region][syst]
            for syst in systHistoMap["validation"][region]:
                nBins = systHistoMap["validation"][region][syst].GetNbinsX()
                error = systHistoMap["validation"][region][syst].GetBinContent(1)
                print "DEBUG: {0}, {1}: nBins = {2}, bin 1 error = {3}".format(region, syst, nBins, error)
        
        # be careful with bin index, which needs to start at 1 in both lowdm and highdm
        b_i = 1
        for b in validationBinMap[region]:
            p = histo["validation"][region][""].GetBinContent(b_i)
            syst_up_sum   = 0.0
            syst_down_sum = 0.0
            if p != 0:
                # syst from p, p_up, p_down
                for syst in systematics:
                    # do not apply SB syst in high dm
                    if region == "highdm" and syst == "eff_sb":
                        continue
                    # syst_histo[systemaitc][bintype][region][direction]
                    h_up    = syst_histo[syst]["validation"][region]["up"]
                    h_down  = syst_histo[syst]["validation"][region]["down"]
                    p_up    = h_up.GetBinContent(b_i)
                    p_down  = h_down.GetBinContent(b_i)
                    syst_up   = (p_up - p  ) / p
                    syst_down = (p - p_down) / p
                    syst_up_sum     += syst_up**2
                    syst_down_sum   += syst_down**2
                # syst from root file
                # symmetric error with up = down
                #systHistoMap[bintype][region][syst]
                for syst in systHistoMap["validation"][region]:
                    error = systHistoMap["validation"][region][syst].GetBinContent(b_i)
                    syst_up   = error  
                    syst_down = error  
                    syst_up_sum     += syst_up**2
                    syst_down_sum   += syst_down**2
            syst_up_total   = np.sqrt(syst_up_sum)
            syst_down_total = np.sqrt(syst_down_sum)
            final_up   = 1.0 + syst_up_total
            final_down = 1.0 - syst_down_total
            print "bin {0}, pred={1}, syst_up={2}, syst_down={3}".format(b_i, p, final_up, final_down)
            validationHistoMap[region]["up"].SetBinContent(     b_i, final_up   )
            validationHistoMap[region]["down"].SetBinContent(   b_i, final_down )
            b_i += 1
        
        # --- write histograms to file
        validationHistoMap[region]["up"].Write()
        validationHistoMap[region]["down"].Write()

        # --- plot histograms
                    
        # run_systematics.plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir)
        mySyst = "total"
        run_systematics.plot(histo["validation"][region][""], histo["validation"][region]["up"], histo["validation"][region]["down"], mySyst, "validation", region, era, syst_dir)

    f_out.Close()


if __name__ == "__main__":
    main()



