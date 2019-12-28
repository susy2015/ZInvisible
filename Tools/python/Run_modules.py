# Run_modules.py (modified for systematics)
import copy
import numpy as np
import json
import argparse
import os
import ROOT
import run_systematics
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from search_bins import  ValidationBins, SearchBins, CRUnitBins
from tools import invert, setupHist, stringifyMap, removeCuts

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# -------------------------------------------------------------
# Combine all systematics to get total syst up/down
# -------------------------------------------------------------
# loop over search or validation bins
# loop over systematics
# syst_up_i       = (p_up - p)   / p
# syst_down_i     = (p - p_down) / p
# syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
# syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 
# -------------------------------------------------------------

# Use histogram which stores systematic errors
# Includes functionality to write multiple syst. nuisances (e.g. Rz syst.) according to selection (e.g. Nb, Nsv)
def writeToConfFromSyst(outFile, binMap, process, syst, h, region, offset, selectionMap = {}, removeCut = ""):
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = binMap[str(sb_i)]
        # add selection to systematic name if selection map given
        final_syst = syst 
        if selectionMap:
            selection = selectionMap[str(sb_i)]["selection"]
            if removeCut:
                # remove cut only if removeCut is given
                # removeCuts(cutString, pattern, delim = "_"):
                selection = removeCuts(selection, removeCut)
            final_syst += "_" + selection
        
        # systematic error is stored in bin content
        # treat error as symmetric
        error   = h.GetBinContent(i)
        r_up   = 1 + error 
        r_down = 1 - error 

        outFile.write("{0}  {1}_Up    {2}  {3}\n".format( sb_name, final_syst, process, r_up   ) )
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, final_syst, process, r_down ) )

# use prediction, syst up/down histograms
def writeToConfFromPred(outFile, binMap, process, syst, h, h_up, h_down, region, offset):
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = binMap[str(sb_i)]
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        # If there is a zero prediction for nominal, up, and down can you set the systematic to 1.
        # The error is a deviation from 1.
        r_up   = 1
        r_down = 1
        # do not apply SB syst in high dm
        if not (region == "highdm" and syst == "ivfunc"):
            if p != 0:
                r_up   = p_up   / p
                r_down = p_down / p
            else:
                print "WARNING: pred = 0 for search bin {0}".format(sb_i)
        # syst is already systForConf (systForConf = systMap[syst]["name"])
        outFile.write("{0}  {1}_Up    {2}  {3}\n".format( sb_name, syst, process, r_up   ) )
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, syst, process, r_down ) )

# symmetrize systematic if up/down variation is in the same direction compared to nominal
# modify histograms passed to function
def symmetrizeSyst(h, h_up, h_down):
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        # symmetrize systematic if up/down variation is in the same direction compared to nominal
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        diff_up   = p_up - p
        diff_down = p_down - p
        same_dir   = diff_up * diff_down > 0
        if same_dir:
            diff_symm = np.mean([abs(diff_up), abs(diff_down)])
            h_up.SetBinContent(   i, p + diff_symm)
            h_down.SetBinContent( i, p - diff_symm)


def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--runs_json",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--syst_json",    "-s", default="",                             help="json file containing systematics")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options                     = parser.parse_args()
    units_json                  = "dc_BkgPred_BinMaps_master.json"
    searchbin_selection_json    = "search_bins_v4.json"
    runs_json                   = options.runs_json
    syst_json                   = options.syst_json
    verbose                     = options.verbose
    # Note: DO NOT apply Rz syst. in CR unit bins
    rz_syst_files   = {}
    rz_syst_files["validation"]   = "RzSyst_ValidationBins.root"
    rz_syst_files["search"]       = "RzSyst_SearchBins.root"
    # Note: apply Z vs. Photon syst. in CR unit bins
    ZvPhoton_syst_files = {}
    ZvPhoton_syst_files["validation"]   = "ZvsPhotonSyst_ValidationBins.root"
    ZvPhoton_syst_files["search"]       = "ZvsPhotonSyst_SearchBins.root" 
    ZvPhoton_syst_files["controlUnit"]  = "ZvsPhotonSyst_CRUnitBins.root"
   
    doSymmetrize            = True
    doUnits                 = True
    draw                    = False
    saveRootFile            = False
    masterBinMap            = {}
    searchBinSelectionMap   = {}
    runMap                  = {}
    systMap                 = {}
    systHistoMap            = {}

    # list required files to check if they exist
    file_list = [] 
    file_list.append(units_json)
    file_list.append(searchbin_selection_json)
    file_list.append(runs_json)
    file_list.append(syst_json)
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
    with open(units_json, "r") as input_file:
        masterBinMap = json.load(input_file) 
    with open(searchbin_selection_json, "r") as input_file:
        searchBinSelectionMap = stringifyMap(json.load(input_file))
    with open(runs_json, "r") as input_file:
        runMap  = json.load(input_file)
    with open(syst_json, "r") as input_file:
        systMap = json.load(input_file) 
    
    # map search bin numbers to string names
    searchBinMap = invert(masterBinMap["binNum"])
    unitBinMap   = invert(masterBinMap["unitCRNum"]["phocr"])

    eras            = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM", "Run2"]
    era             = "Run2"
    runDir          = runMap[era]
    result_file     = "condor/" + runDir + "/result.root"
    conf_file       = "datacard_inputs/zinv_syst_" + era + ".conf"
    total_syst_dir  = "prediction_histos/"
    out_dir         = "syst_plots/" 
    tmp_dir         = "tmp_plots/"
    variable        = "mc"
    regions         = ["lowdm", "highdm"]
    directions      = ["up", "", "down"]
    bintypes        = ["validation", "search", "controlUnit"]
    # including prefire; WARNING: prefire only exists in (2016,2017) and needs to be handled carefully 
    systematics_znunu  = ["jes","btag","eff_restoptag","eff_sb","eff_toptag","eff_wtag","met_trig","pileup","prefire"]
    systematics_phocr  = ["jes","btag","eff_restoptag_photon","eff_sb_photon","eff_toptag_photon","eff_wtag_photon","photon_trig","pileup","prefire","photon_sf"]
    systematics = list(set(systematics_znunu)| set(systematics_phocr))
    # TODO. missing systematics: pdf_weight, MET resolution (uncluster), lepton veto SF, ISR_Weight

    # histo_tmp[region][direction]
    histo_tmp  = {region:dict.fromkeys(directions) for region in regions}

    # histo[bintype][region][direction]
    histo      = {bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} 

    # syst_histo[systemaitc][bintype][region][direction]
    syst_histo = {syst:{bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} for syst in systematics} 

    #-------------------------------------------------------
    # Class instanceses summoning 
    #-------------------------------------------------------

    N   =  Normalization(tmp_dir, verbose)
    S   =  Shape(tmp_dir, draw, doUnits, verbose)
    VB  =  ValidationBins(  N, S, eras, tmp_dir, verbose, draw, saveRootFile)
    SB  =  SearchBins(      N, S, eras, tmp_dir, verbose, draw, saveRootFile)
    CRU =  CRUnitBins(      N, S, eras, tmp_dir, verbose)
    
    #-------------------------------------------------------
    # Normal predictions (no systematics) 
    #-------------------------------------------------------
    
    # total syst. calculated for Run 2
    # however, loop over all eras to fix syst. not present in all eras (e.g. prefire)
    # warning; be careful to not use variables "era" and "result_file" as those are only for total era
    for e in eras:
        d = runMap[e]
        r  = "condor/" + d + "/result.root"
        N.getNormAndError(r, e)
        S.getShape(r, e)
        VB.getValues(r, e)
        SB.getValues(r, e)
        CRU.getValues(r, e)
    
    for region in regions:
        histo["validation"][region][""]     =  VB.histograms[era][region][variable].Clone()
        histo["search"][region][""]         =  SB.histograms[era][region][variable].Clone()
        histo["controlUnit"][region][""]    =  CRU.histograms[era][region]["mc_gjets"].Clone()
    
    #-------------------------------------------------------
    # Calculate normalization and shape factors
    #-------------------------------------------------------
    with open(conf_file, "w") as outFile:
        
        #-------------------------------------------------------
        # Systematics for znunu
        #-------------------------------------------------------
           
        # --- begin loop over systematics
        for syst in systematics_znunu:
            for direction in ["up", "down"]:
                systTag = "_{0}_syst_{1}".format(syst, direction)
        
                # --- calculate variation for this systematic
                # only vary Z nu nu; do not vary Normalization or Shape
                VB.getValues(       result_file,  era, systTag)
                SB.getValues(       result_file,  era, systTag)
                
                for region in regions:
                    histo["validation"][region][direction]    =  VB.histograms[era][region][variable].Clone()
                    histo["search"][region][direction]        =  SB.histograms[era][region][variable].Clone()
                    
                    # fix prefire because it does not have a weight or syst. in 2018, but 2018 hist. was not added
                    if syst == "prefire":
                        for e in ["2018_PreHEM", "2018_PostHEM"]:
                            histo["validation"][region][direction].Add(     VB.histograms[e][region][variable]         )
                            histo["search"][region][direction].Add(         SB.histograms[e][region][variable]         )
                    
                    syst_histo[syst]["validation"][region][direction]    = histo["validation"][region][direction]
                    syst_histo[syst]["search"][region][direction]        = histo["search"][region][direction]
            
            #-------------------------------------------------------
            # Symmetrize systematics which are in the same direction
            #-------------------------------------------------------

            if doSymmetrize:
                for region in regions:
                    # symmetrize systematic if up/down variation is in the same direction compared to nominal
                    # modify histograms passed to function
                    # symmetrizeSyst(h, h_up, h_down)
                    symmetrizeSyst(histo["validation"][region][""],     syst_histo[syst]["validation"][region]["up"],   syst_histo[syst]["validation"][region]["down"])
                    symmetrizeSyst(histo["search"][region][""],         syst_histo[syst]["search"][region]["up"],       syst_histo[syst]["search"][region]["down"])
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------

            # WARNING: offset is starting point for low/high dm bins, use with care
            #writeToConfFromPred(outFile, binMap, process, syst, h, h_up, h_down, region, offset):
            systForConf = systMap[syst]["name"]  
            writeToConfFromPred(outFile, searchBinMap, "znunu", systForConf, histo["search"]["lowdm"][""],  histo["search"]["lowdm"]["up"],  histo["search"]["lowdm"]["down"],  "lowdm",  SB.low_dm_start)
            writeToConfFromPred(outFile, searchBinMap, "znunu", systForConf, histo["search"]["highdm"][""], histo["search"]["highdm"]["up"], histo["search"]["highdm"]["down"], "highdm", SB.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------
            
            for bintype in ["validation", "search"]:
                for region in regions:
                    # run_systematics.plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir)
                    run_systematics.plot(histo[bintype][region][""], histo[bintype][region]["up"], histo[bintype][region]["down"], syst, bintype, region, era, out_dir)
                
        # --- end loop over systematics
        
        #-------------------------------------------------------
        # Systematics for phocr
        #-------------------------------------------------------
        
        # --- begin loop over systematics
        for syst in systematics_phocr:
            for direction in ["up", "down"]:
                systTag = "_{0}_syst_{1}".format(syst, direction)
        
                # --- calculate variation for this systematic
                # Do not vary Normalization
                # Vary shape to get GJets MC variation
                S.getShape(         result_file,  era, systTag)
                CRU.getValues(      result_file,  era)
                
                for region in regions:
                    histo["controlUnit"][region][direction]   =  CRU.histograms[era][region]["mc_gjets"].Clone()
                    
                    # fix prefire because it does not have a weight or syst. in 2018, but 2018 hist. was not added
                    if syst == "prefire":
                        for e in ["2018_PreHEM", "2018_PostHEM"]:
                            histo["controlUnit"][region][direction].Add(    CRU.histograms[e][region]["mc_gjets"]      )
                    
                    syst_histo[syst]["controlUnit"][region][direction]   = histo["controlUnit"][region][direction]
            
            #-------------------------------------------------------
            # Symmetrize systematics which are in the same direction
            #-------------------------------------------------------

            if doSymmetrize:
                for region in regions:
                    # symmetrize systematic if up/down variation is in the same direction compared to nominal
                    # modify histograms passed to function
                    # symmetrizeSyst(h, h_up, h_down)
                    symmetrizeSyst(histo["controlUnit"][region][""],    syst_histo[syst]["controlUnit"][region]["up"],  syst_histo[syst]["controlUnit"][region]["down"])
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------

            # WARNING: offset is starting point for low/high dm bins, use with care
            #writeToConfFromPred(outFile, binMap, process, syst, h, h_up, h_down, region, offset):
            systForConf = systMap[syst]["name"]  
            writeToConfFromPred(outFile, unitBinMap,   "phocr_gjets", systForConf, histo["controlUnit"]["lowdm"][""],  histo["controlUnit"]["lowdm"]["up"],  histo["controlUnit"]["lowdm"]["down"],  "lowdm",  CRU.low_dm_start)
            writeToConfFromPred(outFile, unitBinMap,   "phocr_gjets", systForConf, histo["controlUnit"]["highdm"][""], histo["controlUnit"]["highdm"]["up"], histo["controlUnit"]["highdm"]["down"], "highdm", CRU.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------
            
            for bintype in ["controlUnit"]:
                for region in regions:
                    # run_systematics.plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir)
                    run_systematics.plot(histo[bintype][region][""], histo[bintype][region]["up"], histo[bintype][region]["down"], syst, bintype, region, era, out_dir)
                
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
            # Note: DO NOT apply Rz syst. in CR unit bins
            if bintype != "controlUnit":
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
        # writeToConfFromSyst(outFile, binMap, process, syst, h, region, offset, selectionMap = {}, removeCut = "")
        systForConf = systMap["znunu_rzunc"]["name"]  
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["lowdm"]["znunu_rzunc"],  "lowdm",  SB.low_dm_start,  searchBinSelectionMap, "NJ")
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["highdm"]["znunu_rzunc"], "highdm", SB.high_dm_start, searchBinSelectionMap, "NJ")
        
        # --- Z vs Photon syst --- #
        # writeToConfFromSyst(outFile, binMap, process, syst, h, region, offset, selectionMap = {}, removeCut = "")
        systForConf = systMap["znunu_zgammdiff"]["name"]  
        writeToConfFromSyst(outFile, unitBinMap, "phocr_gjets", systForConf, systHistoMap["controlUnit"]["lowdm"]["znunu_zgammdiff"],  "lowdm",  CRU.low_dm_start)
        writeToConfFromSyst(outFile, unitBinMap, "phocr_gjets", systForConf, systHistoMap["controlUnit"]["highdm"]["znunu_zgammdiff"], "highdm", CRU.high_dm_start)

    #-------------------------------------------------------
    # Calculate total systematic up/down
    #-------------------------------------------------------

    # loop over validation bins
    # loop over systematics
    # syst_up_i       = (p_up - p)   / p
    # syst_down_i     = (p - p_down) / p
    # syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
    # syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 
    
    # --- validation bins ---
    f_out = ROOT.TFile(total_syst_dir + "validationBinsZinv_syst_" + era + ".root", "recreate")
    h_syst_up_lowdm     = ROOT.TH1F("syst_up_lowdm",    "syst_up_lowdm",    VB.low_dm_nbins,  VB.low_dm_start,  VB.low_dm_end  + 1)
    h_syst_up_highdm    = ROOT.TH1F("syst_up_highdm",   "syst_up_highdm",   VB.high_dm_nbins, VB.high_dm_start, VB.high_dm_end + 1)
    h_syst_down_lowdm   = ROOT.TH1F("syst_down_lowdm",  "syst_down_lowdm",  VB.low_dm_nbins,  VB.low_dm_start,  VB.low_dm_end  + 1)
    h_syst_down_highdm  = ROOT.TH1F("syst_down_highdm", "syst_down_highdm", VB.high_dm_nbins, VB.high_dm_start, VB.high_dm_end + 1)
    
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
        # get histograms for this region
        h_total_syst_up   = validationHistoMap[region]["up"]
        h_total_syst_down = validationHistoMap[region]["down"]

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
                for syst in systematics_znunu:
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
                #systHistoMap[bintype][region][syst]
                for syst in systHistoMap["validation"][region]:
                    error = systHistoMap["validation"][region][syst].GetBinContent(b_i)
                    # symmetric error with up = down
                    syst_up   = error  
                    syst_down = error  
                    syst_up_sum     += syst_up**2
                    syst_down_sum   += syst_down**2
            syst_up_total   = np.sqrt(syst_up_sum)
            syst_down_total = np.sqrt(syst_down_sum)
            final_up   = 1.0 + syst_up_total
            final_down = 1.0 - syst_down_total
            #print "bin {0}, pred={1}, syst_up={2}, syst_down={3}".format(b_i, p, final_up, final_down)
            h_total_syst_up.SetBinContent(     b_i, final_up   )
            h_total_syst_down.SetBinContent(   b_i, final_down )
            b_i += 1
        
        # --- write histograms to file
        h_total_syst_up.Write()
        h_total_syst_down.Write()

        #-------------------------------------------------------
        # Plot total systematic up/down
        #-------------------------------------------------------
                    
        # correct plot
        bintype = "validation"
        mySyst = "total"
        
        eraTag = "_" + era
        draw_option = "hist"
        # colors
        color_red    = "vermillion"
        color_blue   = "electric blue"
        color_green  = "irish green" 
        color_purple = "violet"
        color_black  = "black"
        
        # legend: TLegend(x1,y1,x2,y2)
        legend_x1 = 0.6
        legend_x2 = 0.9 
        legend_y1 = 0.7
        legend_y2 = 0.9 
    
        c = ROOT.TCanvas("c", "c", 800, 800)
        name = "{0}_{1}_syst".format(bintype, mySyst)
        
        title = "Z to Invisible: " + name + " in " + region + " for " + era
        x_title = "Validation Bins"
        setupHist(h_total_syst_up,     title, x_title, "total systematic",  color_red,    0.0, 2.0)
        setupHist(h_total_syst_down,   title, x_title, "total systematic",  color_blue,   0.0, 2.0)
        
        # draw histograms
        h_total_syst_up.Draw(draw_option)
        h_total_syst_down.Draw(draw_option + " same")
        
        # legend: TLegend(x1,y1,x2,y2)
        legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
        legend.AddEntry(h_total_syst_up,           "syst up",                  "l")
        legend.AddEntry(h_total_syst_down,         "syst down",                "l")
        legend.Draw()
        
        
        # save histograms
        plot_name = out_dir + name + "_" + region + eraTag
        c.Update()
        c.SaveAs(plot_name + ".png")
        del c

    f_out.Close()


if __name__ == "__main__":
    main()



