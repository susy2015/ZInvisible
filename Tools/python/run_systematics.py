# run_systematics.py (modified for systematics)
import warnings
import copy
import numpy as np
import json
import argparse
import os
import ROOT
from norm_lepton_zmass import Normalization
from shape_photon_met import Shape
from search_bins import  ValidationBins, ValidationBinsMETStudy, SearchBins, CRUnitBins
from tools import invert, setupHist, stringifyMap, removeCuts, ERROR_SYST, isclose

# set numpy to provide warnings instead of just printing errors
np.seterr(all='warn')
warnings.filterwarnings('error')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# ----------------
# Plot
# ----------------

def plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir, mc_label):

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
    c.Divide(1, 2)
    
    name = "{0}_{1}_syst".format(bintype, mySyst)
    
    h_ratio_up      = h_up.Clone("h_ratio_up") 
    h_ratio_down    = h_down.Clone("h_ratio_down") 
    h_ratio_up.Divide(h)
    h_ratio_down.Divide(h)

    fixSystForZeroPred(h, h_ratio_up, h_ratio_down)
    
    title = "Z to Invisible: " + name + " in " + region + " for " + era
    x_title = bintype + " bins"
    setupHist(h,                title, x_title, "Events",               color_black,  10.0 ** -2, 10.0 ** 5)
    setupHist(h_up,             title, x_title, "Events",               color_red,    10.0 ** -2, 10.0 ** 5)
    setupHist(h_down,           title, x_title, "Events",               color_blue,   10.0 ** -2, 10.0 ** 5)
    setupHist(h_ratio_up,       title, x_title, "variation / nominal",  color_red,    0.5, 1.5)
    setupHist(h_ratio_down,     title, x_title, "variation / nominal",  color_blue,   0.5, 1.5)
    
    # histograms
    c.cd(1)
    ROOT.gPad.SetLogy(1) # set log y
    # draw
    h.Draw(draw_option)
    h_up.Draw(draw_option + " same")
    h_down.Draw(draw_option + " same")
    
    # legend: TLegend(x1,y1,x2,y2)
    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
    legend.AddEntry(h,              mc_label,                 "l")
    legend.AddEntry(h_up,           "syst up",                "l")
    legend.AddEntry(h_down,         "syst down",              "l")
    legend.Draw()
    
    # ratios
    pad = c.cd(2)
    pad.SetGrid()
    # draw
    h_ratio_up.Draw(draw_option)
    h_ratio_down.Draw(draw_option + " same")
    
    # save histograms
    plot_name = plot_dir + name + "_" + region + eraTag
    c.Update()
    c.SaveAs(plot_name + ".png")






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
def writeToConfFromSyst(outFile, binMap, process, syst, h, region, offset, selectionMap = {}, removeCut = []):
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        sb_i = i - 1 + offset
        sb_name = binMap[str(sb_i)]
        # add selection to systematic name if selection map given
        final_syst = syst + "_" + region
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
# syst name is "syst for conf" name from systematics.json
def writeToConfFromPred(outFile, infoFile, binMap, process, syst, h, h_up, h_down, region, offset):
    nBins = h.GetNbinsX()
    values_up   = []
    values_down = []
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
                print "WARNING: pred = 0 for bin {0}".format(sb_i)
        
        # syst is already systForConf (systForConf = systMap[syst]["name"])
        # investigate jes
        #if process == "znunu" and syst == "JES":
        #    print "SEARCH_BIN_{0}  {1}  {2}  {3}: nominal={4}, up={5}, down={6}".format(sb_i, sb_name, syst, process, p, p_up, p_down)
        
        # avoid taking log of negative number or 0
        new_up   = r_up
        new_down = r_down
        if new_up <= 0:
            new_up = abs(ERROR_SYST * p)
        if new_down <= 0:
            new_down = abs(ERROR_SYST * p)
        # make both up and down variations >= 1.0 independent of direction
        values_up.append(np.exp(abs(np.log(new_up))))
        values_down.append(np.exp(abs(np.log(new_down))))
        #values_up.append(r_up)
        #values_down.append(r_down)
        
        outFile.write("{0}  {1}_Up    {2}  {3}\n".format( sb_name, syst, process, r_up   ) )
        outFile.write("{0}  {1}_Down  {2}  {3}\n".format( sb_name, syst, process, r_down ) )
    
    # get percentages >= 0.0%
    ave_up   = 100 * (np.average(values_up) - 1)
    ave_down = 100 * (np.average(values_down) - 1)
    min_up   = 100 * (min(values_up) - 1)
    min_down = 100 * (min(values_down) - 1)
    max_up   = 100 * (max(values_up) - 1)
    max_down = 100 * (max(values_down) - 1)
    # record average systematic here
    infoFile.write("{0}  {1}  {2}_Up    ave={3:.3f}%, range=[{4:.3f}%, {5:.3f}%]\n".format(process, region, syst, ave_up,   min_up,   max_up   ))
    infoFile.write("{0}  {1}  {2}_Down  ave={3:.3f}%, range=[{4:.3f}%, {5:.3f}%]\n".format(process, region, syst, ave_down, min_down, max_down ))

# set p_up and p_down to 1.0 if p is 0.0
# modify histograms passed to function
def fixSystForZeroPred(h, h_up, h_down):
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        # set p_up and p_down to 1.0 if p is 0.0
        if p == 0:
            h_up.SetBinContent(   i, 1.0)
            h_down.SetBinContent( i, 1.0)

# symmetrize systematic if up/down variation is in the same direction compared to nominal
# modify histograms passed to function
def symmetrizeSyst(h, h_up, h_down):
    useLogNormal = True
    nBins = h.GetNbinsX()
    for i in xrange(1, nBins + 1):
        p       = h.GetBinContent(i)
        p_up    = h_up.GetBinContent(i)
        p_down  = h_down.GetBinContent(i)
        # symmetrize systematic if up/down variation is in the same direction compared to nominal
        diff_up   = p_up - p
        diff_down = p_down - p
        same_dir   = diff_up * diff_down > 0
        
        # default for p = 0
        Aup   = 1
        Adown = 1
        if p:
            Aup   = p_up / p
            Adown = p_down / p
        
        if same_dir:
            if useLogNormal:
                # using geometric mean
                if Aup * Adown < 0:
                    print "WARNING: Aup * Adown = {0}, about to take square root.".format(Aup * Adown)
                #geometric_mean = np.sqrt(Aup * Adown)
                # HACK to avoid negative values
                geometric_mean = np.sqrt(abs(Aup * Adown))
                Aup     /= geometric_mean
                Adown   /= geometric_mean
                h_up.SetBinContent(   i, p * Aup)
                h_down.SetBinContent( i, p * Adown)
            else:
                # using arithmetic mean
                ave_diff = np.mean([abs(diff_up), abs(diff_down)])
                if p_up >= p_down:
                    h_up.SetBinContent(   i, p + ave_diff)
                    h_down.SetBinContent( i, p - ave_diff)
                else:
                    h_up.SetBinContent(   i, p - ave_diff)
                    h_down.SetBinContent( i, p + ave_diff)

# Total systematics function which can run on validaiton, MET study, search bins, CR unit bins, etc.

def getTotalSystematics(BinObject, bintype, systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2):
    
    #-------------------------------------------------------
    # Calculate total systematic up/down
    #-------------------------------------------------------

    # loop over bins
    # loop over systematics
    # syst_up_i       = (p_up - p)   / p
    # syst_down_i     = (p - p_down) / p
    # syst_up_total   = sqrt ( sum ( syst_up_i ^2 ) ) 
    # syst_down_total = sqrt ( sum ( syst_down_i ^2 ) ) 
    
    verbose = False
    useLogNormal    = True
    total_syst_dir  = "prediction_histos/"
    # histo_tmp[region][direction]
    histo_tmp  = {region:dict.fromkeys(directions) for region in regions}

    # --- bins --- #
    f_out = ROOT.TFile(total_syst_dir + bintype + "BinsZinv_syst_" + era + ".root", "recreate")
    h_syst_up_lowdm     = ROOT.TH1F("syst_up_lowdm",    "syst_up_lowdm",    BinObject.low_dm_nbins,  BinObject.low_dm_start,  BinObject.low_dm_end  + 1)
    h_syst_up_highdm    = ROOT.TH1F("syst_up_highdm",   "syst_up_highdm",   BinObject.high_dm_nbins, BinObject.high_dm_start, BinObject.high_dm_end + 1)
    h_syst_down_lowdm   = ROOT.TH1F("syst_down_lowdm",  "syst_down_lowdm",  BinObject.low_dm_nbins,  BinObject.low_dm_start,  BinObject.low_dm_end  + 1)
    h_syst_down_highdm  = ROOT.TH1F("syst_down_highdm", "syst_down_highdm", BinObject.high_dm_nbins, BinObject.high_dm_start, BinObject.high_dm_end + 1)
    
    myBinMap = {}
    myBinMap["lowdm"]   = BinObject.low_dm_bins
    myBinMap["highdm"]  = BinObject.high_dm_bins
    # use copy.deepcopy() to avoid modifying original
    # histo_tmp[region][direction]
    myHistoMap = copy.deepcopy(histo_tmp)
    myHistoMap["lowdm"]["up"]       = h_syst_up_lowdm
    myHistoMap["lowdm"]["down"]     = h_syst_down_lowdm
    myHistoMap["highdm"]["up"]      = h_syst_up_highdm
    myHistoMap["highdm"]["down"]    = h_syst_down_highdm
    
    h_pred_lowdm    = histo[bintype]["lowdm"][""]
    h_pred_highdm   = histo[bintype]["highdm"][""]

    # bins are list of strings starting at 0
    # loop over regions (lowdm and highdm)
    if verbose:
        print "# --- {0} bin systematics --- #".format(bintype)
    debug = False
    for region in regions:
        # get histograms for this region
        h_total_syst_up   = myHistoMap[region]["up"]
        h_total_syst_down = myHistoMap[region]["down"]

        # DEBUG
        if debug:
            #systHistoMap[bintype][region][syst]
            for syst in systHistoMap[bintype][region]:
                nBins = systHistoMap[bintype][region][syst].GetNbinsX()
                error = systHistoMap[bintype][region][syst].GetBinContent(1)
                print "DEBUG: {0}, {1}: nBins = {2}, bin 1 error = {3}".format(region, syst, nBins, error)
        
        # be careful with bin index, which needs to start at 1 in both lowdm and highdm
        b_i = 1
        for b in myBinMap[region]:
            p = histo[bintype][region][""].GetBinContent(b_i)
            syst_up_sum        = 0.0
            syst_down_sum      = 0.0
            log_syst_up_sum    = 0.0
            log_syst_down_sum  = 0.0
            if p != 0:
                # syst from p, p_up, p_down
                for syst in systematics_znunu:
                    # do not apply SB syst in high dm
                    if region == "highdm" and syst == "eff_sb":
                        continue
                    # syst_histo[systemaitc][bintype][region][direction]
                    h_up    = syst_histo[syst][bintype][region]["up"]
                    h_down  = syst_histo[syst][bintype][region]["down"]
                    p_up    = h_up.GetBinContent(b_i)
                    p_down  = h_down.GetBinContent(b_i)
                    syst_up         = (p_up - p  ) / p
                    syst_down       = (p - p_down) / p
                    log_syst_up     = p_up / p
                    log_syst_down   = p_down / p
                    # avoid taking log of negative number or 0
                    if log_syst_up <= 0:
                        new_value = ERROR_SYST * p
                        print "WARNING: For {0} bin {1}, syst {2}: log_syst_up = {3}. Setting log_syst_up = {4}.".format(bintype, b, syst, log_syst_up, new_value)
                        log_syst_up = new_value
                    if log_syst_down <= 0:
                        new_value = ERROR_SYST * p
                        print "WARNING: For {0} bin {1}, syst {2}: log_syst_down = {3}. Setting log_syst_down = {4}.".format(bintype, b, syst, log_syst_down, new_value)
                        log_syst_down = new_value
                    # sum in quadrature 
                    syst_up_sum     += syst_up**2
                    syst_down_sum   += syst_down**2
                    # If both systematics go the same direction, need to symmetrize
                    # Because all the nuisance parameters are log-normal, symmetrize by dividing by the geometric mean
                    if ((log_syst_up > 1) and (log_syst_down > 1)) or ((log_syst_up < 1) and (log_syst_down < 1)):
                        geometric_mean = np.sqrt(log_syst_up * log_syst_down)
                        log_syst_up   /= geometric_mean
                        log_syst_down /= geometric_mean
                    # Because all the nuisance parameters are log-normal, sum the log of the ratios in quadrature
                    # Sum (the square of the log of) all the ratios that are greater than 1
                    # Sum (the square of the log of) all the ratios that are less than 1
                    # Then at the end, take the exponential of the square root of each sum to get the total systematic ratio
                    try: 
                        if log_syst_up > 1 or log_syst_down < 1:
                            log_syst_up_sum     += np.log(log_syst_up)**2
                            log_syst_down_sum   += np.log(log_syst_down)**2
                        else:
                            log_syst_up_sum     += np.log(log_syst_down)**2
                            log_syst_down_sum   += np.log(log_syst_up)**2
                    except:
                        print "ERROR for np.log(), location 1: syst = {0}, log_syst_up = {1}, log_syst_down = {2}".format(syst, log_syst_up, log_syst_down)
                    # check leading systematic for bin 61
                    # if b == "61":
                    #     print "BIN_61: {0} up={1}, down={2}".format(syst, log_syst_up, log_syst_down)
                # these syst. are only available when doing Run 2
                if doRun2:
                    # syst from root file
                    #systHistoMap[bintype][region][syst]
                    for syst in systHistoMap[bintype][region]:
                        error = systHistoMap[bintype][region][syst].GetBinContent(b_i)
                        # symmetric error with up = down
                        syst_up         = error  
                        syst_down       = error  
                        log_syst_up     = 1.0 + error 
                        log_syst_down   = 1.0 - error
                        # avoid taking log of negative number or 0
                        if log_syst_up <= 0:
                            new_value = ERROR_SYST * p
                            print "WARNING: For {0} bin {1}, syst {2}: log_syst_up = {3}. Setting log_syst_up = {4}.".format(bintype, b, syst, log_syst_up, new_value)
                            log_syst_up = new_value
                        if log_syst_down <= 0:
                            new_value = ERROR_SYST * p
                            print "WARNING: For {0} bin {1}, syst {2}: log_syst_down = {3}. Setting log_syst_down = {4}.".format(bintype, b, syst, log_syst_down, new_value)
                            log_syst_down = new_value
                        # sum in quadrature 
                        syst_up_sum     += syst_up**2
                        syst_down_sum   += syst_down**2
                        # Because all the nuisance parameters are log-normal, sum the log of the ratios in quadrature
                        try: 
                            if log_syst_up > 1 or log_syst_down < 1:
                                log_syst_up_sum     += np.log(log_syst_up)**2
                                log_syst_down_sum   += np.log(log_syst_down)**2
                            else:
                                log_syst_up_sum     += np.log(log_syst_down)**2
                                log_syst_down_sum   += np.log(log_syst_up)**2
                        except:
                            print "ERROR for np.log(), location 2: syst = {0}, log_syst_up = {1}, log_syst_down = {2}".format(syst, log_syst_up, log_syst_down)
            syst_up_total   = np.sqrt(syst_up_sum)
            syst_down_total = np.sqrt(syst_down_sum)
            final_up   = 1.0 + syst_up_total
            final_down = 1.0 - syst_down_total
            log_syst_up_total   = np.exp( np.sqrt(log_syst_up_sum))
            log_syst_down_total = np.exp(-np.sqrt(log_syst_down_sum)) # Minus sign is needed because this is the *down* ratio
            log_final_up   = log_syst_up_total
            log_final_down = log_syst_down_total
            if verbose:
                print "bin {0}, pred={1}, syst_up={2}, syst_down={3}, log_final_up={4}, log_final_down={5}".format(b_i, p, final_up, final_down, log_final_up, log_final_down)
            if useLogNormal:
                h_total_syst_up.SetBinContent(     b_i, log_final_up   )
                h_total_syst_down.SetBinContent(   b_i, log_final_down )
            else:
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
        x_title = bintype + " bins"
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



#-------------------------------------------------------
# Calculate total systematic up/down
# - variation in prediction
# - include both Z nunu search bin and photon CR bin variation
#-------------------------------------------------------

def getTotalSystematicsPrediction(SearchBinObject, CRBinObject, N, S, runMap, systematics_map, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2):
    d = runMap[era]
    r  = "condor/" + d + "/result.root"
    N.getNormAndError(r, era)
    S.getShape(r, era)
    CRBinObject.getValues(r, era)
    SearchBinObject.getValues(r, era, CRunits=CRBinObject)
    
    with open("units.json", "r") as input_file:
        searchToUnitMap = json.load(input_file)
    
    total_syst_dir  = "prediction_histos/"
    # histo_tmp[region][direction]
    histo_tmp  = {region:dict.fromkeys(directions) for region in regions}
    
    # --- bins --- #
    f_out = ROOT.TFile(total_syst_dir + "searchBinsZinv_totalPredSyst_" + era + ".root", "recreate")
    h_syst_up_lowdm     = ROOT.TH1F("syst_up_lowdm",    "syst_up_lowdm",    SearchBinObject.low_dm_nbins,  SearchBinObject.low_dm_start,  SearchBinObject.low_dm_end  + 1)
    h_syst_up_highdm    = ROOT.TH1F("syst_up_highdm",   "syst_up_highdm",   SearchBinObject.high_dm_nbins, SearchBinObject.high_dm_start, SearchBinObject.high_dm_end + 1)
    h_syst_down_lowdm   = ROOT.TH1F("syst_down_lowdm",  "syst_down_lowdm",  SearchBinObject.low_dm_nbins,  SearchBinObject.low_dm_start,  SearchBinObject.low_dm_end  + 1)
    h_syst_down_highdm  = ROOT.TH1F("syst_down_highdm", "syst_down_highdm", SearchBinObject.high_dm_nbins, SearchBinObject.high_dm_start, SearchBinObject.high_dm_end + 1)
    
    myBinMap = {}
    myBinMap["lowdm"]   = SearchBinObject.low_dm_bins
    myBinMap["highdm"]  = SearchBinObject.high_dm_bins
    # use copy.deepcopy() to avoid modifying original
    # histo_tmp[region][direction]
    myHistoMap = copy.deepcopy(histo_tmp)
    myHistoMap["lowdm"]["up"]       = h_syst_up_lowdm
    myHistoMap["lowdm"]["down"]     = h_syst_down_lowdm
    myHistoMap["highdm"]["up"]      = h_syst_up_highdm
    myHistoMap["highdm"]["down"]    = h_syst_down_highdm
    
    for region in regions:
        # offset = 0
        # if region == "highdm":
        #     offset = 53
        # get histograms for this region
        h_total_syst_up   = myHistoMap[region]["up"]
        h_total_syst_down = myHistoMap[region]["down"]
        # be careful with bin index, which needs to start at 1 in both lowdm and highdm
        b_i = 1
        for b in myBinMap[region]:
            norm  = SearchBinObject.binValues[era][b]["norm"]
            shape = SearchBinObject.binValues[era][b]["shape"]
            mc    = histo["search"][region][""].GetBinContent(b_i) 
            
            p1 = shape * norm * mc 
            p2 = SearchBinObject.binValues[era][b]["pred"]
            
            p  = p1

            print "PREDICTION search bin {0}: {1}, {2}, {3}".format(b, p1, p2, isclose(p1, p2, rel_tol=1e-06))
            
            # norm        = SearchBinObject.binValues[era][b]["norm"]
            # znunu       = histo["search"][region][""].GetBinContent(b_i)
            # h_pho_data  = histo["controlUnit_data"][region][""]
            # h_pho_gjets = histo["controlUnit_gjets"][region][""]
            # h_pho_back  = histo["controlUnit_back"][region][""]
    
            # # get CR unit bins for this search bin
            # cr_units = searchToUnitMap["unitBinMapCR_phocr"][b]
            # # get histogram bin numbers
            # cr_units_bin_i = [int(cr) + 1 - offset for cr in cr_units]
            # 
            # # check that prediction calculation is correct
            # data_list   = [h_pho_data[i]           for i in cr_units_bin_i]
            # gjets_list  = [h_pho_gjets[i]          for i in cr_units_bin_i]
            # back_list   = [h_pho_back[i]           for i in cr_units_bin_i]
            
            log_syst_up_sum    = 0.0
            log_syst_down_sum  = 0.0
            
            if p != 0:
                # syst from p, p_up, p_down
                for key in systematics_map:
                    syst_znunu = systematics_map[key]["znunu"]
                    syst_phocr = systematics_map[key]["phocr"]
                    # do not apply SB syst in high dm
                    if region == "highdm" and key == "eff_sb":
                        continue
                    # syst_histo[systemaitc][bintype][region][direction]
                    # some systematics are only for Znunu or only for photon CR
                    if syst_znunu:
                        h_znunu_up    = syst_histo[syst_znunu]["search"][region]["up"]
                        h_znunu_down  = syst_histo[syst_znunu]["search"][region]["down"]
                    else:
                        h_znunu_up   = histo["search"][region][""]
                        h_znunu_down = histo["search"][region][""]
                    if syst_phocr:
                        h_gjets_up    = syst_histo[syst_phocr]["controlUnit_gjets"][region]["up"]
                        h_gjets_down  = syst_histo[syst_phocr]["controlUnit_gjets"][region]["down"]
                        h_back_up     = syst_histo[syst_phocr]["controlUnit_back"][region]["up"]
                        h_back_down   = syst_histo[syst_phocr]["controlUnit_back"][region]["down"]
                    else:
                        h_gjets_up    = histo["controlUnit_gjets"][region][""] 
                        h_gjets_down  = histo["controlUnit_gjets"][region][""] 
                        h_back_up     = histo["controlUnit_back"][region][""] 
                        h_back_down   = histo["controlUnit_back"][region][""] 
                    
                    # p_up    = h_up.GetBinContent(b_i)
                    # p_down  = h_down.GetBinContent(b_i)
                    # log_syst_up     = p_up / p
                    # log_syst_down   = p_down / p

            b_i += 1


def run(era, eras, runs_json, syst_json, doRun2, verbose):
    units_json                  = "dc_BkgPred_BinMaps_master.json"
    searchbin_selection_json    = "search_bins_v4.json"
    # Note: DO NOT apply Rz syst. in CR unit bins
    rz_syst_files   = {}
    rz_syst_files["validation"]         = "RzSyst_ValidationBins.root"
    rz_syst_files["validationMetStudy"] = "RzSyst_ValidationBinsMETStudy.root"
    rz_syst_files["search"]             = "RzSyst_SearchBins.root"
    # Note: apply Z vs. Photon syst. in search bins
    ZvPhoton_syst_files = {}
    ZvPhoton_syst_files["validation"]           = "ZvsPhotonSyst_ValidationBins.root"
    ZvPhoton_syst_files["validationMetStudy"]   = "ZvsPhotonSyst_ValidationBinsMETStudy.root"
    ZvPhoton_syst_files["search"]               = "ZvsPhotonSyst_SearchBins.root" 
   
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

    runDir                  = runMap[era]
    result_file             = "condor/" + runDir + "/result.root"
    conf_file               = "datacard_inputs/zinv_syst_" + era + ".conf"
    info_file_znunu         = "syst_info/systematics_znunu_{0}.txt".format(era)
    info_file_phocr_gjets   = "syst_info/systematics_phocr_gjets_{0}.txt".format(era)
    info_file_phocr_back    = "syst_info/systematics_phocr_back_{0}.txt".format(era)
    out_dir                 = "syst_plots/" 
    tmp_dir                 = "tmp_plots/"
    variable                = "mc"
    regions                 = ["lowdm", "highdm"]
    directions              = ["up", "", "down"]
    bintypes                = ["validation", "validationMetStudy", "search", "controlUnit_gjets", "controlUnit_back"]
    # Systematics which we don't use: MET uncluster in photon CR, lepton veto SF, ISR weight for ttbar
    systematics_znunu  = ["jes","btag","pileup","pdf","eff_restoptag","eff_toptag","eff_wtag","eff_fatjet_veto","met_trig","eff_sb","metres"]
    systematics_phocr  = ["jes","btag","pileup","pdf","eff_restoptag","eff_toptag","eff_wtag","eff_fatjet_veto","photon_trig","eff_sb_photon","photon_sf"]
    systematics_map    = {  "jes"               : {"znunu" : "jes",             "phocr" : "jes"},
                            "btag"              : {"znunu" : "btag",            "phocr" : "btag"},
                            "pileup"            : {"znunu" : "pileup",          "phocr" : "pileup"},
                            "pdf"               : {"znunu" : "pdf",             "phocr" : "pdf"},
                            "eff_restoptag"     : {"znunu" : "eff_restoptag",   "phocr" : "eff_restoptag"},
                            "eff_toptag"        : {"znunu" : "eff_toptag",      "phocr" : "eff_toptag"},
                            "eff_wtag"          : {"znunu" : "eff_wtag",        "phocr" : "eff_wtag"},
                            "eff_fatjet_veto"   : {"znunu" : "eff_fatjet_veto", "phocr" : "eff_fatjet_veto"},
                            "met_trig"          : {"znunu" : "met_trig",        "phocr" : None},
                            "photon_trig"       : {"znunu" : None,              "phocr" : "photon_trig"},
                            "eff_sb"            : {"znunu" : "eff_sb",          "phocr" : "eff_sb_photon"},
                            "metres"            : {"znunu" : "metres",          "phocr" : None},
                            "photon_sf"         : {"znunu" : None,              "phocr" : "photon_sf"}
                         }
    # Include prefire here; WARNING: prefire syst. only exists in (2016,2017) and needs to be handled carefully 
    if era != "2018":
        systematics_znunu.append("prefire")
        systematics_phocr.append("prefire")
    systematics = list(set(systematics_znunu) | set(systematics_phocr))

    varTagMap = {
            "jes"    : {"up" : "_jesTotalUp",   "down" : "_jesTotalDown"},
            "metres" : {"up" : "_METUnClustUp", "down" : "_METUnClustDown"}
    }

    # histo[bintype][region][direction]
    histo      = {bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} 

    # syst_histo[systemaitc][bintype][region][direction]
    syst_histo = {syst:{bintype:{region:dict.fromkeys(directions) for region in regions} for bintype in bintypes} for syst in systematics} 

    #-------------------------------------------------------
    # Initiate class instances
    #-------------------------------------------------------

    N       =  Normalization(           tmp_dir, verbose)
    S       =  Shape(                   tmp_dir, draw, doUnits, verbose)
    VB      =  ValidationBins(          N, S, eras, tmp_dir, verbose, draw, saveRootFile)
    VB_MS   =  ValidationBinsMETStudy(  N, S, eras, tmp_dir, verbose, draw, saveRootFile)
    SB      =  SearchBins(              N, S, eras, tmp_dir, verbose, draw, saveRootFile)
    CRU     =  CRUnitBins(              N, S, eras, tmp_dir, verbose)
    
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
        CRU.getValues(r, e)
        VB.getValues(r, e)
        VB_MS.getValues(r, e)
        SB.getValues(r, e, CRunits=CRU)
    
    # histograms for total era being run
    for region in regions:
        histo["validation"][region][""]             =  VB.histograms[era][region][variable].Clone()
        histo["validationMetStudy"][region][""]     =  VB_MS.histograms[era][region][variable].Clone()
        histo["search"][region][""]                 =  SB.histograms[era][region][variable].Clone()
        histo["controlUnit_gjets"][region][""]      =  CRU.histograms[era][region]["mc_gjets"].Clone()
        histo["controlUnit_back"][region][""]       =  CRU.histograms[era][region]["mc_back"].Clone()
    
    #-------------------------------------------------------
    # Calculate normalization and shape factors
    #-------------------------------------------------------
   
    # info files for systematics
    info_znunu          = open(info_file_znunu, "w")
    info_phocr_gjets    = open(info_file_phocr_gjets, "w")
    info_phocr_back     = open(info_file_phocr_back, "w")
    
    with open(conf_file, "w") as outFile:

        
        #-------------------------------------------------------
        # Systematics for znunu
        #-------------------------------------------------------
           
        # --- begin loop over systematics
        for syst in systematics_znunu:
            for direction in ["up", "down"]:
                systTag = "_{0}_syst_{1}".format(syst, direction)
                varTag = ""
                # variable name changes for JEC and MET unclustering systematics
                if syst in varTagMap:
                    varTag = varTagMap[syst][direction]
        
                # --- calculate variation for this systematic
                # only vary Z nu nu; do not vary Normalization or Shape
                VB.getValues(       result_file,  era, systTag, varTag)
                VB_MS.getValues(    result_file,  era, systTag, varTag)
                SB.getValues(       result_file,  era, systTag, varTag, CRunits=CRU)
                
                for region in regions:
                    histo["validation"][region][direction]          =  VB.histograms[era][region][variable].Clone()
                    histo["validationMetStudy"][region][direction]  =  VB_MS.histograms[era][region][variable].Clone()
                    histo["search"][region][direction]              =  SB.histograms[era][region][variable].Clone()
                    
                    # fix prefire because it does not have a weight or syst. in 2018, but 2018 hist. was not added
                    if era == "Run2" and syst == "prefire":
                        #for e in ["2018_PreHEM", "2018_PostHEM"]:
                        for e in ["2018"]:
                            histo["validation"][region][direction].Add(         VB.histograms[e][region][variable]      )
                            histo["validationMetStudy"][region][direction].Add( VB_MS.histograms[e][region][variable]   )
                            histo["search"][region][direction].Add(             SB.histograms[e][region][variable]      )
                    
                    syst_histo[syst]["validation"][region][direction]           = histo["validation"][region][direction]
                    syst_histo[syst]["validationMetStudy"][region][direction]   = histo["validationMetStudy"][region][direction]
                    syst_histo[syst]["search"][region][direction]               = histo["search"][region][direction]
            
            #-------------------------------------------------------
            # Symmetrize systematics which are in the same direction
            #-------------------------------------------------------

            if doSymmetrize:
                for region in regions:
                    # symmetrize systematic if up/down variation is in the same direction compared to nominal
                    # modify histograms passed to function
                    # symmetrizeSyst(h, h_up, h_down)
                    symmetrizeSyst(histo["validation"][region][""],         syst_histo[syst]["validation"][region]["up"],           syst_histo[syst]["validation"][region]["down"])
                    symmetrizeSyst(histo["validationMetStudy"][region][""], syst_histo[syst]["validationMetStudy"][region]["up"],   syst_histo[syst]["validationMetStudy"][region]["down"])
                    symmetrizeSyst(histo["search"][region][""],             syst_histo[syst]["search"][region]["up"],               syst_histo[syst]["search"][region]["down"])
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------

            # WARNING: offset is starting point for low/high dm bins, use with care
            #writeToConfFromPred(outFile, infoFile, binMap, process, syst, h, h_up, h_down, region, offset)
            systForConf = systMap[syst]["name"]  
            writeToConfFromPred(outFile, info_znunu, searchBinMap, "znunu", systForConf, histo["search"]["lowdm"][""],  histo["search"]["lowdm"]["up"],  histo["search"]["lowdm"]["down"],  "lowdm",  SB.low_dm_start)
            writeToConfFromPred(outFile, info_znunu, searchBinMap, "znunu", systForConf, histo["search"]["highdm"][""], histo["search"]["highdm"]["up"], histo["search"]["highdm"]["down"], "highdm", SB.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------
            
            for bintype in ["validation", "validationMetStudy", "search"]:
                for region in regions:
                    # plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir, mc_label)
                    plot(histo[bintype][region][""], histo[bintype][region]["up"], histo[bintype][region]["down"], syst, bintype, region, era, out_dir, "Z#rightarrow#nu#nu MC")
                
        # --- end loop over systematics
        
        #-------------------------------------------------------
        # Systematics for phocr
        #-------------------------------------------------------
        
        # --- begin loop over systematics
        for syst in systematics_phocr:
            for direction in ["up", "down"]:
                systTag = "_{0}_syst_{1}".format(syst, direction)
                varTag  = ""
                # variable name changes for JEC and MET unclustering systematics
                if syst in varTagMap:
                    varTag = varTagMap[syst][direction]
        
                # --- calculate variation for this systematic
                # Do not vary Normalization
                # Vary shape to get GJets MC variation
                S.getShape(         result_file,  era, systTag, varTag)
                print "Using systTag={0}, varTag={1}".format(systTag, varTag)
                CRU.getValues(      result_file,  era)
                
                for region in regions:
                    histo["controlUnit_gjets"][region][direction]   =  CRU.histograms[era][region]["mc_gjets"].Clone()
                    histo["controlUnit_back"][region][direction]    =  CRU.histograms[era][region]["mc_back"].Clone()
                    
                    # fix prefire because it does not have a weight or syst. in 2018, but 2018 hist. was not added
                    if era == "Run2" and syst == "prefire":
                        #for e in ["2018_PreHEM", "2018_PostHEM"]:
                        for e in ["2018"]:
                            histo["controlUnit_gjets"][region][direction].Add( CRU.histograms[e][region]["mc_gjets"] )
                            histo["controlUnit_back"][region][direction].Add(  CRU.histograms[e][region]["mc_back"]  )
                    
                    syst_histo[syst]["controlUnit_gjets"][region][direction]   = histo["controlUnit_gjets"][region][direction]
                    syst_histo[syst]["controlUnit_back"][region][direction]    = histo["controlUnit_back"][region][direction]
            
            #-------------------------------------------------------
            # Symmetrize systematics which are in the same direction
            #-------------------------------------------------------

            if doSymmetrize:
                for region in regions:
                    # symmetrize systematic if up/down variation is in the same direction compared to nominal
                    # modify histograms passed to function
                    # symmetrizeSyst(h, h_up, h_down)
                    symmetrizeSyst(histo["controlUnit_gjets"][region][""],    syst_histo[syst]["controlUnit_gjets"][region]["up"],  syst_histo[syst]["controlUnit_gjets"][region]["down"])
                    symmetrizeSyst(histo["controlUnit_back"][region][""],     syst_histo[syst]["controlUnit_back"][region]["up"],   syst_histo[syst]["controlUnit_back"][region]["down"])
            
            #-------------------------------------------------------
            # Write to conf
            #-------------------------------------------------------

            # WARNING: offset is starting point for low/high dm bins, use with care
            #writeToConfFromPred(outFile, infoFile, binMap, process, syst, h, h_up, h_down, region, offset)
            systForConf = systMap[syst]["name"]  
            writeToConfFromPred(outFile, info_phocr_gjets, unitBinMap,   "phocr_gjets", systForConf, histo["controlUnit_gjets"]["lowdm"][""],  histo["controlUnit_gjets"]["lowdm"]["up"],  histo["controlUnit_gjets"]["lowdm"]["down"],  "lowdm",  CRU.low_dm_start)
            writeToConfFromPred(outFile, info_phocr_gjets, unitBinMap,   "phocr_gjets", systForConf, histo["controlUnit_gjets"]["highdm"][""], histo["controlUnit_gjets"]["highdm"]["up"], histo["controlUnit_gjets"]["highdm"]["down"], "highdm", CRU.high_dm_start)
            writeToConfFromPred(outFile, info_phocr_back,  unitBinMap,   "phocr_back",  systForConf, histo["controlUnit_back"]["lowdm"][""],   histo["controlUnit_back"]["lowdm"]["up"],   histo["controlUnit_back"]["lowdm"]["down"],  "lowdm",  CRU.low_dm_start)
            writeToConfFromPred(outFile, info_phocr_back,  unitBinMap,   "phocr_back",  systForConf, histo["controlUnit_back"]["highdm"][""],  histo["controlUnit_back"]["highdm"]["up"],  histo["controlUnit_back"]["highdm"]["down"], "highdm", CRU.high_dm_start)
            
            #-------------------------------------------------------
            # Plot
            #-------------------------------------------------------

            labels = {
                "controlUnit_gjets" : "#gamma+Jets MC",
                "controlUnit_back"  : "photon background MC"
            }
            
            for bintype in ["controlUnit_gjets", "controlUnit_back"]:
                for region in regions:
                    # plot(h, h_up, h_down, mySyst, bintype, region, era, plot_dir, mc_label)
                    plot(histo[bintype][region][""], histo[bintype][region]["up"], histo[bintype][region]["down"], syst, bintype, region, era, out_dir, labels[bintype])
                
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
            # Note: DO NOT apply Z vs Photon syst. in CR unit bins
            if bintype != "controlUnit_gjets" and bintype != "controlUnit_back":
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
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["lowdm"]["znunu_rzunc"],  "lowdm",  SB.low_dm_start,  searchBinSelectionMap, ["NJ"])
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["highdm"]["znunu_rzunc"], "highdm", SB.high_dm_start, searchBinSelectionMap, ["NJ"])
       
        # --- Z vs Photon syst --- #
        # writeToConfFromSyst(outFile, binMap, process, syst, h, region, offset, selectionMap = {}, removeCut = "")
        systForConf = systMap["znunu_zgammdiff"]["name"]  
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["lowdm"]["znunu_zgammdiff"],  "lowdm",  SB.low_dm_start)
        writeToConfFromSyst(outFile, searchBinMap, "znunu", systForConf, systHistoMap["search"]["highdm"]["znunu_zgammdiff"], "highdm", SB.high_dm_start)


    # Total systematics function which can run on validaiton, MET study, search bins, CR unit bins, etc.
    
    # getTotalSystematics(BinObject, bintype, systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2)
    getTotalSystematics(VB,     "validation",           systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2)
    getTotalSystematics(VB_MS,  "validationMetStudy",   systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2)
    getTotalSystematics(SB,     "search",               systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2)
    # CR unit bin not supported

    # total systematics using variation on full prediction
    # getTotalSystematicsPrediction(SearchBinObject, CRBinObject, N, S, runMap, systematics_znunu, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2):
    # getTotalSystematicsPrediction(SearchBinObject, CRBinObject, N, S, runMap, systematics_map, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2):
    getTotalSystematicsPrediction(SB, CRU, N, S, runMap, systematics_map, systHistoMap, histo, syst_histo, era, directions, regions, out_dir, doRun2)

    # close files
    info_znunu.close()
    info_phocr_gjets.close()
    info_phocr_back.close()

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--runs_json",    "-j", default="",                             help="json file containing runs")
    parser.add_argument("--syst_json",    "-s", default="systematics.json",             help="json file containing systematics")
    parser.add_argument("--verbose",      "-v", default = False, action = "store_true", help="verbose flag to print more things")
    
    options                     = parser.parse_args()
    runs_json                   = options.runs_json
    syst_json                   = options.syst_json
    verbose                     = options.verbose
    doRun2                      = True
    
    #eras = ["2016", "2017", "2018", "Run2"]
    eras = ["2016", "2017", "2018", "2016and2017", "Run2"]
    for era in eras:
        # run(era, eras, runs_json, syst_json, doRun2, verbose)
        run(era, eras, runs_json, syst_json, doRun2, verbose)


if __name__ == "__main__":
    main()



