# combine_hists.py

import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Example
# input files: searchBinsZinv_2016.root, searchBinsZinv_syst_2016.root
#
# TFile**     searchBinsZinv_2016.root    
#  TFile*     searchBinsZinv_2016.root    
#   KEY: TH1F data_lowdm;1    Z to Invisible Data, MC and Prediction for 2016
#   KEY: TH1F mc_lowdm;1  Z to Invisible Data, MC and Prediction for 2016
#   KEY: TH1F pred_lowdm;1    Z to Invisible Data, MC and Prediction for 2016
#   KEY: TH1F data_highdm;1   Z to Invisible Data, MC and Prediction for 2016
#   KEY: TH1F mc_highdm;1 Z to Invisible Data, MC and Prediction for 2016
#   KEY: TH1F pred_highdm;1   Z to Invisible Data, MC and Prediction for 2016
# 
# TFile**     searchBinsZinv_syst_2016.root   
#  TFile*     searchBinsZinv_syst_2016.root   
#   KEY: TH1F syst_up_lowdm;1 syst_up_lowdm
#   KEY: TH1F syst_down_lowdm;1   syst_down_lowdm
#   KEY: TH1F syst_up_highdm;1    syst_up_highdm
#   KEY: TH1F syst_down_highdm;1  syst_down_highdm


def combine(predFile, systFile, outFile):
    f_pred        = ROOT.TFile(predFile, "read")
    f_syst        = ROOT.TFile(systFile, "read")
    f_out         = ROOT.TFile(outFile,  "recreate")
    # data, mc, pred
    h_data_lowdm  = f_pred.Get("data_lowdm")
    h_mc_lowdm    = f_pred.Get("mc_lowdm")
    h_pred_lowdm  = f_pred.Get("pred_lowdm")
    h_data_highdm = f_pred.Get("data_highdm")
    h_mc_highdm   = f_pred.Get("mc_highdm")
    h_pred_highdm = f_pred.Get("pred_highdm")
    # syst
    h_syst_up_lowdm    = f_syst.Get("syst_up_lowdm")
    h_syst_down_lowdm  = f_syst.Get("syst_down_lowdm")
    h_syst_up_highdm   = f_syst.Get("syst_up_highdm")
    h_syst_down_highdm = f_syst.Get("syst_down_highdm")
    # search bins 
    low_dm_start   = 0
    low_dm_end     = 52
    high_dm_start  = 53
    high_dm_end    = 182
    combined_start = low_dm_start
    combined_end   = high_dm_end
    low_dm_nbins   = low_dm_end - low_dm_start + 1 
    high_dm_nbins  = high_dm_end - high_dm_start + 1 
    combined_nbins = combined_end - combined_start + 1
    
    h_data_combined  = ROOT.TH1F("data",  "data",  combined_nbins,  combined_start,  combined_end + 1) 
    h_mc_combined    = ROOT.TH1F("mc",    "mc",    combined_nbins,  combined_start,  combined_end + 1) 
    h_pred_combined  = ROOT.TH1F("pred",  "pred",  combined_nbins,  combined_start,  combined_end + 1) 

    if not h_data_lowdm:
        print "ERROR: No data lowdm for {0}".format(predFile)
    if not h_data_highdm:
        print "ERROR: No data highdm for {0}".format(predFile)

    bin_i = 1
    for b in xrange(low_dm_start, low_dm_end + 1):
        h_data_combined.SetBinContent(b + 1, h_data_lowdm.GetBinContent(bin_i)) 
        h_data_combined.SetBinError(b + 1,   h_data_lowdm.GetBinError(bin_i)) 
        h_mc_combined.SetBinContent(b + 1,   h_mc_lowdm.GetBinContent(bin_i)) 
        h_mc_combined.SetBinError(b + 1,     h_mc_lowdm.GetBinError(bin_i)) 
        h_pred_combined.SetBinContent(b + 1, h_pred_lowdm.GetBinContent(bin_i)) 
        h_pred_combined.SetBinError(b + 1,   h_pred_lowdm.GetBinError(bin_i)) 
        bin_i += 1
    bin_i = 1
    for b in xrange(high_dm_start, high_dm_end + 1):
        h_data_combined.SetBinContent(b + 1, h_data_highdm.GetBinContent(bin_i)) 
        h_data_combined.SetBinError(b + 1,   h_data_highdm.GetBinError(bin_i)) 
        h_mc_combined.SetBinContent(b + 1,   h_mc_highdm.GetBinContent(bin_i)) 
        h_mc_combined.SetBinError(b + 1,     h_mc_highdm.GetBinError(bin_i)) 
        h_pred_combined.SetBinContent(b + 1, h_pred_highdm.GetBinContent(bin_i)) 
        h_pred_combined.SetBinError(b + 1,   h_pred_highdm.GetBinError(bin_i)) 
        bin_i += 1

    h_data_combined.Write()
    h_mc_combined.Write()
    h_pred_combined.Write()
    f_out.Close()

def main():
    # WARNING: assumes data, pred, and mc exist for each era
    #eras = ["2016", "2017", "2018", "Run2"]
    eras = ["2016"]
    for era in eras:
        predFile = "prediction_histos/searchBinsZinv_{0}.root".format(era)
        systFile = "prediction_histos/searchBinsZinv_syst_{0}.root".format(era)
        outFile  = "prediction_histos/searchBinsZinv_combined_{0}.root".format(era)
        combine(predFile, systFile, outFile)

if __name__ == "__main__":
    main()


