# compare.py
# compare histograms

import ROOT
from tools import setupHist

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# plot comparison of histograms
def plot(h_map_1, h_map_2, ratio_limits, bin_type, era):
    
    debug = False
    plot_dir = "comparison_plots"
    eraTag = "_" + era
    regions = ["lowdm", "highdm"]
    draw_option = "hist"
    color_red    = "vermillion"
    color_blue   = "electric blue"
    label1 = h_map_1["label"]
    label2 = h_map_2["label"]
    label_ratio = "{0}/{1}".format(label2, label1)
    
    ###################
    # Draw Histograms #
    ###################

    # draw histograms
    c = ROOT.TCanvas("c", "c", 800, 800)
    c.Divide(1, 2)
    
    # legend: TLegend(x1,y1,x2,y2)
    legend_x1 = 0.7
    legend_x2 = 0.9 
    legend_y1 = 0.7 
    legend_y2 = 0.9 
    
    for region in regions:
        for value in h_map_1[region]:
            
            if debug:
                print "{0} {1}".format(region, value)
            h1 = h_map_1[region][value]
            h2 = h_map_2[region][value]
            h_ratio = h2.Clone("h_ratio") 
            h_ratio.Divide(h1)
            
            title_main  = "Z Invisible {0} {1}: Compare {2} and {3}".format(value, era, label1, label2)
            title_ratio = "Z Invisible {0}: {1}".format(value, label_ratio)
            x_title = bin_type + " bin"
            
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h1, title_main, x_title, "Events", color_red,  10.0 ** -2, 10.0 ** 4)
            setupHist(h2, title_main, x_title, "Events", color_blue, 10.0 ** -2, 10.0 ** 4)
            setupHist(h_ratio, title_ratio, x_title, label_ratio, color_blue, ratio_limits[0], ratio_limits[1])
            
            # histograms
            c.cd(1)
            ROOT.gPad.SetLogy(1) # set log y
            h1.Draw(draw_option)
            h2.Draw(draw_option + " same")
            # legend: TLegend(x1,y1,x2,y2)
            legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
            legend.AddEntry(h1, "{0}: {1}".format(value, label1), "l")
            legend.AddEntry(h2, "{0}: {1}".format(value, label2), "l")
            legend.Draw()
           
            # ratios
            c.cd(2)
            h_ratio.Draw(draw_option)
                
            # save histograms
            plot_name = "{0}/{1}_{2}_{3}_{4}".format(plot_dir, region, value, label1, label2)
            c.Update()
            c.SaveAs(plot_name + eraTag + ".pdf")
            c.SaveAs(plot_name + eraTag + ".png")

# return histogram map
def getValidationHists(f, label):
    h_map = {}
    h_map["label"]  = label
    h_map["lowdm"]  = {}
    h_map["highdm"] = {}
    
    # --- Angel's histograms --- #
    #
    # nValidationBinLowDM_jetpt30
    #   data: MET_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30Data MET Validation Bin Low DMdata
    #   mc:   ZNuNu_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata
    #   pred: ZNuNu_nValidationBin_LowDM_nj_shape_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata
    #
    # nValidationBinLowDMHighMET_jetpt30
    #   data: MET_nValidationBin_LowDM_HighMET_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30Data MET Validation Bin Low DM High METdata
    #   mc:   ZNuNu_nValidationBin_LowDM_HighMET_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30ZJetsToNuNu Validation Bin Low DM High METdata
    #   pred: ZNuNu_nValidationBin_LowDM_HighMET_nj_shape_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30ZJetsToNuNu Validation Bin Low DM High METdata
    #
    # nValidationBinHighDM_jetpt30
    #   data: MET_nValidationBin_HighDM_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30Data MET Validation Bin High DMdata
    #   mc:   ZNuNu_nValidationBin_HighDM_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30ZJetsToNuNu Validation Bin High DMdata
    #   pred: ZNuNu_nValidationBin_HighDM_nj_shape_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30ZJetsToNuNu Validation Bin High DMdata
    #
    # -------------------------- #
    
    h_lowdm_data            = f.Get("nValidationBinLowDM_jetpt30/MET_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30Data MET Validation Bin Low DMdata")
    h_lowdm_mc              = f.Get("nValidationBinLowDM_jetpt30/ZNuNu_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata")
    h_lowdm_pred            = f.Get("nValidationBinLowDM_jetpt30/ZNuNu_nValidationBin_LowDM_nj_shape_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata")
    h_lowdm_highmet_data    = f.Get("nValidationBinLowDMHighMET_jetpt30/MET_nValidationBin_LowDM_HighMET_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30Data MET Validation Bin Low DM High METdata")
    h_lowdm_highmet_mc      = f.Get("nValidationBinLowDMHighMET_jetpt30/ZNuNu_nValidationBin_LowDM_HighMET_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30ZJetsToNuNu Validation Bin Low DM High METdata")
    h_lowdm_highmet_pred    = f.Get("nValidationBinLowDMHighMET_jetpt30/ZNuNu_nValidationBin_LowDM_HighMET_nj_shape_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30ZJetsToNuNu Validation Bin Low DM High METdata")
    h_highdm_data           = f.Get("nValidationBinHighDM_jetpt30/MET_nValidationBin_HighDM_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30Data MET Validation Bin High DMdata")
    h_highdm_mc             = f.Get("nValidationBinHighDM_jetpt30/ZNuNu_nValidationBin_HighDM_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30ZJetsToNuNu Validation Bin High DMdata")
    h_highdm_pred           = f.Get("nValidationBinHighDM_jetpt30/ZNuNu_nValidationBin_HighDM_nj_shape_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30ZJetsToNuNu Validation Bin High DMdata")
       
    low_dm_start           = 0
    low_dm_normal_end      = 14
    low_dm_highmet_start   = 15
    low_dm_end             = 18
    low_dm_nbins           = low_dm_end - low_dm_start + 1 
    h_lowdm_data_total     = ROOT.TH1F("lowdm_data",    "lowdm_data",    low_dm_nbins,  low_dm_start,  low_dm_end + 1) 
    h_lowdm_mc_total       = ROOT.TH1F("lowdm_mc",      "lowdm_mc",      low_dm_nbins,  low_dm_start,  low_dm_end + 1) 
    h_lowdm_pred_total     = ROOT.TH1F("lowdm_pred",    "lowdm_pred",    low_dm_nbins,  low_dm_start,  low_dm_end + 1) 

    # combine lowdm and lowdm_highmet histograms
    bin_i = 1
    for b in xrange(low_dm_start, low_dm_normal_end + 1):
        h_lowdm_data_total.SetBinContent(b + 1, h_lowdm_data.GetBinContent(bin_i))
        h_lowdm_data_total.SetBinError(b + 1,   h_lowdm_data.GetBinError(bin_i))
        h_lowdm_mc_total.SetBinContent(b + 1,   h_lowdm_mc.GetBinContent(bin_i))
        h_lowdm_mc_total.SetBinError(b + 1,     h_lowdm_mc.GetBinError(bin_i))
        h_lowdm_pred_total.SetBinContent(b + 1, h_lowdm_pred.GetBinContent(bin_i))
        h_lowdm_pred_total.SetBinError(b + 1,   h_lowdm_pred.GetBinError(bin_i))
        bin_i += 1
    bin_i = 1
    for b in xrange(low_dm_highmet_start, low_dm_end + 1):
        h_lowdm_data_total.SetBinContent(b + 1, h_lowdm_highmet_data.GetBinContent(bin_i))
        h_lowdm_data_total.SetBinError(b + 1,   h_lowdm_highmet_data.GetBinError(bin_i))
        h_lowdm_mc_total.SetBinContent(b + 1,   h_lowdm_highmet_mc.GetBinContent(bin_i))
        h_lowdm_mc_total.SetBinError(b + 1,     h_lowdm_highmet_mc.GetBinError(bin_i))
        h_lowdm_pred_total.SetBinContent(b + 1, h_lowdm_highmet_pred.GetBinContent(bin_i))
        h_lowdm_pred_total.SetBinError(b + 1,   h_lowdm_highmet_pred.GetBinError(bin_i))
        bin_i += 1


    h_map["lowdm"]["data"]  = h_lowdm_data_total
    h_map["lowdm"]["mc"]    = h_lowdm_mc_total
    h_map["lowdm"]["pred"]  = h_lowdm_pred_total
    h_map["highdm"]["data"] = h_highdm_data
    h_map["highdm"]["mc"]   = h_highdm_mc
    h_map["highdm"]["pred"] = h_highdm_pred
    
    return h_map


def getSearchHistsSimple(f, h_type, h_name):
    # e.g. h_type = "data", h_name = "hdata"
    h_map = {}
    h_map["lowdm"]  = {}
    h_map["highdm"] = {}
    h = f.Get(h_name)
    
    low_dm_start   = 0
    low_dm_end     = 52
    high_dm_start  = 53
    high_dm_end    = 182
    low_dm_nbins   = low_dm_end - low_dm_start + 1 
    high_dm_nbins  = high_dm_end - high_dm_start + 1 
    h_lowdm        = ROOT.TH1F("lowdm_" + h_type,    "lowdm_" + h_type,    low_dm_nbins,  low_dm_start,  low_dm_end + 1) 
    h_highdm       = ROOT.TH1F("highdm_" + h_type,   "highdm_" + h_type,   high_dm_nbins, high_dm_start, high_dm_end + 1) 
    
    # fill low and high dm histograms
    bin_i = 1
    for b in xrange(low_dm_start, low_dm_end + 1):
        h_lowdm.SetBinContent(bin_i,   h.GetBinContent(b + 1))
        h_lowdm.SetBinError(bin_i,     h.GetBinError(b + 1))
        bin_i += 1
    bin_i = 1
    for b in xrange(high_dm_start, high_dm_end + 1):
        h_highdm.SetBinContent(bin_i,   h.GetBinContent(b + 1))
        h_highdm.SetBinError(bin_i,     h.GetBinError(b + 1))
        bin_i += 1
    
    h_map["lowdm"]  = h_lowdm
    h_map["highdm"] = h_highdm

    return h_map

# return histogram map
def getSearchHists(f, label, h_names):
    h_map = {}
    h_map["label"]  = label
    h_map["lowdm"]  = {}
    h_map["highdm"] = {}
    
    h_data = f.Get("hdata")
    h_pred = f.Get("hznunu")
    
    low_dm_start   = 0
    low_dm_end     = 52
    high_dm_start  = 53
    high_dm_end    = 182
    low_dm_nbins   = low_dm_end - low_dm_start + 1 
    high_dm_nbins  = high_dm_end - high_dm_start + 1 
    h_lowdm_data    = ROOT.TH1F("lowdm_data",    "lowdm_data",    low_dm_nbins,  low_dm_start,  low_dm_end + 1) 
    h_lowdm_pred    = ROOT.TH1F("lowdm_pred",    "lowdm_pred",    low_dm_nbins,  low_dm_start,  low_dm_end + 1) 
    h_highdm_data   = ROOT.TH1F("highdm_data",   "highdm_data",   high_dm_nbins, high_dm_start, high_dm_end + 1) 
    h_highdm_pred   = ROOT.TH1F("highdm_pred",   "highdm_pred",   high_dm_nbins, high_dm_start, high_dm_end + 1) 
    
    # fill low and high dm histograms
    bin_i = 1
    for b in xrange(low_dm_start, low_dm_end + 1):
        h_lowdm_data.SetBinContent(bin_i,   h_data.GetBinContent(b + 1))
        h_lowdm_data.SetBinError(bin_i,     h_data.GetBinError(b + 1))
        h_lowdm_pred.SetBinContent(bin_i,   h_pred.GetBinContent(b + 1))
        h_lowdm_pred.SetBinError(bin_i,     h_pred.GetBinError(b + 1))
        bin_i += 1
    bin_i = 1
    for b in xrange(high_dm_start, high_dm_end + 1):
        h_highdm_data.SetBinContent(bin_i,   h_data.GetBinContent(b + 1))
        h_highdm_data.SetBinError(bin_i,     h_data.GetBinError(b + 1))
        h_highdm_pred.SetBinContent(bin_i,   h_pred.GetBinContent(b + 1))
        h_highdm_pred.SetBinError(bin_i,     h_pred.GetBinError(b + 1))
        bin_i += 1
    
    h_map["lowdm"]["data"]      = h_lowdm_data
    h_map["lowdm"]["pred"]      = h_lowdm_pred
    h_map["highdm"]["data"]     = h_highdm_data
    h_map["highdm"]["pred"]     = h_highdm_pred
    
    return h_map

# return histogram map
def getMyHists(f, label, values):
    regions = ["lowdm", "highdm"]
    h_map = {}
    h_map["label"] = label
    for region in regions:
        h_map[region] = {}
        for value in values:
            name = "{0}_{1}".format(value, region)
            h_map[region][value] = f.Get(name)
    return h_map

# compare validation bin histograms
def validation(file_map, labels, era):
    label1  = labels[0]
    label2  = labels[1]
    values  = ["data", "mc", "pred"]
    file1   = file_map[label1]
    file2   = file_map[label2]
    f1      = ROOT.TFile(file1, "read")
    f2      = ROOT.TFile(file2, "read")
    h_map_1 = getMyHists(           f1, label1 , values )
    h_map_2 = getValidationHists(   f2, label2 )

    # make plots
    print "Running plot() for {0} and {1}, {2}".format(label1, label2, era)
    ratio_limits = [0.5, 1.5]
    plot(h_map_1, h_map_2, ratio_limits, "validation", era)

# compare search bin histograms
def search(file_map, values, labels, h_names, era):
    label1  = labels[0]
    label2  = labels[1]
    file1   = file_map[label1]
    file2   = file_map[label2]
    f1      = ROOT.TFile(file1, "read")
    f2      = ROOT.TFile(file2, "read")
    h_map_1 = getMyHists(       f1, label1, values )
    
    h_map_2 = {}
    h_map_2["label"]  = label2
    h_map_2["lowdm"]  = {}
    h_map_2["highdm"] = {}

    for value in values:
        h_name = value
        if value == "pred":
            h_name = "znunu_pred"
    
        simple_map = getSearchHistsSimple(f2, value, h_names[h_name])
        h_map_2["lowdm"][value]  = simple_map["lowdm"]
        h_map_2["highdm"][value] = simple_map["highdm"] 

    # make plots
    print "Running plot() for {0} and {1}, {2}".format(label1, label2, era)
    ratio_limits = [0.98, 1.02]
    plot(h_map_1, h_map_2, ratio_limits, "search", era)

def searchCompareMyHists(file_map, values, labels, era):
    label1 = labels[0]
    label2 = labels[1]
    file1 = file_map[label1]
    file2 = file_map[label2]
    f1    = ROOT.TFile(file1, "read")
    f2    = ROOT.TFile(file2, "read")
    h_map_1 = getMyHists(       f1, label1, values )
    h_map_2 = getMyHists(       f2, label2, values )

    # make plots
    print "Running plot() for {0} and {1}, {2}".format(label1, label2, era)
    ratio_limits = [0.5, 1.5]
    plot(h_map_1, h_map_2, ratio_limits, "search", era)

if __name__ == "__main__":
    
    eras = ["2016", "2017", "2018", "Run2"]
    
    # ----------------------- #
    # --- validation bins --- #
    # ----------------------- #
    
    # # v6 ntuples, old top/w weights
    # #caleb_date = "2020-03-09"
    # #caleb_date = "2020-03-19"
    # caleb_date = "2020-03-24"
    # 
    # # v6 ntuples, new top/w weights
    # #caleb_date = "2020-04-15"
    # 
    # angel_date = "2020-04-01"
    # 
    # labels = ["Caleb", "Angel"]
    # 
    # for era in eras:
    #     f_map = {}
    #     f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/{0}/prediction_histos/validationBinsZinv_{1}.root".format(caleb_date, era)
    #     f_map["Angel"] = "/uscms/home/caleb/archive/angel/{0}/{1}/result.root".format(angel_date, era)
    #     validation(f_map, labels, era)
    
    
    # ------------------- #
    # --- search bins --- #
    # ------------------- #
    
    # Note: Matt only provided 2016 for now; do not use other eras
    # Matt's files:
    # /uscms/home/caleb/archive/matt/2020-04-22/2016/SumOfBkg_T1tttt_2000_0.root
    # /uscms/home/caleb/archive/matt/2020-04-22/2016/SumOfBkg_T2tt_1000_0.root
    
    # # v6 ntuples, old top/w weights
    # # note: these runs do not have MET data histograms... I am not sure why
    # # note: Matt still uses 09Mar2020_dev_v6 inputs which match my 2020-03-24 inputs
    # #caleb_date = "2020-03-19"
    # caleb_date = "2020-03-24"
    # 
    # # v6 ntuples, new top/w weights
    # #caleb_date = "2020-04-09"
    # #caleb_date = "2020-04-15"
    # 
    # matt_date  = "2020-04-22"
    # 
    # for era in ["2016"]:
    #     f_map = {}
    #     f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/{0}/prediction_histos/searchBinsZinv_{1}.root".format(caleb_date, era)
    #     f_map["Matt"]  = "/uscms/home/caleb/archive/matt/{0}/{1}/SumOfBkg_T1tttt_2000_0.root".format(matt_date, era)
    #     search(f_map, era)
    
    # --- Matt's histograms --- #
    #
    # data:         hdata
    # total pred:   hpred
    # znunu pred:   hznunu
    #
    # -------------------------- #
    
    # v6p5 ntuple unblinding Run2 comparison
    # comparison with Matt
    # Run2
    values = ["data", "pred"]
    labels = ["Caleb", "Matt"]
    h_names = {"data" : "hdata", "znunu_pred" : "hznunu"}
    f_map = {}
    #f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/2020-05-19/prediction_histos/searchBinsZinv_Run2.root"
    #f_map["Matt"]  = "/eos/uscms/store/user/lpcsusyhad/Stop_production/LimitInputs/19May2020_Run2Unblind_dev_v6/SearchBinsPlot/SumOfBkg.root"
    f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/2020-07-06/prediction_histos/searchBinsZinv_Run2.root"
    f_map["Matt"]  = "/uscms/home/mkilpatr/nobackup/CMSSW_9_4_10/src/AnalysisMethods/EstTools/SUSYNano19/getFinalPlot/SumOfBkg.root"
    search(f_map, values, labels, h_names, "Run2")

    # compare systematics with Matt
    values = ["syst_up", "syst_down"]
    labels = ["Caleb", "Matt"]
    h_names = {"syst_up" : "znunu_syst_up", "syst_down" : "znunu_syst_dn"}
    f_map = {}
    f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/2020-07-06/prediction_histos/searchBinsZinv_totalPredSyst_Run2.root"
    f_map["Matt"]  = "/uscms/home/mkilpatr/nobackup/CMSSW_9_4_10/src/AnalysisMethods/EstTools/SUSYNano19/getFinalPlot/SumOfBkg.root"
    search(f_map, values, labels, h_names, "Run2")

    # --- Jon's histograms --- #
    #
    # data:         Data
    # QCD pred:     QCD
    #
    # -------------------------- #
    
    # comparison with Jon
    # /uscms_data/d3/jsw/SusyAna/CMSSW_10_2_9/src/notebooks/Yields_v6_5/2016_SearchBins_QCDResidMET.root
    # /eos/uscms/store/user/lpcsusyhad/Stop_production/LimitInputs/26May2020_2016Unblind_dev_v6/SearchBinsPlot/SearchBins_QCDResidMET.root
    # 2016
    values = ["data"]
    labels = ["Caleb", "Jon"]
    h_names = {"data" : "Data", "qcd_pred" : "QCD"}
    f_map = {}
    f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/2020-05-19/prediction_histos/searchBinsZinv_2016.root"
    f_map["Jon"]   = "/eos/uscms/store/user/lpcsusyhad/Stop_production/LimitInputs/26May2020_2016Unblind_dev_v6/SearchBinsPlot/SearchBins_QCDResidMET.root"
    search(f_map, values, labels, h_names, "2016")

    
    # --- Caleb's comparison --- #
    #
    # compare using Rz per year with using Run 2 Rz
    #
    # -------------------------- #
    values = ["pred"]
    labels = ["RzRun2", "RzPerYear"]
    f_map = {}
    f_map["RzRun2"]    = "/uscms_data/d3/caleb/SusyAnalysis/CMSSW_10_2_9/src/ZInvisible/Tools/prediction_histos/searchBinsZinv_Run2.root"
    f_map["RzPerYear"] = "/uscms_data/d3/caleb/SusyAnalysis/CMSSW_10_2_9/src/ZInvisible/Tools/prediction_histos/searchBinsZinv_useRzPerYear_Run2.root"
    searchCompareMyHists(f_map, values, labels, "Run2")
    


