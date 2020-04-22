# compare.py
# compare histograms

import ROOT
from tools import setupHist

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# return histogram map
def load(f, label, era):
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

def plot(h_map_1, h_map_2, era):
    
    plot_dir = "comparison_plots"
    eraTag = "_" + era
    regions = ["lowdm", "highdm"]
    draw_option = "hist error"
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
            
            h1 = h_map_1[region][value]
            h2 = h_map_2[region][value]
            h_ratio = h2.Clone("h_ratio") 
            h_ratio.Divide(h1)
            
            title_main  = "Z to Invisible {0} {1}: Compare {2} and {3}".format(value, era, label1, label2)
            title_ratio = "Z to Invisible {0}: {1}".format(value, label_ratio)
            x_title = "Validation Bin"
            
            #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
            setupHist(h1, title_main, x_title, "Events", color_red,  10.0 ** -2, 10.0 ** 4)
            setupHist(h2, title_main, x_title, "Events", color_blue, 10.0 ** -2, 10.0 ** 4)
            setupHist(h_ratio, title_ratio, x_title, label_ratio, color_blue, 0.5, 1.5)
            
            # histograms
            c.cd(1)
            ROOT.gPad.SetLogy(1) # set log y
            h1.Draw(draw_option)
            h2.Draw("error same")
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

# compare validation histograms
def validation(file_map, era):
    label1 = "Caleb"
    label2 = "Angel"
    file1 = file_map[label1]
    file2 = file_map[label2]
    f1    = ROOT.TFile(file1, "read")
    f2    = ROOT.TFile(file2, "read")
    # setup histogram maps
    h_map_1 = {}
    h_map_1["label"] = label1
    h_map_1["lowdm"]  = {}
    h_map_1["highdm"] = {}
    # load histograms 
    h_map_1["lowdm"]["data"]  = f1.Get("data_lowdm")
    h_map_1["lowdm"]["mc"]    = f1.Get("mc_lowdm")
    h_map_1["lowdm"]["pred"]  = f1.Get("pred_lowdm")
    h_map_1["highdm"]["data"] = f1.Get("data_highdm")
    h_map_1["highdm"]["mc"]   = f1.Get("mc_highdm")
    h_map_1["highdm"]["pred"] = f1.Get("pred_highdm")
    h_map_2 = load(f2, label2, era)

    # make plots
    plot(h_map_1, h_map_2, era)


if __name__ == "__main__":
    eras = ["2016", "2017", "2018", "Run2"]
    #caleb_date = "2020-03-09"
    #caleb_date = "2020-03-19"
    caleb_date = "2020-03-24"
    #caleb_date = "2020-04-15"
    angel_date = "2020-04-01"
    for era in eras:
        f_map = {}
        f_map["Caleb"] = "/uscms/home/caleb/archive/zinv_results/{0}/prediction_histos/validationBinsZinv_{1}.root".format(caleb_date, era)
        f_map["Angel"] = "/uscms/home/caleb/archive/angel/{0}/{1}/result.root".format(angel_date, era)
        validation(f_map, era)

