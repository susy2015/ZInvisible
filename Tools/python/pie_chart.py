# pie_chart.py

import ROOT
import json
import numpy as np
import matplotlib.pyplot as plt

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def pieChart(inputFile):
    f = ROOT.TFile(inputFile, "read")
    h_ttbar = f.Get("httbar_pred")
    h_znunu = f.Get("znunu_pred")
    h_qcd   = f.Get("qcd_pred")
    h_ttz   = f.Get("ttZ_pred")
    h_rare  = f.Get("Rare_pred")
    ttbar_lowdm  = 0
    znunu_lowdm  = 0
    qcd_lowdm    = 0
    ttz_lowdm    = 0
    rare_lowdm   = 0
    ttbar_highdm = 0
    znunu_highdm = 0
    qcd_highdm   = 0
    ttz_highdm   = 0
    rare_highdm  = 0
    lowdm_range  = [0,  52]
    highdm_range = [53, 182]
    
    #print "--- lowdm ---"
    for b in xrange(lowdm_range[0], lowdm_range[1] + 1):
        i = b + 1
        #print "i = {0}, b = {1}".format(i, b)
        ttbar_lowdm += h_ttbar.GetBinContent(i)
        znunu_lowdm += h_znunu.GetBinContent(i)
        qcd_lowdm   += h_qcd.GetBinContent(i)
        ttz_lowdm   += h_ttz.GetBinContent(i)
        rare_lowdm  += h_rare.GetBinContent(i)
    
    #print "--- highdm ---"
    for b in xrange(highdm_range[0], highdm_range[1] + 1):
        i = b + 1
        #print "i = {0}, b = {1}".format(i, b)
        ttbar_highdm += h_ttbar.GetBinContent(i)
        znunu_highdm += h_znunu.GetBinContent(i)
        qcd_highdm   += h_qcd.GetBinContent(i)
        ttz_highdm   += h_ttz.GetBinContent(i)
        rare_highdm  += h_rare.GetBinContent(i)

    ttbar_total = ttbar_lowdm + ttbar_highdm
    znunu_total = znunu_lowdm + znunu_highdm
    qcd_total   = qcd_lowdm + qcd_highdm
    ttz_total   = ttz_lowdm + ttz_highdm
    rare_total  = rare_lowdm + rare_highdm

    print "low dm:  ttbar = {0:.3f}, znunu = {1:.3f}, qcd = {2:.3f}, ttz = {3:.3f}, rare = {4:.3f}".format(ttbar_lowdm, znunu_lowdm, qcd_lowdm, ttz_lowdm, rare_lowdm)
    print "high dm: ttbar = {0:.3f}, znunu = {1:.3f}, qcd = {2:.3f}, ttz = {3:.3f}, rare = {4:.3f}".format(ttbar_highdm, znunu_highdm, qcd_highdm, ttz_highdm, rare_highdm)
    print "total:   ttbar = {0:.3f}, znunu = {1:.3f}, qcd = {2:.3f}, ttz = {3:.3f}, rare = {4:.3f}".format(ttbar_total, znunu_total, qcd_total, ttz_total, rare_total)


    # Pie chart, where the slices will be ordered and plotted counter-clockwise

    # xkcd colors
    with open("rgb.json", "r") as j:
        colorMap = json.load(j) 

    labels = ["ttbar", "znunu", "qcd", "ttz", "rare"]
    colors = [colorMap["blue"], colorMap["orange"], colorMap["green"], colorMap["red"], colorMap["purple"]]
    #colors = ['xkcd:blue', 'xkcd:orange', 'xkcd:green', 'xkcd:red', 'xkcd:purple']
    #colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    #colors = ['b', 'y', 'g', 'r', 'm']
    
    # lowdm
    sizes_lowdm  = [ttbar_lowdm, znunu_lowdm, qcd_lowdm, ttz_lowdm, rare_lowdm]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes_lowdm, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title("Total Background Predictions: Low $\Delta m$")
    plt.tight_layout()
    plt.savefig("more_plots/backgrounds_lowdm.pdf", bbox_inches="tight")
    plt.savefig("more_plots/backgrounds_lowdm.png", bbox_inches="tight")
    
    # highdm
    sizes_highdm  = [ttbar_highdm, znunu_highdm, qcd_highdm, ttz_highdm, rare_highdm]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes_highdm, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title("Total Background Predictions: High $\Delta m$")
    plt.tight_layout()
    plt.savefig("more_plots/backgrounds_highdm.pdf", bbox_inches="tight")
    plt.savefig("more_plots/backgrounds_highdm.png", bbox_inches="tight")
    
    # total
    sizes_total  = [ttbar_total, znunu_total, qcd_total, ttz_total, rare_total]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes_total, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title("Total Background Predictions: Full Search Region")
    plt.tight_layout()
    plt.savefig("more_plots/backgrounds_total.pdf", bbox_inches="tight")
    plt.savefig("more_plots/backgrounds_total.png", bbox_inches="tight")

def main():
    inputFile = "/uscms/home/mkilpatr/nobackup/CMSSW_9_4_10/src/AnalysisMethods/EstTools/SUSYNano19/getFinalPlot_allMethods/pred_binnum_getFinalPlot_Nano.root"
    pieChart(inputFile)

if __name__ == "__main__":
    main()


