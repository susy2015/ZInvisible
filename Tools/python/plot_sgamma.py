# plot_sgamma.py

import ROOT
import json
from tools import setupHist

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def run(era):
    verbose = 1
    print "---------- Running {0} ----------".format(era)
    # for datacard, the allowed Sgamma range is [0.01, 5]
    minSgamma = 0.01
    maxSgamma = 5.0
    
    fileSB = "results/SearchBinResults.json"
    fileYields = "datacard_inputs/zinv_yields_" + era + ".json"
    fileBinMap = "dc_BkgPred_BinMaps_master.json"
    with open(fileSB, "r") as f:
        sbResults = json.load(f)
    with open(fileYields, "r") as f:
        yieldResults = json.load(f)
    with open(fileBinMap, "r") as f:
        binMap = json.load(f)

    # histograms
    h_sgamma_searchbins_1       = ROOT.TH1F("h_sgamma_searchbins_1",        "h_sgamma_searchbins_1",        80, -20, 20)
    h_sgamma_crunits_1          = ROOT.TH1F("h_sgamma_crunits_1",           "h_sgamma_crunits_1",           80, -20, 20)
    h_sgamma_searchbins_2       = ROOT.TH1F("h_sgamma_searchbins_2",        "h_sgamma_searchbins_2",        60, -1, 5)
    h_sgamma_crunits_2          = ROOT.TH1F("h_sgamma_crunits_2",           "h_sgamma_crunits_2",           60, -1, 5)

    # get bin and sgamma values for each search bin
    sgammaSearchBins = list((int(b), sbResults[era][b]["shape"]) for b in sbResults[era])
    sgammaSearchBins.sort(key = lambda x: x[0])
    for x in sgammaSearchBins:
        b       = x[0]
        sgamma  = x[1]
        h_sgamma_searchbins_1.Fill(sgamma)
        h_sgamma_searchbins_2.Fill(sgamma)
        if verbose > 1 and (sgamma < minSgamma or sgamma > maxSgamma):
            print "search bin {0}, sgamma = {1}".format(b, sgamma)
    # get bin and sgamma values for each control region bin
    sgammaCRUnits = []
    CRBinNames = {}
    for binName in binMap["unitCRNum"]["phocr"]:
        b =  int(binMap["unitCRNum"]["phocr"][binName])
        CRBinNames[b] = binName
    
    # keys for yieldResults["yieldsMap"]:
    # znunu
    # phocr_back
    # phocr_data
    # phocr_gjets
    
    for b in xrange(len(CRBinNames)):
        binName = CRBinNames[b]
        phocr_data  = yieldResults["yieldsMap"]["phocr_data"][binName][0]
        phocr_gjets = yieldResults["yieldsMap"]["phocr_gjets"][binName][0]
        phocr_back  = yieldResults["yieldsMap"]["phocr_back"][binName][0]
        sgamma = -999
        den = phocr_gjets + phocr_back
        if den > 0.0:
            sgamma = phocr_data / den
        else:
            print "WARNING: CR bin {0}, denominator = {1}".format(b, den)
        h_sgamma_crunits_1.Fill(sgamma)
        h_sgamma_crunits_2.Fill(sgamma)
        
        if verbose > 1 and (sgamma < minSgamma or sgamma > maxSgamma):
            print "CR bin {0}, sgamma = {1}".format(b, sgamma)
        
        if verbose > 2 and phocr_data < 2.0:
            print "CR bin {0}: phocr_data = {1}".format(b, phocr_data)

    histograms_1 = [h_sgamma_searchbins_1, h_sgamma_crunits_1]
    histograms_2 = [h_sgamma_searchbins_2, h_sgamma_crunits_2]
    labels = ["sgamma_searchbins", "sgamma_crunits"]
    
    # plot
    # plot(histograms, labels, name, title, x_title, y_min, y_max, era)
    plot(histograms_1, labels, "sgamma_binning1", "Sgamma", "Sgamma", 0.0, 200.0, era) 
    plot(histograms_2, labels, "sgamma_binning2", "Sgamma", "Sgamma", 0.0, 60.0,  era) 
    
    del h_sgamma_searchbins_1
    del h_sgamma_searchbins_2
    del h_sgamma_crunits_1
    del h_sgamma_crunits_2

def plot(histograms, labels, name, title, x_title, y_min, y_max, era):
    eraTag = "_" + era
    draw_option = "hist"
            
    # colors
    color_red    = "vermillion"
    color_blue   = "electric blue"
    color_green  = "irish green" 
    color_purple = "violet"
    color_black  = "black"
    colors = [color_red, color_blue, color_green, color_purple, color_black]
    
    # legend: TLegend(x1,y1,x2,y2)
    legend_x1 = 0.6
    legend_x2 = 0.9 
    legend_y1 = 0.7
    legend_y2 = 0.9 
    legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
    
    c = ROOT.TCanvas("c", "c", 800, 800)

    y_title = "Events"
    
    for i in xrange(len(histograms)):
        # setupHist(hist, title, x_title, y_title, color, y_min, y_max, adjust=False)
        setupHist(histograms[i], title, x_title, y_title, colors[i], y_min, y_max, adjust=False)
        # draw
        if i == 0:
            histograms[i].Draw(draw_option)
        else:
            histograms[i].Draw(draw_option + " same")
        legend.AddEntry(histograms[i],    labels[i],   "l")

    legend.Draw()
    
    # save histograms
    plot_dir = "more_plots/"
    plot_name = plot_dir + name + eraTag
    c.Update()
    c.SaveAs(plot_name + ".png")

def main():
    eras = ["2016", "2017", "2018", "Run2"]
    for era in eras:
        run(era)

if __name__ == "__main__":
    main()



