# search_bins.py
import re
import ROOT
import copy
import json
from tools import setupHist, getMultiplicationErrorList, removeCuts, getBinError, ERROR_ZERO, getTexSelection, getTexMultiCut, stringifyMap

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# classes for search bins and validation bins

# Common class is parent class
# ValidationBins and SearchBins inherit from Common
class Common:
    def __init__(self):
        # colors
        self.color_red    = "vermillion"
        self.color_blue   = "electric blue"
        self.color_green  = "irish green" 
        self.color_purple = "violet"
        self.color_black  = "black"
        
        # common atributes
        # WARNING: be careful
        # it is safest to only put constant general attributes here
        self.values = ["norm", "shape", "mc", "pred"]

        self.results_dir = "prediction_histos/"
    
    # ---------------------------------------------------------------------- #
    # setBinValues():                                                        #
    #    - set bin values and errors                                         #
    # ---------------------------------------------------------------------- #
    def setBinValues(self, b_map, h_map, era):
        debug = False
        # Note: h_map should have regions (lowdm, highdm) and h_type (mc, data)
        # h_map[region][h_type]
        
        # Note: bin_i and b are different
        # bin_i is histogram bin number
        # b is validation bin number
        for region in h_map:
            for h_type in h_map[region]:
                bin_i = 1
                for b in b_map[region]:
                    h = h_map[region][h_type]
                    if debug:
                        print "era {0}; b={1} b_i={2} {3} {4}, types {5}, {6}, {7}, {8}".format(era, b, bin_i, region, h_type, type(b), type(bin_i), type(region), type(h_type))
                    value       = h.GetBinContent(bin_i)
                    value_error = h.GetBinError(bin_i)
                    final_error = getBinError(value, value_error, ERROR_ZERO) 
                    self.binValues[era][b][h_type]            = value
                    self.binValues[era][b][h_type + "_error"] = final_error
                    self.binValues[era][b][h_type + "_tex"] = "${0:.3f} \pm {1:.3f}$".format(value, final_error)
                    bin_i += 1

    def makeJson(self, map_, file_):
        with open(file_, "w") as output:
            json.dump(map_, output, sort_keys=True, indent=4, separators=(',', ' : '))

    def writeLine(self, line):
        self.output_file.write(line + "\n")

    def makeTexFile(self, caption, output_name, total_era=""):
        # write latex file with table
        with open(output_name, "w+") as f:
            self.output_file = f
            # begin document
            self.writeLine("\\documentclass{article}")
            self.writeLine("\\usepackage[utf8]{inputenc}")
            self.writeLine("\\usepackage{geometry}")
            self.writeLine("\\usepackage{longtable}")
            self.writeLine("\\usepackage{cancel}")
            self.writeLine("\\geometry{margin=0.1cm}")
            self.writeLine("")
            self.writeLine("\\begin{document}")
            self.writeLine("\\footnotesize")
            self.writeLine("\\tabcolsep=0.01cm")
            # make one table will total Run 2 predictions
            if total_era:
                total_era_tex = total_era.replace("_", " ")
                # begin table
                self.writeLine("\\centering")
                # *n{} syntax with vertical lines for n columns; put last | in expression: *n{...|}
                # make first column for bin numbers small
                self.writeLine("\\begin{longtable}{|p{0.05\\textwidth}|p{0.5\\textwidth}|*2{p{0.2\\textwidth}|}}")
                # column headers
                self.writeLine("\\hline Bin & Selection & $N_{MC}$ & $N_{p}$ \\\\")
                # write values to table
                for b in self.all_bins:
                    total_selection = self.bins[b]["total_selection"]
                    mc              = self.binValues[total_era][b]["mc_tex"]
                    pred            = self.binValues[total_era][b]["pred_tex"]
                    self.writeLine("\\hline {0} & {1} & {2} & {3} \\\\".format(b, total_selection, mc, pred))
                self.writeLine("\\hline")
                # for longtable, caption must go at the bottom of the table... it is not working at the top
                self.writeLine("\\caption{{{0} ({1})}}".format(caption, total_era_tex))
                # end table
                self.writeLine("\\end{longtable}")
            # make a table for each era
            else:
                for era in self.eras:
                    era_tex = era.replace("_", " ")
                    # begin table
                    self.writeLine("\\centering")
                    # *n{} syntax with vertical lines for n columns; put last | in expression: *n{...|}
                    # make first column for bin numbers small
                    self.writeLine("\\begin{longtable}{|p{0.03\\textwidth}|p{0.3\\textwidth}|*6{p{0.1\\textwidth}|}}")
                    # column headers
                    self.writeLine("\\hline Bin & Selection & $R_{Z}$ & $S_{\\gamma}$ & $N_{MC}$ & $N_{p}$ & $\\langle w \\rangle$ & $N_{eff}$ \\\\")
                    # write values to table
                    for b in self.all_bins:
                        total_selection = self.bins[b]["total_selection"]
                        norm            = self.binValues[era][b]["norm_tex"]
                        shape           = self.binValues[era][b]["shape_tex"]
                        mc              = self.binValues[era][b]["mc_tex"]
                        pred            = self.binValues[era][b]["pred_tex"]
                        avg_w           = self.binValues[era][b]["avg_w_tex"]
                        n_eff           = self.binValues[era][b]["n_eff_tex"]
                        avg_w_final     = self.binValues[era][b]["avg_w_final_tex"]
                        n_eff_final     = self.binValues[era][b]["n_eff_final_tex"]
                        self.writeLine("\\hline {0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} \\\\".format(b, total_selection, norm, shape, mc, pred, avg_w, n_eff))
                    self.writeLine("\\hline")
                    # for longtable, caption must go at the bottom of the table... it is not working at the top
                    self.writeLine("\\caption{{{0} ({1})}}".format(caption, era_tex))
                    # end table
                    self.writeLine("\\end{longtable}")
            # end document
            self.writeLine("\\end{document}")

    # ---------------------------------------------------------------------- #
    # makeTotalPred():                                                       #
    #    - add data, mc, and pred for all eras                               #
    #    - plot and save histograms                                          #
    #    - save relevant values to map                                       #
    # ---------------------------------------------------------------------- #
    def makeTotalPred(self, output_file, x_title, name, total_era):
        debug = False
        
        # bin map
        b_map = {}
        b_map["lowdm"]   = self.low_dm_bins
        b_map["highdm"]  = self.high_dm_bins
        
        if debug:
            print "makeTotalPred(): {0} {1} {2} {3}".format(output_file, x_title, name, total_era)
        
        # setup bin values
        self.binValues[total_era] = {}
        for b in self.all_bins:
            self.binValues[total_era][b] = {}
        
        # histogram map
        h_map = {}
        for i, era in enumerate(self.eras):
            # region: lowdm, highdm
            for region in self.histograms[era]:
                if i == 0:
                    h_map[region] = {}
                # h_type: mc, pred, data
                for h_type in self.histograms[era][region]:
                    h = self.histograms[era][region][h_type]
                    if debug:
                        print "DEBUG: {0} {1} {2}, h={3}".format(era, region, h_type, h)
                    if i == 0:
                        h_map[region][h_type] = h.Clone()
                    else:
                        h_map[region][h_type].Add(h)
        
        
        self.histograms[total_era] = h_map
        
        if debug:
            regions = list(r for r in self.histograms[total_era])
            print "regions: {0}".format(regions) 
        
        # set bin values 
        self.setBinValues(b_map, h_map, total_era)
        self.makeHistos(output_file, x_title, name, total_era)

    
    # ---------------------------------------------------------------------- #
    # makeHistos():                                                          #
    #    - make, plot, and save histograms                                   #
    # ---------------------------------------------------------------------- #
    def makeHistos(self, output_file, x_title, name, era):
        debug = False
        eraTag = "_" + era
        draw_option = "hist error"
        if self.saveRootFile:
            #print "Running makeHistos(). Saving output to {0}".format(output_file)
            f_out = ROOT.TFile(output_file, "recreate")
        
        # define histograms 
        if (self.unblind):
            h_data_lowdm    = ROOT.TH1F("data_lowdm",    "data_lowdm",    self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1) 
            h_data_highdm   = ROOT.TH1F("data_highdm",   "data_highdm",   self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1) 
        h_mc_lowdm    = ROOT.TH1F("mc_lowdm",    "mc_lowdm",    self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1) 
        h_mc_highdm   = ROOT.TH1F("mc_highdm",   "mc_highdm",   self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1) 
        h_pred_lowdm  = ROOT.TH1F("pred_lowdm",  "pred_lowdm",  self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1) 
        h_pred_highdm = ROOT.TH1F("pred_highdm", "pred_highdm", self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1) 

        # setup histograms
        #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
        title = "Z to Invisible Data, MC and Prediction for " + era
        if (self.unblind):
            setupHist(h_data_lowdm,    title, x_title, "Events", self.color_black,  10.0 ** -2, 10.0 ** 5)
            setupHist(h_data_highdm,   title, x_title, "Events", self.color_black,  10.0 ** -2, 10.0 ** 5)
        setupHist(h_mc_lowdm,    title, x_title, "Events", self.color_red,  10.0 ** -2, 10.0 ** 5)
        setupHist(h_mc_highdm,   title, x_title, "Events", self.color_red,  10.0 ** -2, 10.0 ** 5)
        setupHist(h_pred_lowdm,  title, x_title, "Events", self.color_blue, 10.0 ** -2, 10.0 ** 5)
        setupHist(h_pred_highdm, title, x_title, "Events", self.color_blue, 10.0 ** -2, 10.0 ** 5)
                
        # set histogram content and error
        bin_i = 1
        for b in self.low_dm_bins:
            if (self.unblind):
                h_data_lowdm.SetBinContent(bin_i,   self.binValues[era][b]["data"])
                h_data_lowdm.SetBinError(bin_i,     self.binValues[era][b]["data_error"])
            h_mc_lowdm.SetBinContent(bin_i,     self.binValues[era][b]["mc"])
            h_mc_lowdm.SetBinError(bin_i,       self.binValues[era][b]["mc_error"])
            h_pred_lowdm.SetBinContent(bin_i,   self.binValues[era][b]["pred"])
            h_pred_lowdm.SetBinError(bin_i,     self.binValues[era][b]["pred_error_mc_only"])
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            if (self.unblind):
                h_data_highdm.SetBinContent(bin_i,  self.binValues[era][b]["data"])
                h_data_highdm.SetBinError(bin_i,    self.binValues[era][b]["data_error"])
            h_mc_highdm.SetBinContent(bin_i,    self.binValues[era][b]["mc"])
            h_mc_highdm.SetBinError(bin_i,      self.binValues[era][b]["mc_error"])
            h_pred_highdm.SetBinContent(bin_i,  self.binValues[era][b]["pred"])
            h_pred_highdm.SetBinError(bin_i,    self.binValues[era][b]["pred_error_mc_only"])
            bin_i += 1

        # histogram map
        h_map = {}
        
        h_map["lowdm"] = {}
        h_map["lowdm"]["mc"]   = h_mc_lowdm
        h_map["lowdm"]["pred"] = h_pred_lowdm
        
        h_map["highdm"] = {}
        h_map["highdm"]["mc"]   = h_mc_highdm
        h_map["highdm"]["pred"] = h_pred_highdm
        
        if (self.unblind):
            h_map["lowdm"]["data"]   = h_data_lowdm
            h_map["highdm"]["data"]  = h_data_highdm
        
        # WARNING
        # - histograms will be deleted when TFile is closed
        # - histograms need to be copied to use them later on 
        self.histograms[era] = copy.deepcopy(h_map)
        
        ###################
        # Draw Histograms #
        ###################
        if self.draw:
            
            # draw histograms
            c = ROOT.TCanvas("c", "c", 800, 800)
            c.Divide(1, 2)
            
            # legend: TLegend(x1,y1,x2,y2)
            legend_x1 = 0.7
            legend_x2 = 0.9 
            legend_y1 = 0.7 
            legend_y2 = 0.9 

            for region in h_map:
                if (self.unblind):
                    h_data   = h_map[region]["data"]
                h_mc   = h_map[region]["mc"]
                h_pred = h_map[region]["pred"]
                h_ratio = h_pred.Clone("h_ratio")
                h_ratio.Divide(h_mc)
            
                #setupHist(hist, title, x_title, y_title, color, y_min, y_max)
                setupHist(h_ratio, "Z to Invisible Prediction / MC", x_title, "Pred / MC", self.color_blue, 0.0, 2.0)

                # histograms
                c.cd(1)
                ROOT.gPad.SetLogy(1) # set log y
                # ZInv MC and Prediction
                h_mc.Draw(draw_option)
                h_pred.Draw("error same")
                if (self.unblind):
                    h_data.Draw("error same")
                
                # legend: TLegend(x1,y1,x2,y2)
                legend = ROOT.TLegend(legend_x1, legend_y1, legend_x2, legend_y2)
                if (self.unblind):
                    legend.AddEntry(h_data,   "MET Data",   "l")
                legend.AddEntry(h_mc,   "Z#rightarrow#nu#nu MC",   "l")
                legend.AddEntry(h_pred, "Z#rightarrow#nu#nu Pred", "l")
                legend.Draw()
               
                # ratios
                c.cd(2)
                h_ratio.Draw(draw_option)
                    
                # save histograms
                plot_name = self.plot_dir + name + "_" + region + eraTag
                c.Update()
                c.SaveAs(plot_name + ".pdf")
                c.SaveAs(plot_name + ".png")
                # debug
                if debug:
                    print "region: {0}".format(region)
                    if (self.unblind):
                        print "h_data[1] = {0}".format(h_data.GetBinContent(1))
                    print "h_mc[1] = {0}".format(h_mc.GetBinContent(1))
                    print "h_pred[1] = {0}".format(h_pred.GetBinContent(1))
                
                # write histograms to file
                if self.saveRootFile:
                    if (self.unblind):
                        h_data.Write()
                    h_mc.Write()
                    h_pred.Write()
        
        if self.saveRootFile:
            f_out.Close()
    
    # ---------------------------------------------------------------------- #
    # calcPrediction():                                                      #
    #    - calculate final prediction                                        #
    #    - save relevant values to map                                       #
    # ---------------------------------------------------------------------- #
    def calcPrediction(self, output_file, x_title, name, era):
        # save values in map
        if self.verbose:
            print era
        for b in self.all_bins:
            region          = self.bins[b]["region"]
            selection       = self.bins[b]["selection"]
            met             = self.bins[b]["met"]
            selection_tex   = getTexSelection(region + "_" + selection)
            met_tex         = getTexMultiCut(met)
            total_selection = "{0}, {1}".format(selection_tex, met_tex)
            self.bins[b]["total_selection"] = total_selection
            
            n       = self.binValues[era][b]["norm"]
            n_error = self.binValues[era][b]["norm_error"]
            s       = self.binValues[era][b]["shape"]
            s_error = self.binValues[era][b]["shape_error"]
            m       = self.binValues[era][b]["mc"]
            m_error = self.binValues[era][b]["mc_error"]
            
            # prediction:                   p     = bin value
            # uncertainty:                  sigma = bin error
            # average weight:               avg_w = sigma^2 / p 
            # effective number of events:   N_eff = p / avg_w
            
            # For data card and validation bin histograms:
            # - do not propagate Rz statistical error because total Rz uncertainty will be included in Rz systematic histogram
            # - do not propagate Sgamma statistical error because Higgs Combine will do shape factor for us
            p_error_mc_only    = m_error
            # For prediction table in analysis note:
            # - propagate Rz and Sg errors
            p                  = n * s * m
            x_list             = [n, s, m]
            dx_list            = [n_error, s_error, m_error]
            p_error_propagated = getMultiplicationErrorList(p, x_list, dx_list)
            
            # error < 0.0 due to error code
            if p_error_propagated < 0.0:
                p_error_propagated = ERROR_ZERO 
            if p == 0:
                if self.verbose:
                    print "WARNING: bin {0}, pred = {1}; seting avg weight to {2}".format(b, p, ERROR_ZERO)
                avg_w   = ERROR_ZERO
            else:
                avg_w   = (p_error_propagated ** 2) / p
            n_eff = p / avg_w
            n_eff_final = int(n_eff)
            if n_eff_final == 0:
                if self.verbose:
                    print "WARNING: bin {0}, n_eff_final = {1}; leaving avg weight unchanged".format(b, n_eff_final)
                avg_w_final = avg_w
            else:
                avg_w_final = p / n_eff_final

            self.binValues[era][b]["pred"]                  = p
            self.binValues[era][b]["pred_error_mc_only"]    = p_error_mc_only
            self.binValues[era][b]["pred_error_propagated"] = p_error_propagated
            self.binValues[era][b]["avg_w"]                 = avg_w
            self.binValues[era][b]["n_eff"]                 = n_eff
            self.binValues[era][b]["avg_w_final"]           = avg_w_final
            self.binValues[era][b]["n_eff_final"]           = n_eff_final
            
            for value in self.values:
                error_suffix = "_error"
                if value == "pred":
                    error_suffix = "_error_propagated"
                self.binValues[era][b][value + "_tex"] = "${0:.3f} \pm {1:.3f}$".format(self.binValues[era][b][value], self.binValues[era][b][value + error_suffix])
                
            for value in ["avg_w", "n_eff", "avg_w_final", "n_eff_final"]:
                self.binValues[era][b][value + "_tex"] = "${0:.3f}$".format(self.binValues[era][b][value])

            if self.verbose:
                print "bin {0}: N = {1:.3f} +/- {2:.3f} S = {3:.3f} +/- {4:.3f} M = {5:.3f} +/- {6:.3f} P = {7:.3f} +/- {8:.3f}".format(
                            b, n, n_error, s, s_error, m, m_error, p, p_error_propagated 
                        )

    # ---------------------------------------------------------------------- #
    # getRzSyst():                                                           #
    #    - get Rz systematic in validation and search bins                   #
    #    - save systematic to root file                                      #
    # ---------------------------------------------------------------------- #
    def getRzSyst(self, rz_syst_map, bin_type, output_file):
        # output root file
        f_out = ROOT.TFile(output_file, "recreate")
        
        rz_syst_low_dm  = ROOT.TH1F("rz_syst_low_dm",  "rz_syst_low_dm",  self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1)
        rz_syst_high_dm = ROOT.TH1F("rz_syst_high_dm", "rz_syst_high_dm", self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1)
        
        regions = ["LowDM", "HighDM"]
        systMap = {}
        systMap["LowDM"]  = {}
        systMap["HighDM"] = {}
        systMap["LowDM"]["bins"]   = self.low_dm_bins
        systMap["HighDM"]["bins"]  = self.high_dm_bins
        systMap["LowDM"]["hist"]   = rz_syst_low_dm
        systMap["HighDM"]["hist"]  = rz_syst_high_dm
       
        for region in regions:
            i = 1
            for b in systMap[region]["bins"]:
                selection = self.bins[b]["selection"]
                selection_norm = removeCuts(selection, ["NJ"])
                
                #print "DEBUG: b={0}, selection={1}, selection_norm={2}, type={3}".format(b, selection, selection_norm, type(selection_norm))
                #print "DEBUG: rz_syst_map[{0}] keys = {1}".format(region, rz_syst_map[bin_type][region].keys())
                
                # rz_syst_map[bin_type][region][selection]
                syst = rz_syst_map[bin_type][region][selection_norm]
                
                #print "--- syst = {0}".format(syst)
                systMap[region]["hist"].SetBinContent(i, syst)
                systMap[region]["hist"].SetBinError(i, 0)
                i += 1
            # write histo to file
            systMap[region]["hist"].Write()
        
        f_out.Close()
    
    # ---------------------------------------------------------------------- #
    # getZvsPhotonSyst():                                                    #
    #    - get Z vs Photon systematic in validation, search, or CR unit bins #
    #    - save systematic to root file                                      #
    # ---------------------------------------------------------------------- #
    def getZvsPhotonSyst(self, h_map_syst, output_file):
        # output root file
        f_out = ROOT.TFile(output_file, "recreate")
        ZvsPhoton_syst_low_dm  = ROOT.TH1F("ZvsPhoton_syst_low_dm",  "ZvsPhoton_syst_low_dm",  self.low_dm_nbins,  self.low_dm_start,  self.low_dm_end + 1)
        ZvsPhoton_syst_high_dm = ROOT.TH1F("ZvsPhoton_syst_high_dm", "ZvsPhoton_syst_high_dm", self.high_dm_nbins, self.high_dm_start, self.high_dm_end + 1)
        
        regions = ["LowDM", "HighDM"]
        systMap = {}
        systMap["LowDM"]  = {}
        systMap["HighDM"] = {}
        systMap["LowDM"]["bins"]   = self.low_dm_bins
        systMap["HighDM"]["bins"]  = self.high_dm_bins
        systMap["LowDM"]["hist"]   = ZvsPhoton_syst_low_dm
        systMap["HighDM"]["hist"]  = ZvsPhoton_syst_high_dm
        
        for region in regions:
            i = 1
            for b in systMap[region]["bins"]:
                # get MET bin edges
                met_cut = self.bins[b]["met"]
                m = re.search("met_(.*)to(.*)", met_cut)
                value_1 = float(m.group(1))
                value_2 = float(m.group(2))
                # use same systematic for validation and search bins
                # find correct MET bin to use
                # h_map_syst[region]
                met_bin_num = h_map_syst[region].FindBin(value_1)
                syst = h_map_syst[region].GetBinContent(met_bin_num)
                
                systMap[region]["hist"].SetBinContent(i, syst)
                systMap[region]["hist"].SetBinError(i, 0)
                i += 1
            # write histo to file
            systMap[region]["hist"].Write()
        
        f_out.Close()

# vadliation bins
class ValidationBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose, draw, saveRootFile):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.draw = draw
        self.saveRootFile = saveRootFile
        self.binValues = {}
        self.histograms = {}
        self.unblind = True
        # SBv3
        self.low_dm_start           = 0
        self.low_dm_normal_end      = 14
        self.low_dm_highmet_start   = 15
        self.low_dm_end             = 18
        self.high_dm_start          = 19
        self.high_dm_end            = 42
        self.low_dm_nbins           = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins          = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins            = list(str(b) for b in range( self.low_dm_start,         self.low_dm_end + 1)) 
        self.high_dm_bins           = list(str(b) for b in range( self.high_dm_start,        self.high_dm_end + 1)) 
        self.low_dm_bins_normal     = list(str(b) for b in range( self.low_dm_start,         self.low_dm_normal_end + 1)) 
        self.low_dm_bins_highmet    = list(str(b) for b in range( self.low_dm_highmet_start, self.low_dm_end + 1))
        self.all_bins               = self.low_dm_bins + self.high_dm_bins
        with open("validation_bins_v3.json", "r") as j:
            self.bins = stringifyMap(json.load(j))

    def getValues(self, file_name, era, systTag=""):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            # remove cuts from selection for norm and shape
            selection_norm  = removeCuts(selection, ["NJ"])
            selection_shape = removeCuts(selection, ["NSV", "MET"])
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]                  = self.N.norm_map[era]["validation"]["Combined"][region][selection_norm]["R_Z"]
            self.binValues[era][b]["norm_error"]            = self.N.norm_map[era]["validation"]["Combined"][region][selection_norm]["R_Z_error"]
            self.binValues[era][b]["shape"]                 = self.S.shape_map[era]["validation"][region][selection_shape][met]
            self.binValues[era][b]["shape_error"]           = self.S.shape_map[era]["validation"][region][selection_shape][met + "_error"]
            self.binValues[era][b]["photon_data_mc_norm"]   = self.S.shape_map[era]["validation"][region][selection_shape]["photon_data_mc_norm"]

        # Z to NuNu histograms
        # central value
        # TDirectoryFile*     nValidationBinLowDM_jetpt30 nValidationBinLowDM_jetpt30
        # KEY: TH1D MET_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30Data MET Validation Bin Low DMdata;1  nValidationBinLowDM_jetpt30
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1 nValidationBinLowDM_jetpt30
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_njetWeight_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1  nValidationBinLowDM_jetpt30
        # systematics
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_jes_syst_down_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1   nValidationBinLowDM_jetpt30
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_jes_syst_up_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1 nValidationBinLowDM_jetpt30
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_btag_syst_down_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1  nValidationBinLowDM_jetpt30
        # KEY: TH1D ZNuNu_nValidationBin_LowDM_btag_syst_up_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata;1    nValidationBinLowDM_jetpt30

        # WARNING: only apply systematics to MC, not Data
        
        f_in                   = ROOT.TFile(file_name, "read")
        if (self.unblind):
            h_data_lowdm           = f_in.Get("nValidationBinLowDM_jetpt30/MET_nValidationBin_LowDM_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30Data MET Validation Bin Low DMdata") 
            h_data_lowdm_highmet   = f_in.Get("nValidationBinLowDMHighMET_jetpt30/MET_nValidationBin_LowDM_HighMET_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30Data MET Validation Bin Low DM High METdata")
            h_data_highdm          = f_in.Get("nValidationBinHighDM_jetpt30/MET_nValidationBin_HighDM_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30Data MET Validation Bin High DMdata")
        h_mc_lowdm             = f_in.Get("nValidationBinLowDM_jetpt30/ZNuNu_nValidationBin_LowDM"                + systTag + "_jetpt30nValidationBinLowDM_jetpt30nValidationBinLowDM_jetpt30ZJetsToNuNu Validation Bin Low DMdata")
        h_mc_lowdm_highmet     = f_in.Get("nValidationBinLowDMHighMET_jetpt30/ZNuNu_nValidationBin_LowDM_HighMET" + systTag + "_jetpt30nValidationBinLowDMHighMET_jetpt30nValidationBinLowDMHighMET_jetpt30ZJetsToNuNu Validation Bin Low DM High METdata")
        h_mc_highdm            = f_in.Get("nValidationBinHighDM_jetpt30/ZNuNu_nValidationBin_HighDM"              + systTag + "_jetpt30nValidationBinHighDM_jetpt30nValidationBinHighDM_jetpt30ZJetsToNuNu Validation Bin High DMdata")
        
        # bin map
        b_map = {}
        b_map["lowdm"]          = self.low_dm_bins_normal
        b_map["lowdm_highmet"]  = self.low_dm_bins_highmet
        b_map["highdm"]         = self.high_dm_bins
        # histogram map
        h_map                           = {}
        h_map["lowdm"]                  = {}
        h_map["lowdm_highmet"]          = {}
        h_map["highdm"]                 = {}
        h_map["lowdm"]["mc"]            = h_mc_lowdm
        h_map["lowdm_highmet"]["mc"]    = h_mc_lowdm_highmet
        h_map["highdm"]["mc"]           = h_mc_highdm
        if (self.unblind):
            h_map["lowdm"]["data"]          = h_data_lowdm
            h_map["lowdm_highmet"]["data"]  = h_data_lowdm_highmet
            h_map["highdm"]["data"]         = h_data_highdm
        
        # set bin values 
        self.setBinValues(b_map, h_map, era)

        # new root file to save validation bin histograms
        new_file = self.results_dir + "validationBinsZinv_" + era + ".root"
        self.calcPrediction(    new_file, "Validation Bin", "validation", era   )
        self.makeHistos(        new_file, "Validation Bin", "validation", era   )
        f_in.Close()

# vadliation bins for MET study
class ValidationBinsMETStudy(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose, draw, saveRootFile):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.draw = draw
        self.saveRootFile = saveRootFile
        self.binValues = {}
        self.histograms = {}
        self.unblind = True
        self.low_dm_start  = 0
        self.low_dm_end    = 3
        self.high_dm_start = 4
        self.high_dm_end   = 8
        self.low_dm_nbins  = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins   = list(str(b) for b in range( self.low_dm_start,   self.low_dm_end + 1)) 
        self.high_dm_bins  = list(str(b) for b in range( self.high_dm_start,  self.high_dm_end + 1)) 
        self.all_bins      = self.low_dm_bins + self.high_dm_bins
        with open("validation_bins_metStudy.json", "r") as j:
            self.bins = stringifyMap(json.load(j))

    def getValues(self, file_name, era, systTag=""):
        doIntegral = False
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            # remove cuts from selection for norm and shape
            selection_norm  = removeCuts(selection, ["NJ"])
            selection_shape = removeCuts(selection, ["NSV", "MET"])
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]                  = self.N.norm_map[era]["validationMetStudy"]["Combined"][region][selection_norm]["R_Z"]
            self.binValues[era][b]["norm_error"]            = self.N.norm_map[era]["validationMetStudy"]["Combined"][region][selection_norm]["R_Z_error"]
            self.binValues[era][b]["shape"]                 = self.S.shape_map[era]["validationMetStudy"][region][selection_shape][met]
            self.binValues[era][b]["shape_error"]           = self.S.shape_map[era]["validationMetStudy"][region][selection_shape][met + "_error"]
            self.binValues[era][b]["photon_data_mc_norm"]   = self.S.shape_map[era]["validationMetStudy"][region][selection_shape]["photon_data_mc_norm"]

        # --- new (for MET study) --- #
        # "nValidationBinLowDM_METStudy_jetpt30/MET_nValidationBin_LowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30Data MET Validation Bin Low DM MET Studydata"
        # "nValidationBinHighDM_METStudy_jetpt30/MET_nValidationBin_HighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30Data MET Validation Bin High DM MET Studydata"
        # "nValidationBinLowDM_METStudy_jetpt30/ZNuNu_nValidationBin_LowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30ZJetsToNuNu Validation Bin Low DM MET Studydata"
        # "nValidationBinHighDM_METStudy_jetpt30/ZNuNu_nValidationBin_HighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30ZJetsToNuNu Validation Bin High DM MET Studydata"

        # WARNING: only apply systematics to MC, not Data
        
        f_in                   = ROOT.TFile(file_name, "read")
        if (self.unblind):
            h_data_lowdm           = f_in.Get("nValidationBinLowDM_METStudy_jetpt30/MET_nValidationBin_LowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30Data MET Validation Bin Low DM MET Studydata") 
            h_data_highdm          = f_in.Get("nValidationBinHighDM_METStudy_jetpt30/MET_nValidationBin_HighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30Data MET Validation Bin High DM MET Studydata")
        h_mc_lowdm             = f_in.Get("nValidationBinLowDM_METStudy_jetpt30/ZNuNu_nValidationBin_LowDM_METStudy"   + systTag + "_jetpt30nValidationBinLowDM_METStudy_jetpt30nValidationBinLowDM_METStudy_jetpt30ZJetsToNuNu Validation Bin Low DM MET Studydata")
        h_mc_highdm            = f_in.Get("nValidationBinHighDM_METStudy_jetpt30/ZNuNu_nValidationBin_HighDM_METStudy" + systTag + "_jetpt30nValidationBinHighDM_METStudy_jetpt30nValidationBinHighDM_METStudy_jetpt30ZJetsToNuNu Validation Bin High DM MET Studydata")

        if doIntegral:
            # Integrate to count number of events. 
            nDataLowDM  = h_data_lowdm.Integral(1, h_data_lowdm.GetNbinsX()) 
            nDataHighDM = h_data_highdm.Integral(1, h_data_highdm.GetNbinsX())
            nMCLowDM    = h_mc_lowdm.Integral(1, h_mc_lowdm.GetNbinsX())
            nMCHighDM   = h_mc_highdm.Integral(1, h_mc_highdm.GetNbinsX())
            print " --- Era: {0} --- ".format(era)
            print "    nDataLowDM = {0}, nDataHighDM = {1}".format(nDataLowDM, nDataHighDM)
            print "    nMCLowDM = {0},   nMCHighDM = {1}".format(nMCLowDM, nMCHighDM)
        
        # bin map
        b_map = {}
        b_map["lowdm"]          = self.low_dm_bins
        b_map["highdm"]         = self.high_dm_bins
        # histogram map
        h_map                           = {}
        h_map["lowdm"]                  = {}
        h_map["highdm"]                 = {}
        h_map["lowdm"]["mc"]            = h_mc_lowdm
        h_map["highdm"]["mc"]           = h_mc_highdm
        if (self.unblind):
            h_map["lowdm"]["data"]          = h_data_lowdm
            h_map["highdm"]["data"]         = h_data_highdm
        
        # set bin values 
        self.setBinValues(b_map, h_map, era)

        # new root file to save validation bin histograms
        new_file = self.results_dir + "validationBinsMETStudyZinv_" + era + ".root"
        self.calcPrediction(    new_file, "Validation Bin", "validationMetStudy", era   )
        self.makeHistos(        new_file, "Validation Bin", "validationMetStudy", era   )
        f_in.Close()
  
# search bins 
class SearchBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose, draw, saveRootFile):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        self.unblind = False
        self.draw = draw
        self.saveRootFile = saveRootFile
        # SBv4
        self.low_dm_start   = 0
        self.low_dm_end     = 52
        self.high_dm_start  = 53
        self.high_dm_end    = 182
        self.low_dm_nbins   = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins  = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins    = list(str(b) for b in range( self.low_dm_start,  self.low_dm_end + 1)) 
        self.high_dm_bins   = list(str(b) for b in range( self.high_dm_start, self.high_dm_end + 1)) 
        self.all_bins       = self.low_dm_bins + self.high_dm_bins
        self.binValues      = {}
        self.histograms     = {}
        with open("search_bins_v4.json", "r") as j:
            self.bins = stringifyMap(json.load(j))
        with open("units.json", "r") as input_file:
            self.unitMap = json.load(input_file)
    
    def getValues(self, file_name, era, systTag="", CRunits=0):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            region      = self.bins[b]["region"]
            selection   = self.bins[b]["selection"]
            met         = self.bins[b]["met"]
            # remove cuts from selection for norm and shape
            selection_norm  = removeCuts(selection, ["NJ"])
            selection_shape = removeCuts(selection, ["NSV", "MET"])
            if self.verbose:
                print "{0}: {1} {2} {3} {4}".format(b, region, selection_norm, selection_shape, met)
            self.binValues[era][b] = {}
            self.binValues[era][b]["norm"]                  = self.N.norm_map[era]["search"]["Combined"][region][selection_norm]["R_Z"]
            self.binValues[era][b]["norm_error"]            = self.N.norm_map[era]["search"]["Combined"][region][selection_norm]["R_Z_error"]
            photon_data_mc_norm                             = self.S.shape_map[era]["search"][region][selection_shape]["photon_data_mc_norm"]
            self.binValues[era][b]["photon_data_mc_norm"]   = photon_data_mc_norm
            
            if CRunits:
                # get SR and CR unit bins for this search bin
                cr_units = self.unitMap["unitBinMapCR_phocr"][b]
                # loop over units for this search bin
                total_data = 0
                total_mc   = 0
                for cr in cr_units:
                    data                = CRunits.binValues[era][cr]["data"]
                    data_error          = CRunits.binValues[era][cr]["data_error"]
                    mc_gjets            = CRunits.binValues[era][cr]["mc_gjets"]
                    mc_gjets_error      = CRunits.binValues[era][cr]["mc_gjets_error"]
                    mc_back             = CRunits.binValues[era][cr]["mc_back"]
                    mc_back_error       = CRunits.binValues[era][cr]["mc_back_error"]
                    total_data += data
                    total_mc   += mc_gjets
                    total_mc   += mc_back
                # check for 0 data
                if total_data <= 0:
                    print "WARNING: Search bin {0}: NO DATA: total_data = {1} and total_mc = {2}".format(b, total_data, total_mc)
                # avoid dividing by 0
                shape_cr = -999
                if photon_data_mc_norm * total_mc:
                    # S = sum(data) / (Q * sum(MC))
                    shape_cr = total_data / (photon_data_mc_norm * total_mc)
                else:
                    print "WARNING: Search bin {0}: NO MC: total_data = {1} and total_mc = {2}".format(b, total_data, total_mc)
                self.binValues[era][b]["shape"]                 = shape_cr 
                self.binValues[era][b]["shape_error"]           = -999 
            else:
                self.binValues[era][b]["shape"]                 = self.S.shape_map[era]["search"][region][selection_shape][met]
                self.binValues[era][b]["shape_error"]           = self.S.shape_map[era]["search"][region][selection_shape][met + "_error"]
        
        # Z to NuNu MC histograms
        f_in            = ROOT.TFile(file_name, "read")
        if (self.unblind):
            h_data_lowdm    = f_in.Get("nSearchBinLowDM_jetpt30/MET_nSearchBin_LowDM_jetpt30nSearchBinLowDM_jetpt30nSearchBinLowDM_jetpt30Data MET Search Bin Low DMdata") 
            h_data_highdm   = f_in.Get("nSearchBinHighDM_jetpt30/MET_nSearchBin_HighDM_jetpt30nSearchBinHighDM_jetpt30nSearchBinHighDM_jetpt30Data MET Search Bin High DMdata")
        h_mc_lowdm      = f_in.Get("nSearchBinLowDM_jetpt30/ZNuNu_nSearchBin_LowDM"    + systTag + "_jetpt30nSearchBinLowDM_jetpt30nSearchBinLowDM_jetpt30ZJetsToNuNu Search Bin Low DMdata")
        h_mc_highdm     = f_in.Get("nSearchBinHighDM_jetpt30/ZNuNu_nSearchBin_HighDM"  + systTag + "_jetpt30nSearchBinHighDM_jetpt30nSearchBinHighDM_jetpt30ZJetsToNuNu Search Bin High DMdata")
        
        # bin map
        b_map = {}
        b_map["lowdm"]   = self.low_dm_bins
        b_map["highdm"]  = self.high_dm_bins
        # histogram map
        h_map                 = {}
        h_map["lowdm"]        = {}
        h_map["highdm"]       = {}
        h_map["lowdm"]["mc"]  = h_mc_lowdm
        h_map["highdm"]["mc"] = h_mc_highdm
        if (self.unblind):
            h_map["lowdm"]["data"]  = h_data_lowdm
            h_map["highdm"]["data"] = h_data_highdm
        
        # set bin values 
        self.setBinValues(b_map, h_map, era)

        # new root file to save search bin histograms
        new_file = self.results_dir + "searchBinsZinv_" + era + ".root"
        self.calcPrediction(    new_file, "Search Bin", "search", era   )
        self.makeHistos(        new_file, "Search Bin", "search", era   )
        f_in.Close()

# search region unit bins
class SRUnitBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        # SBv4
        self.low_dm_start   = 0
        self.low_dm_end     = 52
        self.high_dm_start  = 53
        self.high_dm_end    = 528
        self.low_dm_nbins   = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins  = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins    = list(str(b) for b in range( self.low_dm_start,  self.low_dm_end + 1)) 
        self.high_dm_bins   = list(str(b) for b in range( self.high_dm_start, self.high_dm_end + 1)) 
        self.all_bins       = self.low_dm_bins + self.high_dm_bins
        self.binValues      = {}
        self.histograms     = {}
    
    # Z to LL normalization is applied to SR unit bins in units.py
    def getValues(self, file_name, era):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            self.binValues[era][b] = {}

        # Z to NuNu MC histograms
        f_in            = ROOT.TFile(file_name, "read")
        h_mc_lowdm      = f_in.Get("nSRUnitLowDM_jetpt30/ZNuNu_nSRUnit_LowDM_jetpt30nSRUnitLowDM_jetpt30nSRUnitLowDM_jetpt30ZJetsToNuNu Search Region Unit Low DMdata")  
        h_mc_highdm     = f_in.Get("nSRUnitHighDM_jetpt30/ZNuNu_nSRUnit_HighDM_jetpt30nSRUnitHighDM_jetpt30nSRUnitHighDM_jetpt30ZJetsToNuNu Search Region Unit High DMdata")
         
        # bin map
        b_map                       = {}
        b_map["lowdm"]              = self.low_dm_bins
        b_map["highdm"]             = self.high_dm_bins
        # histogram map
        h_map                       = {}
        h_map["lowdm"]              = {}
        h_map["highdm"]             = {}
        h_map["lowdm"]["mc"]        = h_mc_lowdm
        h_map["highdm"]["mc"]       = h_mc_highdm
        
        # set bin values 
        self.setBinValues(b_map, h_map, era)

        f_in.Close()

# control region unit bins
class CRUnitBins(Common):
    def __init__(self, normalization, shape, eras, plot_dir, verbose):
        # run parent init function
        Common.__init__(self)
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.plot_dir = plot_dir
        self.verbose = verbose
        # SBv4
        self.low_dm_start   = 0
        self.low_dm_end     = 52
        self.high_dm_start  = 53
        self.high_dm_end    = 111 
        self.low_dm_nbins   = self.low_dm_end - self.low_dm_start + 1 
        self.high_dm_nbins  = self.high_dm_end - self.high_dm_start + 1 
        self.low_dm_bins    = list(str(b) for b in range( self.low_dm_start,  self.low_dm_end + 1)) 
        self.high_dm_bins   = list(str(b) for b in range( self.high_dm_start, self.high_dm_end + 1)) 
        self.all_bins       = self.low_dm_bins + self.high_dm_bins
        self.binValues      = {}
        self.histograms     = {}
        with open("control_region_unit_bins_v4.json", "r") as j:
            self.bins = stringifyMap(json.load(j))
    
    # TODO: normalize photon MC to Data in CR unit bins
    def getValues(self, file_name, era):
        self.binValues[era] = {}
        
        for b in self.all_bins:
            self.binValues[era][b] = {}

        h_data_lowdm        = self.S.cr_unit_histos_summed[era]["LowDM"]["data"]
        h_data_highdm       = self.S.cr_unit_histos_summed[era]["HighDM"]["data"] 
        h_mc_back_lowdm     = self.S.cr_unit_histos_summed[era]["LowDM"]["mc_back"] 
        h_mc_back_highdm    = self.S.cr_unit_histos_summed[era]["HighDM"]["mc_back"] 
        h_mc_gjets_lowdm    = self.S.cr_unit_histos_summed[era]["LowDM"]["mc_gjets"] 
        h_mc_gjets_highdm   = self.S.cr_unit_histos_summed[era]["HighDM"]["mc_gjets"] 
        
        # bin map
        b_map                           = {}
        b_map["lowdm"]                  = self.low_dm_bins
        b_map["highdm"]                 = self.high_dm_bins
        # histogram map
        h_map                           = {}
        h_map["lowdm"]                  = {}
        h_map["highdm"]                 = {}
        h_map["lowdm"]["data"]          = h_data_lowdm
        h_map["highdm"]["data"]         = h_data_highdm
        h_map["lowdm"]["mc_gjets"]      = h_mc_gjets_lowdm
        h_map["highdm"]["mc_gjets"]     = h_mc_gjets_highdm
        h_map["lowdm"]["mc_back"]       = h_mc_back_lowdm
        h_map["highdm"]["mc_back"]      = h_mc_back_highdm
        
        # A temporary solution from Angel to Caleb
        self.histograms[era] = copy.deepcopy(h_map)
        
        # set bin values 
        self.setBinValues(b_map, h_map, era)


