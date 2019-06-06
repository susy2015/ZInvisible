# searchBins.py

import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# search bins and validation bins

# search bins 
class SearchBins:
    def __init__(self, normalization, shape, eras, verbose):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.verbose = verbose
        self.low_dm_bins  = range( 0, 1)
        self.high_dm_bins = range( 1, 2)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.binValues = {}
        self.eras = []

# vadliation bins
class ValidationBins:
    def __init__(self, normalization, shape, eras, verbose):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.verbose = verbose
        self.binValues = {}
        # bins 0 to 45; bins 19, 20, and 21 are not included
        self.low_dm_bins_normal     = range( 0, 15)
        self.low_dm_bins_highmet    = range(15, 19)
        self.low_dm_bins            = range( 0, 19)
        self.high_dm_bins           = range(22, 46)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.low_dm_bins_met_250to300  = [4, 5, 8, 9, 10, 11]
        self.low_dm_bins_met_300toINF  = [17, 18]
        self.low_dm_bins_met_250to400  = [0, 1, 2, 3, 6, 7, 12, 13, 14]
        self.low_dm_bins_met_400toINF  = [15, 16]
        self.high_dm_bins_met_250to400 = range(22, 46, 2)
        self.high_dm_bins_met_400toINF = range(23, 46, 2)

    def getValues(self, file_name, era):
        year = era[0:4]
        self.binValues[era] = {}
        for b in self.all_bins:
            self.binValues[era][b] = {}
        # normalization
        for b in self.low_dm_bins:
            self.binValues[era][b]["norm"]       = self.N.norm_map[era]["Combined"]["LowDM"]["R_Z"]
            self.binValues[era][b]["norm_error"] = self.N.norm_map[era]["Combined"]["LowDM"]["R_Z_error"]
        for b in self.high_dm_bins:
            self.binValues[era][b]["norm"]       = self.N.norm_map[era]["Combined"]["HighDM"]["R_Z"]
            self.binValues[era][b]["norm_error"] = self.N.norm_map[era]["Combined"]["HighDM"]["R_Z_error"]
        # shape
        # LowDM
        for b in self.low_dm_bins_met_250to300: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["LowDM"]["met_250to300"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["LowDM"]["met_250to300_error"]
        for b in self.low_dm_bins_met_300toINF: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["LowDM"]["met_300toINF"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["LowDM"]["met_300toINF_error"]
        for b in self.low_dm_bins_met_250to400: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["LowDM"]["met_250to400"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["LowDM"]["met_250to400_error"]
        for b in self.low_dm_bins_met_400toINF: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["LowDM"]["met_400toINF"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["LowDM"]["met_400toINF_error"]
        # HighDM
        for b in self.high_dm_bins_met_250to400: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["HighDM"]["met_250to400"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["HighDM"]["met_250to400_error"]
        for b in self.high_dm_bins_met_400toINF: 
            self.binValues[era][b]["shape"]       = self.S.shape_map[era]["HighDM"]["met_400toINF"]
            self.binValues[era][b]["shape_error"] = self.S.shape_map[era]["HighDM"]["met_400toINF_error"]

        # Z to NuNu MC values 
        # nValidationBinLowDM: ZNuNu_nValidationBinLowDM_2016nValidationBinLowDMnValidationBinLowDMZJetsToNuNu Validation Bin Low DMdata
        # nValidationBinLowDMHighMET: ZNuNu_nValidationBinLowDMHighMET_2016nValidationBinLowDMHighMETnValidationBinLowDMHighMETZJetsToNuNu Validation Bin Low DM High METdata
        # nValidationBinHighDM: ZNuNu_nValidationBinHighDM_2016nValidationBinHighDMnValidationBinHighDMZJetsToNuNu Validation Bin High DMdata
        f = ROOT.TFile(file_name)
        h_validation_lowdm          = f.Get("nValidationBinLowDM/ZNuNu_nValidationBinLowDM_" + year + "nValidationBinLowDMnValidationBinLowDMZJetsToNuNu Validation Bin Low DMdata")
        h_validation_lowdm_highmet  = f.Get("nValidationBinLowDMHighMET/ZNuNu_nValidationBinLowDMHighMET_" + year + "nValidationBinLowDMHighMETnValidationBinLowDMHighMETZJetsToNuNu Validation Bin Low DM High METdata")
        h_validation_highdm         = f.Get("nValidationBinHighDM/ZNuNu_nValidationBinHighDM_" + year + "nValidationBinHighDMnValidationBinHighDMZJetsToNuNu Validation Bin High DMdata")
        
        # Note: bin_i and b are different
        # bin_i is histogram bin number
        # b is validation bin number
        bin_i = 1
        for b in self.low_dm_bins_normal:
            self.binValues[era][b]["mc"]       = h_validation_lowdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_lowdm.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.low_dm_bins_highmet:
            self.binValues[era][b]["mc"]       = h_validation_lowdm_highmet.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_lowdm_highmet.GetBinError(bin_i)
            bin_i += 1
        bin_i = 1
        for b in self.high_dm_bins:
            self.binValues[era][b]["mc"]       = h_validation_highdm.GetBinContent(bin_i)
            self.binValues[era][b]["mc_error"] = h_validation_highdm.GetBinError(bin_i)
            bin_i += 1

        # print values 
        if self.verbose:
            print era
        for b in self.all_bins:
            n       = self.binValues[era][b]["norm"]
            n_error = self.binValues[era][b]["norm_error"]
            s       = self.binValues[era][b]["shape"]
            s_error = self.binValues[era][b]["shape_error"]
            m       = self.binValues[era][b]["mc"]
            m_error = self.binValues[era][b]["mc_error"]
            p       = n * s * m
            x_list = [n, s, m]
            dx_list = [n_error, s_error, m_error]
            p_error = self.N.getMultiplicationErrorList(p, x_list, dx_list)
            if self.verbose:
                print "bin {0}: N = {1:.3f} +/- {2:.3f} S = {3:.3f} +/- {4:.3f} M = {5:.3f} +/- {6:.3f} P = {7:.3f} +/- {8:.3f}".format(
                            b, n, n_error, s, s_error, m, m_error, p, p_error 
                        )


