# searchBins.py

# search bins and validation bins

# search bins 
class SearchBins:
    def __init__(self, normalization, shape, eras):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.low_dm_bins  = range( 0, 1)
        self.high_dm_bins = range( 1, 2)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.binValues = {}
        self.eras = []

# vadliation bins
class ValidationBins:
    def __init__(self, normalization, shape, eras):
        self.N = normalization
        self.S = shape
        self.eras = eras
        self.binValues = {}
        # bins 0 to 45; bins 19, 20, and 21 are not included
        self.low_dm_bins  = range( 0, 19)
        self.high_dm_bins = range(22, 46)
        self.all_bins     = self.low_dm_bins + self.high_dm_bins
        self.low_dm_bins_met_250to300  = [4, 5, 8, 9, 10, 11]
        self.low_dm_bins_met_300toINF  = [17, 18]
        self.low_dm_bins_met_250to400  = [0, 1, 2, 3, 6, 7, 12, 13, 14]
        self.low_dm_bins_met_400toINF  = [15, 16]
        self.high_dm_bins_met_250to400 = range(22, 46, 2)
        self.high_dm_bins_met_400toINF = range(23, 46, 2)

    def getValues(self):
        for era in self.eras:
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

            # print values 
            print era
            for b in self.all_bins:
                print "bin {0}: N = {1:.3f} +/- {2:.3f} S = {3:.3f} +/- {4:.3f}".format(b, 
                        self.binValues[era][b]["norm"],  self.binValues[era][b]["norm_error"], 
                        self.binValues[era][b]["shape"], self.binValues[era][b]["shape_error"])



