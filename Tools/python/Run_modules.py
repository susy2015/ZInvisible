# run_modules.py

import os
from Norm_lepton_zmass import Normalization
from Shape_photon_met import Shape
from Search_bins import  ValidationBins
    
era     = "2018_PreHEM"
verbose  = False 
doUnits  = False
draw     = False

out_dir     = "validation/" 
result_file = "result.root"

#-------------------------------------------------------
# Class instanceses construction
#-------------------------------------------------------

N = Normalization(out_dir, verbose)
S = Shape(out_dir, draw, doUnits, verbose)

VB = ValidationBins(N, S, era, out_dir, verbose)

#-------------------------------------------------------
# Calculate normalization and shape factors
#-------------------------------------------------------

#syst = ""
#syst = "btag_syst_up"
syst = "btag_syst_down"
#syst = "dielec_sf_syst_down"
#syst = "dielec_trig_syst_up"
#syst = "eff_sb_syst_up"
#syst = "prefire_syst_up"
#syst = "dimu_sf_syst_down"


N.getNormAndError(result_file, syst, era)
S.getShape(result_file, syst, era)

VB.getValues(result_file, syst, era)

