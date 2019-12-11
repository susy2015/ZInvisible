# run_modules.py

import argparse
import json
import os
from norm_lepton_zmass import Normalization  # needed
from shape_photon_met import Shape           # needed
from systematics import Systematic
from validation_bins_systematics import ValidationBins       # needed


era = "2016"
plot_dir = "more_plots"
verbose  = False
doUnits = False
draw     = False
json_file = "/uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/runs/submission_2019-12-02_15-21-40.json"


N = Normalization(plot_dir, verbose)
S = Shape(plot_dir, draw, doUnits, verbose)

with open(json_file, "r") as input_file:

    runMap = json.load(input_file)

    N.getNormAndError(result_file, era)
    S.getShape(result_file, era)
    VB = ValidationBins(N, S, eras, plot_dir, "",  verbose)

    result_file = "condor/" + runMap[era] + "/result.root"

    VB.getValues("condor/" + runMap[era] + "/result.root", era)

