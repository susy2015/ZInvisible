#Shaper

import os
import numpy as np
import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Get root file explicitly

root_file = ROOT.TFile.Open("/uscms/home/arosado/nobackup/YOURWORKINGAREA/CMSSW_10_2_9/src/ZInvisible/Tools/condor/2018_AB/result.root")

in_root = 'nJets_drLeptonCleaned'

# Get histograms

year       = '2018'
particles  =  ['Electron', 'Muon']
regions    =  ['HighDM', 'LowDM']
metcuts    =  ["", 'Mid', 'Loose']


for particle in particles:
    histos = {}          # to recycle memory space?
    histos[particle] = {}
    for region in regions:
        histos[particle][region] = {}
        for metcut in metcuts: 
            prefix = 'DataMC_' + particle + '_' + region + '_' + metcut

            histos[particle][region][metcut] = {

		'DY'    :  prefix + '_nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'DYstack',
		'Sine'  :  prefix + '_nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Single tstack',
		'TTbar' :  prefix + '_nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 't#bar{t}stack',
		'Rare'  :  prefix + '_nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Rarestack',
                'Data'  :  prefix + '_nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Datadata'
		}

	    print('We are now in: ' + particle + ' ' + region + ' ' + metcut)

# Fixing bad name selection

	    if not metcut:
		histos[particle][region][metcut]['DY'] = prefix + 'nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'DYstack'
		histos[particle][region][metcut]['Sine'] = prefix + 'nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Single tstack'
		histos[particle][region][metcut]['TTbar'] = prefix + 'nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 't#bar{t}stack'
		histos[particle][region][metcut]['Rare'] = prefix + 'nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Rarestack'
		histos[particle][region][metcut]['Data'] = prefix + 'nj_' + year + 'nJets_drLeptonCleanednJets_drLeptonCleaned' + 'Datadata'

# Retriving histograms

	    h_Data   =  root_file.Get(in_root + '/' + histos[particle][region][metcut]['Data'])
	    h_DY     =  root_file.Get(in_root + '/' + histos[particle][region][metcut]['DY'])
	    h_Sine   =  root_file.Get(in_root + '/' + histos[particle][region][metcut]['Sine'])
	    h_TTbar  =  root_file.Get(in_root + '/' + histos[particle][region][metcut]['TTbar'])
	    h_Rare   =  root_file.Get(in_root + '/' + histos[particle][region][metcut]['Rare'])

# Verify if histograms exist

	    if not h_Data:
		print(in_root + '/' + histos[particle][region][metcut]['Data'] + " doesn't exist!")
	    if not h_DY:
		print(in_root + '/' + histos[particle][region][metcut]['DY'] + " doesn't exist!")
	    if not h_Sine:
		print(in_root + '/' + histos[particle][region][metcut]['Sine'] + " doesn't exist!")
	    if not h_TTbar:
		print(in_root + '/' + histos[particle][region][metcut]['TTbar'] + " doesn't exist!")
	    if not h_Rare:
		print(in_root + '/' + histos[particle][region][metcut]['Rare'] + " doesn't exist!")

# Shape factor calculation

	    h_Shape = h_Data.Clone()
	    h_Shape.GetXaxis().SetRange(2,10)
	    h_Shape.Sumw2()
	    h_Shape.Add(h_Sine, -1)
	    h_Shape.Add(h_TTbar, -1)
	    h_Shape.Add(h_Rare, -1)
	    h_Shape.Divide(h_DY)

	    bin_0 = h_Shape.GetBinContent(4)
	    print('bin 0: ' + str(bin_0))


# Create a canvas

	    canvas = ROOT.TCanvas("c", "c", 800, 800)

#Draw legend

	    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
	    legend.AddEntry(h_Shape, 'Shape', '1')
	    legend.AddEntry(h_Data, 'Data', '1')
	    legend.Draw()

#Draw histogram

	    h_Shape.Draw('error')
	#   h_Data.Draw('hist error same')

# Save new Histogram in a new root file
	    canvas.Update()
	    canvas.SaveAs(prefix + '_nj_' + year + '_AB_drLeptonCleaned_Shape.png')

