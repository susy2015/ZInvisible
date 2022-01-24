#!/opt/homebrew/bin/bash

# getMoreResults.sh

# rsync plots and latex files
cmssw="CMSSW_10_2_9"
user=caleb
host=cmslpc-sl7.fnal.gov
# nobackup path
path=/uscms/home/caleb/nobackup/SusyAnalysis/$cmssw/src/ZInvisible/Tools
# scratch path
#path=/uscmst1b_scratch/lpc1/3DayLifetime/caleb/$cmssw/src/ZInvisible/Tools
rsync -az $user@$host:$path/more_plots .
rsync -az $user@$host:$path/syst_plots .
rsync -az $user@$host:$path/latex_files .
rsync -az $user@$host:$path/comparison_plots .

