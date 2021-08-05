#!/opt/homebrew/bin/bash

# getResults.sh
# Author: Caleb Smith 
# Date: 16-Aug-2018

# get results from cmslpc
# - make directory with data set and date in name
# - rsync plots

dataSet=$1
cmssw="CMSSW_10_2_9"
user=caleb
host=cmslpc-sl7.fnal.gov
# nobackup path
path=/uscms/home/caleb/nobackup/SusyAnalysis/$cmssw/src/ZInvisible/Tools/plots
# scratch path
#path=/uscmst1b_scratch/lpc1/3DayLifetime/caleb/$cmssw/src/ZInvisible/Tools/plots

if [ -z "$dataSet" ]; then
    echo "- ERROR: Please provide a data set for the directory name."
    echo "  Examples: DYJetsToLL, GJets, ZJetsToNuNu, etc."
    exit 1
fi

# data directory
today=$(date '+%d_%b_%Y')
i=1
dataDir="histos_"$dataSet"_"$today"_"$i""

# check if data directory exists
while [[ -d $dataDir ]]
do
    echo "- Found directory $dataDir"
    i=$[$i+1]
    dataDir="histos_"$1"_"$today"_"$i""
done

# create unique data directory
echo "- Creating directory $dataDir"
mkdir $dataDir

# rsync plots
echo "- Copying plots"
cd $dataDir
rsync -az $user@$host:$path/ .

