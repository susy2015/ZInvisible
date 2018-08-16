#!/bin/bash

# processResults.sh
# Author: Caleb Smith
# Date: 15-Aug-2018

# process root files returned by condor job
# - place files into directory named with data set and date
# - add histograms to produce summed result
# - if there are errors adding histograms, move bad file to another directory 
#   and attampt to add histograms again 
# - make plots using summed result


dataSet=$1

if [ -z "$dataSet" ]; then
    echo "ERROR: Please provide a data set for the directory name."
    exit 1
fi

# ZInvisible directory
zinvDir=$CMSSW_BASE/src/ZInvisible/Tools

# condor directory
condorDir=$zinvDir/condor

# data directory
today=$(date '+%d_%b_%Y')
i=1
dataDir="histos_"$dataSet"_"$today"_"$i""

# directory for broken files
brokenDir="broken_files"

# result file
resultFile="result.root"

# go to condor directory
cd $condorDir

# check if data directory exists
while [[ -d $dataDir ]]
do
    echo "Found directory $dataDir"
    i=$[$i+1]
    dataDir="histos_"$1"_"$today"_"$i""
done

# create unique data directory
echo "Creating directory $dataDir"
mkdir $dataDir

# move root files to data directory
echo "Moving root files to $dataDir"
mv *.root $dataDir

# add histograms
cd $dataDir
error=1
# while there are errors, attempt hadd
while [[ $error != 0 ]]
do
    hadd $resultFile *.root &> hadd.log
    error=$?
    if [[ $error != 0 ]]; then
        echo "ERROR: hadd command failed"
        if [[ ! -d $brokenDir ]]; then
            mkdir $brokenDir
        fi
        substring="hadd exiting due to error in"
        message=$(grep "$substring" hadd.log)
        file=$(expr match "$message" "$substring \(.*\)")
        echo "message: $message"
        echo "file: $file"
        mv $file $brokenDir
        rm $resultFile
    else
        echo "hadd command succeeded"
    fi
done

echo "Copying $resultFile to $zinvDir"
if [[ -f $resultFile ]]; then
    cp $resultFile $zinvDir
else
    echo "ERROR: the file $resultFile does not exit"
    exit 1
fi

echo "Removing old plots"
cd $zinvDir
rm plots/*

echo "Compiling MakePlots"
make -j8

echo "Running MakePlots"
./makePlots -f -I $resultFile

if [[ $? == 0 ]]; then
    echo "MakePlots was successful"
else
    echo "ERROR: MakePlots failed"
fi

