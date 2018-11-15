#!/bin/bash

# processResults.sh
# Author: Caleb Smith
# Date: 15-Aug-2018

# process root files returned by condor job
# - place files into directory named with data set and date
# - add histograms using hadd to produce summed result
# - if there are errors adding histograms due to a specific file, 
#   we move the bad file to another directory and attempt to add histograms again 
# - repeat this process until hadd succeeds
# - make plots using summed result


# options
dataSet=$1
executableOption=$2

exacutable=
resultFile=
dirName=

zinvDir=$CMSSW_BASE/src/ZInvisible/Tools
condorDir=$zinvDir/condor

if [ -z "$dataSet" ]; then
    echo "- ERROR: Please provide a data set for the directory name."
    echo "  Examples: DYJetsToLL, GJets, ZJetsToNuNu, etc."
    exit 1
fi

if [ "$executableOption" = "-c" ]; then
    resultFile="effhists_"$dataSet".root"
    exacutable="echo nothing to do"
    dirName="effhists"
else
    resultFile="result.root"
    exacutable="./makePlots -f -I $resultFile"
    dirName="histos"
fi

echo "- Running processResults.sh for the data set $dataSet"

# data directory
today=$(date '+%d_%b_%Y')
i=1
dataDir=""$dirName"_"$dataSet"_"$today"_"$i""

# directory for broken files
brokenDir="broken_files"

# go to condor directory
cd $condorDir

# check if there are any root files here
if [[ ! $(ls *.root) ]]; then
    echo "- ERROR: There are no root files in the directory $condorDir to process."
    exit 1
fi

# check if data directory exists
while [[ -d $dataDir ]]
do
    echo "- Found directory $dataDir"
    i=$[$i+1]
    dataDir=""$dirName"_"$1"_"$today"_"$i""
done

# create unique data directory
echo "- Creating directory $dataDir"
mkdir $dataDir

# move root files to data directory
echo "- Moving root files to $dataDir"
mv *.root $dataDir

# add histograms
cd $dataDir
error=1
# while there are errors, attempt hadd
while [[ $error != 0 ]]
do
    echo "- Executing hadd to add histograms to create $resultFile; see hadd.log for stdout and stderr"
    hadd $resultFile *.root &> hadd.log
    error=$?
    if [[ $error != 0 ]]; then
        echo "  ERROR: hadd command failed"
        if [[ ! -d $brokenDir ]]; then
            mkdir $brokenDir
        fi
        substring="hadd exiting due to error in"
        message=$(grep "$substring" hadd.log)
        file=$(expr match "$message" "$substring \(.*\)")
        echo "    message: $message"
        echo "    file: $file"
        echo "- Moving broken file $file to $brokenDir"
        mv $file $brokenDir
        echo "- Removing $resultFile"
        rm $resultFile
    else
        echo "  hadd command succeeded; $resultFile should be ready to use"
    fi
done

echo "- Copying $resultFile to $zinvDir"
if [[ -f $resultFile ]]; then
    cp $resultFile $zinvDir
else
    echo "  ERROR: the file $resultFile does not exit"
    exit 1
fi

# go to zinv area
cd $zinvDir

# remove old plots if they exist 
if [ -z "$(ls -A plots)" ]; then
    echo "- There are no old plots to remove"
else
    echo "- Removing old plots"
    rm plots/*
fi

echo "- Compiling Exacutable"
make -j8

#echo "- Running Exacutable with -g for photons"
#./makePlots -f -g -I $resultFile

echo "- Running Exacutable"
#./makePlots -f -I $resultFile
$exacutable

if [[ $? == 0 ]]; then
    echo "  Exacutable was successful"
else
    echo "  ERROR: Exacutable failed"
fi

