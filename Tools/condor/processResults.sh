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
dirName=$1
year=$2
executableOption=$3

executable=
resultFile=
dirPrefix=

dataDir="output"
brokenDir="broken_files" # directory for broken root files
zinvDir=$CMSSW_BASE/src/ZInvisible/Tools
condorDir=$zinvDir/condor

if [ -z "$dirName" ]; then
    echo "Please provide a directory name as the first argument."
    echo "  Example: DataMC_2016_submission_2019-05-05_21-57-41"
    exit 1
fi

if [[ "$year" != "2016" && "$year" != "2017" && "$year" != "2018" && "$year" != "2018_AB" && "$year" != "2018_CD" ]]
then
    echo "Please enter 2016, 2017, 2018, 2018_AB or 2018_CD for the year as the second argument."
    exit 1
fi

if [ "$executableOption" = "-c" ]; then
    resultFile="effhists_"$dirName".root"
    executable="echo nothing to do"
    dirPrefix="effhists"
else
    resultFile="result.root"
    executable="./makePlots -f -I $resultFile -Y $year | grep -v LHAPDF"
    dirPrefix="histos"
fi

echo "- Running processResults.sh for the data set $dirName"

# old version: make new data directory
#today=$(date '+%d_%b_%Y')
#i=1
#dataDir=""$dirPrefix"_"$dirName"_"$today"_"$i""

# go to directory to process
cd $condorDir/$dirName

# check if there are any root files here
if [[ ! $(ls *.root) ]]; then
    echo "- ERROR: There are no root files in the directory $condorDir/$dirName to process."
    exit 1
fi

# check if data directory exists
#while [[ -d $dataDir ]]
#do
#    echo "- Found directory $dataDir"
#    i=$[$i+1]
#    dataDir=""$dirPrefix"_"$dirName"_"$today"_"$i""
#done

# create unique data directory
echo "- Creating directory $dataDir"
mkdir $dataDir

# move root files to data directory
echo "- Moving root files to $dataDir"
mv *.root $dataDir

# number of files
numFiles=$(ls $dataDir/*.root | wc -l | cut -f1 -d " ")
echo "Number of files returned from condor: $numFiles"

# add histograms
cd $dataDir
error=1
# while there are errors, attempt hadd
while [[ $error != 0 ]]
do
    echo "- Executing hadd to add histograms to create $resultFile; see hadd.log for stdout and stderr"
    #time hadd $resultFile *.root &> hadd.log
    time ahadd.py $resultFile *.root &> hadd.log
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
        echo "- Moving broken file $file and hadd.log to $brokenDir"
        mv $file $brokenDir
        mv hadd.log $brokenDir
        echo "- Removing $resultFile"
        rm $resultFile
    else
        echo "  hadd command succeeded; $resultFile should be ready to use"
    fi
done

echo "- Copying $resultFile to $zinvDir and moving $resultFile to $condorDir/$dirName"
if [[ -f $resultFile ]]; then
    cp $resultFile $zinvDir
    mv $resultFile $condorDir/$dirName
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

echo "- Running Exacutable"
# use eval for running command that pipes to another command
echo $executable
eval $executable

if [[ $? == 0 ]]; then
    echo "  Exacutable was successful"
else
    echo "  ERROR: Exacutable failed"
fi

