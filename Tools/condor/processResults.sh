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

dataDir="output"
brokenDir="broken_files" # directory for broken root files
zinvDir=$CMSSW_BASE/src/ZInvisible/Tools
condorDir=$zinvDir/condor

if [ -z "$dirName" ]; then
    echo "Please provide a directory name as the first argument."
    echo "  Example: DataMC_2016_submission_2019-05-05_21-57-41"
    exit 1
fi

if [[ "$year" != "2016" && "$year" != "2017" && "$year" != "2017_BE" && "$year" != "2017_F" && "$year" != "2018" && "$year" != "2018_PreHEM" && "$year" != "2018_PostHEM" ]]
then
    echo "Please enter 2016, 2017, 2017_BE, 2017_F, 2018, 2018_PreHEM or 2018_PostHEM for the year as the second argument."
    exit 1
fi

if [ "$executableOption" = "-c" ]; then
    resultFile="effhists_${dirName}.root"
    executable="echo nothing to do"
elif [ "$executableOption" = "-n" ]; then
    resultFile="result.root"
    executable="echo nothing to do"
else
    resultFile="result.root"
    executable="./makePlots -q -f -I ${resultFile} -Y ${year} -R Data_MET_${year} | grep -v LHAPDF"
fi

echo "- Running processResults.sh for the data set ${dirName}"

# go to directory to process
cd $condorDir/$dirName

# check if there are any root files here
if [[ ! $(ls *.root) ]]; then
    echo "- ERROR: There are no root files in the directory $condorDir/$dirName to process."
    exit 1
fi

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

