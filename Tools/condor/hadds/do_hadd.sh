#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#get the release setup and in place
tar -xzf CMSSW_10_2_9.tar.gz
cd CMSSW_10_2_9/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
tar -xzf ${_CONDOR_SCRATCH_DIR}/gmp.tar.gz
cd WORLDSWORSTSOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}

pwd

ls -lhtr

year=$1
syst=$2

echo YEAR: ${year}
echo SYST: ${syst}

syst_files=`xrdfs root://cmseos.fnal.gov ls -u /store/user/jsw/v6p5/${year}_${syst} | grep '\.root$'`
data_files=`xrdfs root://cmseos.fnal.gov ls -u /store/user/jsw/v6p5/${year}_Data | grep '\.root$'`

echo SYST_FILES:
echo ${syst_files}

echo DATA_FILES:
echo ${data_files}

time hadd -v -f ${year}_${syst}.root ${syst_files} ${data_files} && \
time xrdcp ${year}_${syst}.root root://cmseos.fnal.gov//store/user/jsw/v6p5/ && \
rm ${year}_${syst}.root
