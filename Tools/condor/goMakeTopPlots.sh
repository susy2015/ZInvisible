#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#get the release setup and in place
tar -xzf $2.tar.gz
cd $2/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
tar -xzf ${_CONDOR_SCRATCH_DIR}/gtp.tar.gz
cd WORLDSWORSTSOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls -lhrt

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${_CONDOR_SCRATCH_DIR}

./makeTopPlots -s --condor -D $1 -N $3 -M $4 -L $5

ls -lhrt

mv topStudyOutput_* ${_CONDOR_SCRATCH_DIR}

rm $(echo $6 | sed 's|.*/||')
rm -r ${_CONDOR_SCRATCH_DIR}/$2

