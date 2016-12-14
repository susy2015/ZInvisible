#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls

./makePlots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_59_2016

ls

rm $(echo $6 | sed 's|.*/||')

