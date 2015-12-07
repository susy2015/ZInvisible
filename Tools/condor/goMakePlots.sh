#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

xrdcp root://cmseos.fnal.gov/$(echo $5 | sed 's|/eos/uscms||') .

./makePlots --condor -D $1 -N $3 -M $4 -L $5

rm $(echo $5 | sed 's|.*/||')

# Run 2015B
#./makePlots --condor -D $1 -N $3 -M $4 -L 40.838

# Run 2015C
#./makePlots --condor -D $1 -N $3 -M $4 -L 15.522

# Run 2015D
#./makePlots --condor -D $1 -N $3 -M $4 -L 209.797

# Run 2015CD
#./makePlots --condor -D $1 -N $3 -M $4 -L 225.319
#./makePlots --condor -D $1 -N $3 -M $4 -L 569.495

