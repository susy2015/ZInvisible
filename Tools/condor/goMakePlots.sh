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
tar -xzf ${_CONDOR_SCRATCH_DIR}/gmp.tar.gz
cd WORLDSWORSESOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls -lhrt

#./makeDataMCplots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_v1_2017
./makePlots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_v1_2017 

ls -lhrt

# declare array of patterns
declare -a patterns=("histoutput_" "minituple_histoutput_")

# for each pattern, check that files beginning with pattern exist and move them if they do
for pattern in "${patterns[@]}"
do
    if ls $pattern* &> /dev/null
    then
        mv $pattern* ${_CONDOR_SCRATCH_DIR}
    else
        echo "There were no files found beginning with $pattern"
    fi
done

#mv histoutput_* ${_CONDOR_SCRATCH_DIR}
#mv minituple_histoutput_* ${_CONDOR_SCRATCH_DIR}

rm $(echo $6 | sed 's|.*/||')
rm -r ${_CONDOR_SCRATCH_DIR}/$2

