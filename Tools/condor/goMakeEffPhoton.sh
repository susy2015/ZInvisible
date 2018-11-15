#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#get the release setup and in place
CMSSW_TARBALL="$2.tag.gz"
if [ -f $CMSSW_TARBALL ]; then
    tar -xzf $CMSSW_TARBALL
else
    echo "ERROR in goMakeEffPhoton.sh: The tarball $CMSSW_TARBALL does not exist. Exiting now."
    #exit 1
fi

cd $2/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
FILES_TARBALL="${_CONDOR_SCRATCH_DIR}/gmep.tar.gz"
if [ -f $FILES_TARBALL ]; then
    tar -xzf $FILES_TARBALL
else
    echo "ERROR in goMakeEffPhoton.sh: The tarball $FILES_TARBALL does not exist. Exiting now."
    #exit 1
fi

cd WORLDSWORSESOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls -lhrt

./calcEffPhoton --condor -D $1 -N $3 -M $4 

ls -lhrt

# declare array of patterns
declare -a patterns=("effhists_")

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

rm $(echo $6 | sed 's|.*/||')
rm -r ${_CONDOR_SCRATCH_DIR}/$2

