#!/bin/bash

_PWD=${PWD}
#_CONDOR_SCRATCH_DIR=${PBS_O_WORKDIR}

printenv

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#move the files from /home to /var/tmp
mv ${_PWD}/$2.tar.gz ${_CONDOR_SCRATCH_DIR}
mv ${_PWD}/TT.tar.gz ${_CONDOR_SCRATCH_DIR}

cd ${_CONDOR_SCRATCH_DIR}

#get the release setup and in place
tar -xzf $2.tar.gz
cd $2/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
tar -xzf ${_CONDOR_SCRATCH_DIR}/TT.tar.gz
cd WORLDSWORSESOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}:/cvmfs/sft.cern.ch/lcg/releases/xgboost/0.60-04a4f/x86_64-slc6-gcc62-opt/lib

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

#if [ ! -f ${6##*/} ]; then
#    echo "File not found! Try xrdcp again"
#    export X509_USER_PROXY=/cms/data/$USER/.x509_user_proxy
#    xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .
#fi

pwd
ls -lhrt

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${_CONDOR_SCRATCH_DIR}

./simpleAnalyzer -c -D $1 -N $3 -M $4 -O TT  

ls -lhrt

mv TT*.root ${_CONDOR_SCRATCH_DIR}
#mv topStudyOutput_* /cms/data/pastika/condorOutput

rm $(echo $6 | sed 's|.*/||')
rm -r ${_CONDOR_SCRATCH_DIR}/$2

