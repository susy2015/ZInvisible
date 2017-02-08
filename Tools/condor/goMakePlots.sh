#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#cd $2/src
#eval `scramv1 runtime -sh`

#cd ${_CONDOR_SCRATCH_DIR}
#cd /uscms_data/d3/snorberg/CMSSW_8_0_23_patch1/src/ZInvisible/Tools/

tar -xzf CMSSW_8_0_23_patch1.tar.gz
cd CMSSW_8_0_23_patch1/src/ZInvisible/Tools/
scram b ProjectRename
eval `scramv1 runtime -sh`
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

source /cvmfs/cms.cern.ch/cmsset_default.sh
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_23/src/SusyAnaTools/Tools/obj:${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_23/src/TopTagger/TopTagger/test:${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_23/src/opencv/lib:${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_23/lib/slc6_amd64_gcc530
# cmsenv
#eval `scramv1 runtime -sh`

#cd ${_CONDOR_SCRATCH_DIR}
#cd CMSSW_8_0_23_patch1/src/ZInvisible/Tools 
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls

./makePlots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_v1_2017

ls
mv histoutput_* ${_CONDOR_SCRATCH_DIR}
mv minituple_histoutput_* ${_CONDOR_SCRATCH_DIR} 
rm $(echo $6 | sed 's|.*/||')

