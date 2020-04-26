# ZInvisible

## Overview

These instructions will walk through downloading and setting up the following repositories that are required to run the ZInvisible framework.
- [ZInvisible](https://github.com/susy2015/ZInvisible)
- [SusyAnaTools](https://github.com/susy2015/SusyAnaTools)
- [TopTagger](https://github.com/susy2015/TopTagger)
- [TopTaggerTools](https://github.com/susy2015/TopTaggerTools)

The ZInvisible framework also uses configuration files from the StopCfg and TopTaggerCfg repositories. Here are links to the available releases.
- [StopCfg](https://github.com/susy2015/StopCfg/releases)
- [TopTaggerCfg](https://github.com/susy2015/TopTaggerCfg/releases)

## Setup TopTagger and SusyAnaTools 

First follow the instructions [here](https://github.com/susy2015/SusyAnaTools#instructions). Once you are done, you should have a CMSSW area setup that contains the TopTagger, TopTaggerTools, and SusyAnaTools repositories.

<details> 

<summary> TopTagger for Standalone (edm free) </summary>

If you want to install the TopTagger for Standalone (edm free), follow the instructions within CMSSW [here](https://github.com/susy2015/TopTagger/tree/master/TopTagger#standalone-edm-free-install-instructions-within-cmssw), but exclude the commands you have already done (don't repeat CMSSW setup and cloning TopTagger repository).

</details>

## Setup ZInvisible

Go your CMSSW area which you should have already setup (see instructions [here](https://github.com/susy2015/SusyAnaTools#instructions). We recommend using CMSSW_10_2_9, which has support for uproot. The command cmsenv will need to be run during every new terminal session.
```
cd CMSSW_10_2_9
cmsenv
```

**WARNING:** It is unwise to rename the path of your CMSSW area (e.g. any directory in your path before the CMSSW_10_2_9 directory) after checking out a CMSSW release. It will cause things to break because CMSSW_BASE will still be set to the old directory name.


Checkout the ZInvisible repository.
```
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/ZInvisible.git
cd ZInvisible/Tools
```

Checkout the NanoSUSY-tools repository and make softlinks of the trigger efficiency root files.
```
cd $CMSSW_BASE/src
git clone -b postpro_v3.0 git@github.com:susy2015/NanoSUSY-tools.git PhysicsTools/NanoSUSYTools
cd $CMSSW_BASE/src/ZInvisible/Tools
ln -s $CMSSW_BASE/src/PhysicsTools/NanoSUSYTools/data/trigger_eff/*.root .
```

## Get Configuration Files


Go to the `ZInvisible/Tools` directory and checkout the config files using getTaggerCfg.sh and getStopCfg.sh.

The top tagger implementation has changed between v5 and v6 ntuples. Use the config files described bellow.

Configuration files for using v5 ntuples:

```
cd $CMSSW_BASE/src/ZInvisible/Tools
mkdir ../../myTopTaggerCfgs
mkdir ../../myStopCfgs
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2016_V1.1.0 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2016.cfg -o 
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2017_V1.1.0 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2017.cfg -o 
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2018_V1.1.0 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2018.cfg -o 
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2016_V2.0.0 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2016.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2017_V2.0.0 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2017.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2018_V2.0.0 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2018.cfg -o
getStopCfg.sh -d ../../myStopCfgs -t PostProcessed_StopNtuple_v5.0.1 -o
```

Configuration files for using v6 ntuples:

```
cd $CMSSW_BASE/src/ZInvisible/Tools
mkdir ../../myTopTaggerCfgs
mkdir ../../myStopCfgs
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2016_V1.1.2 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2016.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2017_V1.1.2 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2017.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2018_V1.1.2 -d ../../myTopTaggerCfgs/ -f TopTagger_Tensorflow_2018.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2016_V2.0.1 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2016.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2017_V2.0.1 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2017.cfg -o
getTaggerCfg.sh -t DeepCombined_DeepCSV_GR_fromStop0lPostProc_2018_V2.0.1 -d ../../myTopTaggerCfgs/ -f TopTagger_DiscriminatorFilter_2018.cfg -o
getStopCfg.sh -d ../../myStopCfgs -t PostProcessed_StopNtuple_v6.1.1 -o
```




There are more detailed instructions that you can reference [here](https://github.com/susy2015/SusyAnaTools#get-configuration-files).

Make sure that you checkout the configuration files in the `$CMSSW_BASE/src/ZInvisible/Tools` directory (with softlinks if you use getTaggerCfg.sh and getStopCfg.sh). You may specify a different directory for the area where the release is downloaded, as the softlinks will point to that location.

## Get More Dependencies

Run the following commands to get additional dependencies.

Get the top tagger releases used in post-processing and create soft links to top tagger efficiency root files.
```
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6
ln -s TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6/tTagEff_2016.root .
ln -s TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6/tTagEff_2017.root .
ln -s TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6/tTagEff_2018.root .
```

Get json library to use json files in C++.
```
cd $CMSSW_BASE/src
git clone git@github.com:nlohmann/json.git
cd json
git checkout v3.7.3
```

Get EstToolsSUSY repository which has the unit bin json file dc_BkgPred_BinMaps_master.json.
```
cd $CMSSW_BASE/src
git clone git@github.com:mkilpatr/EstToolsSUSY.git
cd EstToolsSUSY
git checkout SBv4
```

Make a softlink to dc_BkgPred_BinMaps_master.json.
```
cd $CMSSW_BASE/src/ZInvisible/Tools
ln -s $CMSSW_BASE/src/EstToolsSUSY/SUSYNano19/dc_BkgPred_BinMaps_master.json .
```

Copy these root files.
```
cd $CMSSW_BASE/src/ZInvisible/Tools
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/2016_trigger_eff.root .
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/2017_trigger_eff.root .
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/2018_trigger_eff.root .
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/shapes_njets.root .
```

## Setup

Go to your working area and then setup the TopTagger environment. This will need to be run after running `cmsenv` in every new terminal session. Use the command for your shell (bash or tcsh). Type "echo $SHELL" to check your shell if you don't know it.

Go to working area:
```
cd $CMSSW_BASE/src/ZInvisible/Tools
```

For bash users:
```
source $CMSSW_BASE/src/TopTagger/TopTagger/test/taggerSetup.sh
```

For tcsh users:
```
source $CMSSW_BASE/src/TopTagger/TopTagger/test/taggerSetup.csh
```

Now compile. 
```
cd $CMSSW_BASE/src/ZInvisible/Tools
mkdir obj
mkdir plots
make -j8
```

You will need to compile after making changes to the source code (.cc, .h, Makefile, etc). If you change the Makefile, you will need to run `make clean` and then `make`.

Copy the following root files from Caleb's area. These files will need to be updated for our full Run 2 analysis (if we still use them).
```
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/effhists_GJets.root .
cp /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/syst_all.root .
```

## Running Scripts

Now try running makePlots.
```
time ./makePlots -D DYJetsToLL_HT_400to600_2016 -Y 2016 -E 1000 -R Data_MET_2016 | grep -v LHAPDF
```
The `-D` option is for the dataset. The `-E` option is for number of events to process.

This script should output some pdf/png plots. The `-st` options can be used for only saving the root file and not making pdf/png files.

```
time ./makePlots -st -D DYJetsToLL_HT_400to600_2016 -Y 2016 -E 1000 -R Data_MET_2016 | grep -v LHAPDF
```

If `makePlots` succeeds it will create a file named `histoutput.root`. You can open this file with a TBrowser either on cmslpc or by copying it to your machine.

Here is an example of using rsync to copy a root file to your computer. Replace USERNAME with your cmslpc username. Replace WORKINGAREA the cmslpc path to the root file which you can get with `pwd`.
```
rsync -avz USERNAME@cmslpc-sl6.fnal.gov:WORKINGAREA/histoutput.root .
```
Here is an example of this command.
```
rsync -avz caleb@cmslpc-sl6.fnal.gov:/uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/histoutput.root .
```

Open the root file in TBrowser and click on the directories to view histograms.
```
root histoutput.root
TBrowser b
```

You can also replace ZJetsToNuNu with DYJetsToLL, TTbarDiLep, TEST, and other datasets.
The "-l" option will make lepton plots, and the "-g" option will make photon plots.
```
./makePlots -D DYJetsToLL -E 1000 -l
./makePlots -D GJets -E 1000 -g
```

There is also a script for running over multiple samples and a small number of events without using condor. The year is the only argument.
```
./quickPlot.sh 2016
./quickPlot.sh 2017
./quickPlot.sh 2018
```

Different samples and number of events can be specified by editing quickPlot.sh.

You can also try running the moneyplot script.
```
./moneyplot
```
The moneyplot script should create some text, pdf, and png files that contain moneyplot in the name. The plots show the the central value of the ZInsibile background estimation for each of the 84 search bins and the accosiated uncertainties. The text files show the contents of each bin and the associated unertainties.

## Condor

If you have reached this point, you probably want to try running on condor. You can use mutliple MC and Data sets on condor, and you can run many jobs to run over a lot of data. Condor will return multiple root files. These can been added together and then plotted.

First you need to setup your CMS voms proxy. Here is the command for setting up a one week long proxy (168 hours).
```
voms-proxy-init --valid 168:00 -voms cms
```
Here is the command for checking your current proxy and remaining time.
```
voms-proxy-info
```

Here is an example of submitting a condor job.
```
cd $CMSSW_BASE/src/ZInvisible/Tools/condor

# list available sample sets
python condorSubmit.py -l

# submit specific sample sets
python condorSubmit.py -d GJets,ZJetsToNuNu
```

Here are some useful condor commands.
```
condor_q                    check condor jobs
condor_rm USERNAME          remove all condor jobs (replace USERNAME with your username)
condor_userprio             check condor priority (lower values are better)
watch "condor_q | tail"     watch job status; use CRTL-C to stop
```

Your condor jobs should produce log, stdout, and stderr files for each job in the logs directory. You can check these for errors.

### Submitting Data and MC Samples

- All Data and MC for Run 2: Electron, Muon, and Photon Control Regions
- Use TTbarNoHad (leptonic ttabr) for the lepton control region and TTbar (inclusive ttbar) for the photon control region. 
- 2017 is split into Periods B-E and F
- 2018 is split into PreHem and PostHEM

```
year="2016"; period="";         python condorSubmit.py -d Data_SingleElectron_${year}${period},Data_SingleMuon_${year}${period},Data_SinglePhoton_${year}${period},DYJetsToLL_${year},TTbarNoHad_${year},SingleTopZinv_${year},Rare_${year},TTZ_${year},Diboson_${year},GJets_${year},QCD_Photon_${year},WJetsToLNu_${year},TTbar_${year},tW_${year},ZJetsToNuNu_${year} -n 2 -y ${year}${period} -r Data_MET_${year}${period}
year="2017"; period="_BE";      python condorSubmit.py -d Data_SingleElectron_${year}${period},Data_SingleMuon_${year}${period},Data_SinglePhoton_${year}${period},DYJetsToLL_${year},TTbarNoHad_${year},SingleTopZinv_${year},Rare_${year},TTZ_${year},Diboson_${year},GJets_${year},QCD_Photon_${year},WJetsToLNu_${year},TTbar_${year},tW_${year},ZJetsToNuNu_${year} -n 2 -y ${year}${period} -r Data_MET_${year}${period}
year="2017"; period="_F";       python condorSubmit.py -d Data_SingleElectron_${year}${period},Data_SingleMuon_${year}${period},Data_SinglePhoton_${year}${period},DYJetsToLL_${year},TTbarNoHad_${year},SingleTopZinv_${year},Rare_${year},TTZ_${year},Diboson_${year},GJets_${year},QCD_Photon_${year},WJetsToLNu_${year},TTbar_${year},tW_${year},ZJetsToNuNu_${year} -n 2 -y ${year}${period} -r Data_MET_${year}${period}
year="2018"; period="_PreHEM";  python condorSubmit.py -d Data_EGamma_${year}${period},Data_SingleMuon_${year}${period},DYJetsToLL_${year},TTbarNoHad_${year},SingleTopZinv_${year},Rare_${year},TTZ_${year},Diboson_${year},GJets_${year},QCD_Photon_${year},WJetsToLNu_${year},TTbar_${year},tW_${year},ZJetsToNuNu_${year} -n 2 -y ${year}${period} -r Data_MET_${year}${period}
year="2018"; period="_PostHEM"; python condorSubmit.py -d Data_EGamma_${year}${period},Data_SingleMuon_${year}${period},DYJetsToLL_${year},TTbarNoHad_${year},SingleTopZinv_${year},Rare_${year},TTZ_${year},Diboson_${year},GJets_${year},QCD_Photon_${year},WJetsToLNu_${year},TTbar_${year},tW_${year},ZJetsToNuNu_${year} -n 2 -y ${year}${period} -r Data_MET_${year}${period}
```

There is now a script to submit all Run 2 Data and MC with one command.
```
python multiSubmit.py
```

The script multiSubmit.py will create a json file in the runs directory that lists the condor submission directories. This json file will be used again to process the output of the condor jobs.

### Process Output from Condor

First get ahadd.py. This is a script that uses multithreading to add histograms together by applying hadd in parallel.
```
mkdir ~/bin
cd ~/bin
git@github.com:pastika/Alex-Hadd.git
ln -s Alex-Hadd/ahadd.py .
```

Add `~/bin` to your path. I also have `SusyAnaTools/Tools/scripts/` in my path. You can add this line to your ~/.bash_profile using your USERNAME and PATH_TO_WORKING_AREA.
```
export PATH="$PATH:/uscms/home/USERNAME/bin/:/uscms/home/USERNAME/PATH_TO_WORKING_AREA/SusyAnaTools/Tools/scripts/"
```

Then source your ~/.bash_profile to apply this change to your PATH.
```
source ~/.bash_profile
```

If your condor jobs complete successfully, the jobs will output root files to the condor directory. We have a script called processResults.sh to hadd the output files, which adds the root histograms together and produces one file.

**WARNING:** The processResults.sh script moves all root files from your current directory to a new directory and adds them together. Make sure your condor directory only contains root files that were output from condor and that you want to add together.

The script requires two arguments in this order: name, year . The name will be used to create the directory that stores all the root files. The year should be 2016, 2017, or 2018 corresponding to the Data/MC year.
```
cd $CMSSW_BASE/src/ZInvisible/Tools/condor
./processResults.sh NAME_FOR_DIRECTORY YEAR
```

Here is an example.
```
cd $CMSSW_BASE/src/ZInvisible/Tools/condor
./processResults.sh ZJetsToNuNu 2016
```

The script processResults.sh will 
- create a new directory with the name you provide, the date, and a version number for that date (1,2,3...).
- move all root files from current directory to the new directory and hadd them (checking for broken files)
- copy the resulting root file to ZInvisible/Tools
- run MakePlots using this root file

This should generate some pdf and png files in the ZInvisible/Tools/plots directory, which you may rsync and view as desired.

Process 2016 results:
```
./processResults.sh PhotonAndMuonControlRegionSelection_2016 2016
```

Process 2017 results:
```
./processResults.sh PhotonAndMuonControlRegionSelection_2017 2017
```

There is now a script to process all Run 2 output from condor. It accepts a json file containing the condor submission directories (and created by multiSubmit.py) as an input.

```
python python/process.py -j runs/submission_2019-08-28_14-24-48.json
```

## Running Prediction Modules
Use the run_modules.py script to run modules to calculate the normalization and shape factors, as well as the Z to invisible prediction in the validation and search bins.
```
python python/run_modules.py -j runs/submission_2019-08-28_14-24-48.json
```

## La Fin

If everything has worked up to this point, you have arrived at The End... for now. You will go far, my friend. Otherwise, don't worry. Keep trying and contact an expert if you cannot resolve the issues. Feel free to post issues on the issues page of this repository. Also, you are welcome to help solve the current issues if you have time.



