# ZInvisible

## Overview

These instructions will walk through downloading and setting up the following repositories that are required to run the ZInvisible framework.
- [ZInvisible](https://github.com/susy2015/ZInvisible/tree/NanoAOD)
- [SusyAnaTools](https://github.com/susy2015/SusyAnaTools/tree/NanoAOD)
- [TopTagger](https://github.com/susy2015/TopTagger)
- [TopTaggerTools](https://github.com/susy2015/TopTaggerTools)

The ZInvisible framework also uses configuration files from the StopCfg and TopTaggerCfg repositories. Here are links to the available releases.
- [StopCfg](https://github.com/susy2015/StopCfg/releases)
- [TopTaggerCfg](https://github.com/susy2015/TopTaggerCfg/releases)

## Setup TopTagger and SusyAnaTools 

First follow the instructions [here](https://github.com/susy2015/SusyAnaTools/tree/NanoAOD#instructions). Once you are done, you should have a CMSSW area setup that contains the TopTagger, TopTaggerTools, and SusyAnaTools repositories.

<details> 

<summary> TopTagger for Standalone (edm free) </summary>

If you want to install the TopTagger for Standalone (edm free), follow the instructions within CMSSW [here](https://github.com/susy2015/TopTagger/tree/master/TopTagger#standalone-edm-free-install-instructions-within-cmssw), but exclude the commands you have already done (don't repeat CMSSW setup and cloning TopTagger repository).

</details>

## Setup ZInvisible

Go your CMSSW area which you should have already setup (see instructions [here](https://github.com/susy2015/SusyAnaTools/tree/NanoAOD#instructions). We recommend using CMSSW_10_2_9, which has support for uproot. The command cmsenv will need to be run during every new terminal session.
```
cd CMSSW_10_2_9
cmsenv
```

WARNING: It is unwise to rename the path of your CMSSW area (e.g. any directory in your path before the CMSSW_10_2_9 directory) after checking out a CMSSW release. It will cause things to break because CMSSW_BASE will still be set to the old directory name.

Checkout the ZInvisible repository.
```
cd $CMSSW_BASE/src
git clone git@github.com:susy2015/ZInvisible.git
cd ZInvisible/Tools
git checkout NanoAOD
```

## Get Configuration Files

Go to the `ZInvisible/Tools` directory.
```
cd $CMSSW_BASE/src/ZInvisible/Tools
mkdir ../../myTopTaggerCfgs
mkdir ../../myStopCfgs
getTaggerCfg.sh -t DeepCombined_fromNanoAOD_RES_T_DeepAK8_T_v1.0.1 -d ../../myTopTaggerCfgs -o
getStopCfg.sh -t PostProcessed_StopTuple_V2.9.0 -d ../../myStopCfgs -o
```

There are more detailed instructions that you can reference [here](https://github.com/susy2015/SusyAnaTools/tree/NanoAOD#get-configuration-files).

Make sure that you checkout the configuration files in the `$CMSSW_BASE/src/ZInvisible/Tools` directory (with softlinks if you use getTaggerCfg.sh and getStopCfg.sh). You may specify a different directory for the area where the release is downloaded, as the softlinks will point to that location.

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
./makePlots -D ZJetsToNuNu_2016 -Y 2016 -E 1000 | grep -v LHAPDF
```

The `-D` option is for the dataset (ZJetsToNuNu). The `-E` option is for number of events to process (1000).

This script should output some pdf/png plots. The `-s` option can be used for only saving the root file and not making pdf/png files.

You can also run over a specific HT sample range.
```
./makePlots -D ZJetsToNuNu_HT_100to200_2016 -E 1000 -Y 2016
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

There is also a script for running over multiple samples and a small number of events without using condor. The year is the only argument. Here are the commands for 2016 and 2017.
```
./quickPlot.sh 2016
./quickPlot.sh 2017
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

### Data/MC Samples

Submit 2016 MC only:
```
python condorSubmit.py -d DYJetsToLL_2016,TTbarNoHad_2016,SingleTopZinv_2016,Rare_2016,TTZ_2016,Diboson_2016,GJets_2016,QCD_2016,WJetsToLNu_2016,TTbarAll_2016,tW_2016,ZJetsToNuNu_2016 -n 1 -y 2016
```

Submit 2016 Data only:
```
python condorSubmit.py -d Data_SingleMuon_2016,Data_SinglePhoton_2016 -n 1 -y 2016
```

Submit 2016 Data and MC:
```
python condorSubmit.py -d Data_SingleMuon_2016,DYJetsToLL_2016,TTbarNoHad_2016,SingleTopZinv_2016,Rare_2016,TTZ_2016,Diboson_2016,Data_SinglePhoton_2016,GJets_2016,QCD_2016,WJetsToLNu_2016,TTbarAll_2016,tW_2016,ZJetsToNuNu_2016 -n 1 -y 2016
```

Submit 2017 MC only:
```
python condorSubmit.py -d DYJetsToLL_2017,TTbarNoHad_2017,SingleTopZinv_2017,Rare_2017,TTZ_2017,Diboson_2017,GJets_2017,QCD_2017,WJetsToLNu_2017,TTbarAll_2017,tW_2017,ZJetsToNuNu_2017 -n 1 -y 2017
```

Submit 2017 Data only:
```
python condorSubmit.py -d Data_SingleMuon_2017,Data_SinglePhoton_2017 -n 1 -y 2017
```

Submit 2017 Data and MC:
```
python condorSubmit.py -d Data_SingleMuon_2017,DYJetsToLL_2017,TTbarNoHad_2017,SingleTopZinv_2017,Rare_2017,TTZ_2017,Diboson_2017,Data_SinglePhoton_2017,GJets_2017,QCD_2017,WJetsToLNu_2017,TTbarAll_2017,tW_2017,ZJetsToNuNu_2017 -n 1 -y 2017
```

Submit 2017 Muon Data and MC:
```
python condorSubmit.py -d Data_SingleMuon_2017,DYJetsToLL_2017,TTbarNoHad_2017,SingleTopZinv_2017,Rare_2017,TTZ_2017,Diboson_2017,ZJetsToNuNu_2017 -n 1 -y 2017
```

Submit 2017 Photon Data and MC:
```
python condorSubmit.py -d Rare_2017,TTZ_2017,Diboson_2017,Data_SinglePhoton_2017,GJets_2017,QCD_2017,WJetsToLNu_2017,TTbarAll_2017,tW_2017,ZJetsToNuNu_2017 -n 1 -y 2017
```


### Process Output from Condor

If your jobs complete successfully, the jobs will output root files to the condor directory. We have a script called processResults.sh to hadd the output files, which adds the root histograms together and produces one file.

WARNING: The processResults.sh script moves all root files from your current directory to a new directory and adds them together. Make sure your condor directory only contains root files that were output from condor and that you want to add together.

The script requires two arguments in this order: name, year . The name will be used to create the directory that stores all the root files. The year should be 2016, 2017, or 2018 corresponding to the Data/MC year.
```
./processResults.sh NAME_FOR_DIRECTORY YEAR
```

Here is an example.
```
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

## La Fin

If everything has worked up to this point, you have arrived at The End... for now. You will go far, my friend. Otherwise, don't worry. Keep trying and contact an expert if you cannot resolve the issues. Feel free to post issues on the issues page of this repository. Also, you are welcome to help solve the current issues if you have time.



