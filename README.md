# ZInvisible


## Download

Checkout a CMSSW 9_X_X release of 9_3_3 or later (including 9_4_X, but not 10_X_X_).
```
cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src
cmsenv
```

The command `cmsenv` will need to be run from this area every time that you login. Make sure that the command `cmsenv` worked by checking the environment variable CMSSW_BASE. 
```
echo $CMSSW_BASE
```
The result should be the path of your CMSSW area.
```
cmslpc23.fnal.gov src $ echo $CMSSW_BASE
/uscms_data/d3/caleb/SUSY_Analysis/CMSSW_9_3_3
```
If CMSSW_BASE is empty, then `cmsenv` did not work. You can try doing `cmsenv` again from the directory CMSSW_9_3_3/src. Also, as Chris and Caleb have discovered, it is unwise to rename the directories of your working area after checking out a CMSSW release. It will cause things to break because CMSSW_BASE will still be set to the old directory name.

Clone the following repositories

```
git clone git@github.com:susy2015/ZInvisible.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:susy2015/TopTagger.git
```

For the Zinvible estimation, we are currently using the calebGJets branch of the ZInvisible repository. We are using the master branch of the SusyAnaTools repository (though the recent change of double to float is not supported yet). We are using the master branch of the TopTagger repository.




## Setup


### TopTagger

Follow the TopTagger instructions for Standalone (edm free) install instructions within CMSSW [here](https://github.com/susy2015/TopTagger/tree/master/TopTagger#standalone-edm-free-install-instructions-within-cmssw), but exclude the commands you have already done (don't repeat CMSSW setup and cloning TopTagger repository).
```
cd $CMSSW_BASE/srcTopTagger/TopTagger/test
./configure
make -j8
```

### Zinvible

Checkout the branch calebGJets.
```
cd $CMSSW_BASE/src/Zinvible/Tools
git fetch origin
git checkout calebGJets
```
Setup the TopTagger environment. This will need to be run after running `cmsenv` everytime 
```
source $CMSSW_BASE/src/TopTagger/TopTagger/test/taggerSetup.sh
```
Checkout a compatible version of the TopTagger. This will only need to be done once in this working area (unless you want to use a different TopTagger version). The commands `tcsh` and `bash` should only be used if you are using bash. This is because there is a `setup.csh` script, but there is no `setup.sh` available.

If you are using bash, do run `tcsh` to swich to tcsh. Then do the following.
```
source $CMSSW_BASE/src/SusyAnaTools/Tools/setup.csh
getTaggerCfg.sh -t MVAAK8_Tight_v1.2.1
```
If you are using bash run `bash` to return to bash.

The flag `-t` is for tag. There is a `-o` flag for override for overriding an existing version of the TopTagger in your area.

Check the TopTagger.cfg soft link with `ls -l TopTagger.cfg`. You should see the same version of the TopTagger that you checked out.

```
cmslpc23.fnal.gov Tools $ ls -l TopTagger.cfg
lrwxrwxrwx 1 caleb us_cms 119 Jun 21 18:52 TopTagger.cfg -> /uscms/home/caleb/nobackup/SusyAnalysis/CMSSW_9_4_4/src/ZInvisible/Tools/TopTaggerCfg-MVAAK8_Tight_v1.2.1/TopTagger.cfg
```

Now compile.
```
cd ../Zinvisible/Tools
mkdir obj
make -j8
```

Copy some files from Caleb.


Now try running makePlots.
```
./makePlots -D ZJetsToNuNu -E 1000
```


Now you are ready to go!


