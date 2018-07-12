# ZInvisible


### Recipe

Checkout a CMSSW 9_X_X release of 9_3_3 or later (including 9_4_X, but not 10_X_X_).
```
cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src
cmsenv
```

Checkout following repositories

```
git clone git@github.com:susy2015/ZInvisible.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:susy2015/TopTagger.git
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

```

For the Zinvible estimation, we are currently using the calebGJets branch of the ZInvisible repository. We are using the master branch of the SusyAnaTools repository (though the recent change of double to float is not supported yet). We are using the master branch of the TopTagger repository.

Now compile the relevant pieces
```
cd recipeAUX
scram b

cd ../Zinvisible/Tools
mkdir obj
make
```

Now you are ready to go!


