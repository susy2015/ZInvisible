# ZInvisible


### Recipe

Checkout a CMSSW release, 7_4_8 seems to work. 
```
cmsrel CMSSW_7_4_15
cd CMSSW_7_4_15/src
cmsenv
```

Checkout following repositories

```
git clone git@github.com:susy2015/ZInvisible.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

```

For the Zinvible estimation, we are currently using the master branch of the SusyAnaTools repository, and the master branch of the ZInvisible repository. 

Now compile the relevant pieces
```
cd recipeAUX
scram b

cd ../Zinvisible/Tools
mkdir obj
make
```

Now you are ready to go!


