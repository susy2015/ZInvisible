# ZInvisible


### Recipe

Checkout a CMSSW release, 7_4_8 seems to work. 
```
cmsrel CMSSW_7_4_8
cd CMSSW_7_4_8/src
cmsenv
```

Checkout following repositories

```
git clone git@github.com:susy2015/ZInvisible.git
git clone git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:susy2015/recipeAUX.git

```

For the Zinvible estimation, we are currently using this branch of the SusyAnaTools repository:
```
cd SusyAnaTools
git checkout newPlotterDevel
```

Now compile the relevant pieces
```
cd recipeAUX
scram b

cd ../Zinvisible/Tools
mkdir obj
make
```

Now you are ready to go!


