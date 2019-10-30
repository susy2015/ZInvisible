# --------------------- #
#   multiSubmit.py      #
#   Caleb J. Smith      #
#   August 20, 2019     #
# --------------------- #

# ------------------------------------ #
#          _____ __  __  _____         #
#         / ____|  \/  |/ ____|        #   
#        | |    | \  / | (___          # 
#        | |    | |\/| |\___ \         # 
#        | |____| |  | |____) |        # 
#         \_____|_|  |_|_____/         # 
#                                      # 
# ------------------------------------ #

# 1. submit multiple condor jobs
# 2. create json file stories directories

import json
import datetime
from condorSubmit import submit

def main():
    now = datetime.datetime.now()
    jsonName = "../runs/submission_%s.json" % (now.strftime("%Y-%m-%d_%H-%M-%S"))
    
    # different datasets for different eras; note the order
    dataSetsList = []
    dataSetsList.append("Data_MET_{0},Data_SingleElectron_{0},Data_SingleMuon_{0},Data_SinglePhoton_{0},DYJetsToLL_{1},TTbarNoHad_{1},SingleTopZinv_{1},Rare_{1},TTZ_{1},Diboson_{1},GJets_{1},QCD_Photon_{1},WJetsToLNu_{1},TTbar_{1},tW_{1},ZJetsToNuNu_{1}")
    dataSetsList.append("Data_MET_{0},Data_EGamma_{0},Data_SingleMuon_{0},DYJetsToLL_{1},TTbarNoHad_{1},SingleTopZinv_{1},Rare_{1},TTZ_{1},Diboson_{1},GJets_{1},QCD_Photon_{1},WJetsToLNu_{1},TTbar_{1},tW_{1},ZJetsToNuNu_{1}")
    dataSetsMap = {
            "2016"          : 0,
            "2017_BE"       : 0,
            "2017_F"        : 0,
            "2018_PreHEM"   : 1,
            "2018_PostHEM"  : 1
    }
    eras = ["2016", "2017_BE", "2017_F", "2018_PreHEM", "2018_PostHEM"]
    #eras = ["2016"]
    dirMap = {}
    # submit jobs for each era 
    for era in eras:
        year = era[0:4] 
        datasets = dataSetsList[dataSetsMap[era]].format(era, year)
        refLumi = "Data_MET_{0}".format(era)
        dirMap[era] = submit(datasets, refLumi, era, numfile=5, noSubmit=False, verbose=False)
    
    # write directory map to json file 
    with open (jsonName, "w") as j:
        json.dump(dirMap, j, sort_keys=True, indent=4, separators=(',', ' : '))
    print "json file containing submission directories: {0}".format(jsonName)

if __name__ == "__main__":
    main()


