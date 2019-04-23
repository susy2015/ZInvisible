#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import subprocess
import sys
import os
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse 
import subprocess


# main()
def main():
    
    parser = optparse.OptionParser("usage: %prog [options]\n")
    
    parser.add_option ('-n',  dest='numfile',               type='int',          default = 5,     help="number of files per job")
    parser.add_option ('-d',  dest='datasets',              type='string',       default = '',    help="List of datasets 'ZJetsToNuNu,GJets,DYJetsToLL'")
    parser.add_option ('-l',  dest='dataCollections',       action='store_true', default = False, help="List all datacollections")
    parser.add_option ('-L',  dest='dataCollectionslong',   action='store_true', default = False, help="List all datacollections and sub collections")
    parser.add_option ('-r',  dest='refLumi',               type='string',       default = None,  help="Data collection to define lumi (uses default lumi if no reference data collection is defined)")
    parser.add_option ('-y',  dest='year',                  type='string',       default = None,  help="Year data or MC to analyze.")
    parser.add_option ('-c',  dest='noSubmit',              action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('-e',  dest='goMakeEff',             action='store_true', default = False, help="Run calcEff instead of makePlots.")
    parser.add_option ('-p',  dest='goMakeEffPhoton',       action='store_true', default = False, help="Run calcEffPhoton instead of makePlots.")
    parser.add_option ('-b',  dest='goMakeBeff',            action='store_true', default = False, help="Run beffCalc instead of makePlots.")
    parser.add_option ('-s',  dest='goMakeSigEff',          action='store_true', default = False, help="Run makeSignalHistograms instead of makePlots.")
    parser.add_option ('-t',  dest='goMakeTopPlots',        action='store_true', default = False, help="Run makeTopPlots instead of makePlots.")
    parser.add_option ('-m',  dest='goTTPlots',             action='store_true', default = False, help="Run TTPlots instead of makePlots.")
    
    options, args = parser.parse_args()
    year = options.year
    
    datasets = []
    if options.datasets:
        datasets = options.datasets.split(',')
    else:
        print "Please use -d to specify one or more datasets (separate with commas)."
        exit(1)

    # sample config files
    if year != "2016" and year != "2017" and year != "2018":
        print "Please use -y to enter year (2016, 2017, or 2018)."
        exit(1)

    sampleSetsFile = "sampleSets_PostProcessed_" + year + ".cfg"
    sampleCollectionsFile = "sampleCollections_" + year + ".cfg"
    
    # TopTagger.cfg
    mvaFileName = ""
    with file(environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg") as meowttcfgFile:
        for line in meowttcfgFile:
            line = line.split("#")[0]
            if "modelFile" in line:
                mvaFileName = line.split("=")[1].strip().strip("\"")
                break
    
    # TopTagger_Deep.cfg
    mvaFileName_Deep = ""
    with file(environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_Deep.cfg") as meowttcfgFile:
        for line in meowttcfgFile:
            line = line.split("#")[0]
            if "modelFile" in line:
                mvaFileName_Deep = line.split("=")[1].strip().strip("\"")
                break
    
    # TopTagger_DeepCombined.cfg
    mvaFileName_DeepCombined = ""
    with file(environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_DeepCombined.cfg") as meowttcfgFile:
        for line in meowttcfgFile:
            line = line.split("#")[0]
            if "modelFile" in line:
                mvaFileName_DeepCombined = line.split("=")[1].strip().strip("\"")
                break
    
    print "mvaFileName = {0}".format(mvaFileName)
    print "mvaFileName_Deep = {0}".format(mvaFileName_Deep)
    print "mvaFileName_DeepCombined = {0}".format(mvaFileName_DeepCombined)
    
    
    #here I hack in the tarball for GMP, this needs to be generalized to the other options 
    
    filestoTransferGMP = [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/makePlots", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/lepEffHists.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/effhists_GJets.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/njetWgtHists.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/dataMCweights.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/dataMCreweight.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/dataMCreweight_allJets.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/CSVv2_Moriond17_B_H.csv", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/puppiCorr.root",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_Deep.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_DeepCombined.cfg",
                          environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_core.so.3.1",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_ml.so.3.1",
                          #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/allINone_bTagEff.root", 
                          #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/ISRWeights.root", 
                          #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/allINone_ISRJets.root", 
                          #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/PileupHistograms_0121_69p2mb_pm4p6.root",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/allINone_bTagEff.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/ISRWeights.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/allINone_ISRJets.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/PileupHistograms_0121_69p2mb_pm4p6.root"
                         ]
    
    if mvaFileName:                 filestoTransferGMP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName}]
    if mvaFileName_Deep:            filestoTransferGMP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName_Deep}]
    if mvaFileName_DeepCombined:    filestoTransferGMP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName_DeepCombined}]
    if sampleSetsFile:              filestoTransferGMP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleSetsFile]
    if sampleCollectionsFile:       filestoTransferGMP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleCollectionsFile]
    
    filestoTransferGMEP = [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/calcEffPhoton", 
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/CSVv2_Moriond17_B_H.csv", 
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/puppiCorr.root",
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg",
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_Deep.cfg",
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_DeepCombined.cfg",
                           environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                           "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_core.so.3.1",
                           "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_ml.so.3.1",
                           #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/allINone_bTagEff.root", 
                           #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/ISRWeights.root", 
                           #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/allINone_ISRJets.root", 
                           #environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/PileupHistograms_0121_69p2mb_pm4p6.root",
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/allINone_bTagEff.root", 
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/ISRWeights.root", 
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/allINone_ISRJets.root", 
                           environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/PileupHistograms_0121_69p2mb_pm4p6.root"
                          ]
    
    if mvaFileName:                 filestoTransferGMEP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName}]
    if mvaFileName_Deep:            filestoTransferGMEP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName_Deep}]
    if mvaFileName_DeepCombined:    filestoTransferGMEP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName_DeepCombined}]
    if sampleSetsFile:              filestoTransferGMEP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleSetsFile]
    if sampleCollectionsFile:       filestoTransferGMEP += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleCollectionsFile]
    
    #go make plots!
    submitFileGMP = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/gmp.tar.gz,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/$ENV(CMSSW_VERSION).tar.gz
    Output = logs/makePlots_$(Process).stdout
    Error = logs/makePlots_$(Process).stderr
    Log = logs/makePlots_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    +maxWallTime = 2880
    
    """
    
    #Here is the configuration for the Data/MC validation of the TopTagger 
    filestoTransferTT  = [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/simpleAnalyzer",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleCollectionsFile,
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/" + sampleSetsFile,
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/Legacy_TopTagger.cfg",
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_AllComb.cfg",
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg_trijetOnly.cfg",
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg-MVAAK8_Tight_v1.2.1_dijetOnly.cfg",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_core.so.3.1",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_ml.so.3.1",
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"TrainingOutput_dR20_pt30_depth12_500tree_2017_Feb16.model"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"weights-t2tt850-sm-baseline-nodphi-nomtb-hqu-08112016.xml"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"sdWTag_ttbarTraining_v0.xml"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"tfModel_frozen.pb"},
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/puppiCorr.root",
                          environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/allINone_bTagEff.root", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/ISRWeights.root", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/allINone_ISRJets.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/CSVv2_Moriond17_B_H.csv", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/PileupHistograms_0121_69p2mb_pm4p6.root", 
                          ]
    
    if mvaFileName: filestoTransferTT += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName}]
    
    #go make TTopTagger plots!
    submitFileTT = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goTTplots.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goTTplots.sh,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/TT.tar.gz,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/$ENV(CMSSW_VERSION).tar.gz
    Output = logs/TT_$(Process).stdout
    Error = logs/TT_$(Process).stderr
    Log = logs/TT_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    +maxWallTime = 2880
    
    """
    
    #go make top plots!
    filestoTransferGTP = [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/makeTopPlots",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/Legacy_TopTagger.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger_AllComb.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg_trijetOnly.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg-MVAAK8_Tight_v1.2.1_dijetOnly.cfg",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_core.so.3.1",
                          "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_ml.so.3.1",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName},
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"TrainingOutput_dR20_pt30_depth12_500tree_2017_Feb16.model"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"weights-t2tt850-sm-baseline-nodphi-nomtb-hqu-08112016.xml"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"sdWTag_ttbarTraining_v0.xml"},
                          #environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":"tfModel_frozen.pb"},
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/puppiCorr.root",
                          environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/allINone_bTagEff.root", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/ISRWeights.root", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/ISR_Root_Files/allINone_ISRJets.root", 
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/CSVv2_Moriond17_B_H.csv", 
                          environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/data/PileupHistograms_0121_69p2mb_pm4p6.root", 
                          ]
    
    if mvaFileName: filestoTransferTT += [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/%(trainingFile)s"%{"trainingFile":mvaFileName}]
    
    #print filestoTransferGTP
    
    submitFileGTP = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeTopPlots.sh
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/gtp.tar.gz,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/$ENV(CMSSW_VERSION).tar.gz 
    Output = logs/makePlots_$(Process).stdout
    Error = logs/makePlots_$(Process).stderr
    Log = logs/makePlots_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    +maxWallTime = 2880
    
    """
    
    #go make lepton efficiency
    submitFileGME = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEff.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/calcEff, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEff.sh, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/zRes.root
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    
    """
    
    #go make photon efficiency
    submitFileGMEP = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEffPhoton.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEffPhoton.sh,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/gmep.tar.gz,$ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/$ENV(CMSSW_VERSION).tar.gz
    Output = logs/calcEffPhoton_$(Process).stdout
    Error = logs/calcEffPhoton_$(Process).stderr
    Log = logs/calcEffPhoton_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    +maxWallTime = 2880
    
    """
    
    #go B Efficiency calc
    submitFileGBE = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeBeff.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/beffCalc, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeBeff.sh
    Output = logs/makePlots_$(Process).stdout
    Error = logs/makePlots_$(Process).stderr
    Log = logs/makePlots_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    
    """
    
    #go Signal Efficiency
    submitFileGSE = """universe = vanilla
    Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeSigEff.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/makeSignalHistograms, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeSigEff.sh
    Output = logs/makePlots_$(Process).stdout
    Error = logs/makePlots_$(Process).stderr
    Log = logs/makePlots_$(Process).log
    notify_user = ${LOGNAME}@FNAL.GOV
    x509userproxy = $ENV(X509_USER_PROXY)
    
    """
    
    
    submitFile = ""
    exeName = ""
    
    def makeExeAndFriendsTarrball(filestoTransfer, fname):
        if not options.dataCollections and not options.dataCollectionslong:
            #WORLDSWORSESOLUTIONTOAPROBLEM
            system("mkdir -p WORLDSWORSESOLUTIONTOAPROBLEM")
            for fn in filestoTransfer:
                #print "fn = {0}".format(fn)
                system("cd WORLDSWORSESOLUTIONTOAPROBLEM; ln -s %s" % fn)
            
            print "Create tarball {0}.tag.gz".format(fname)
            tarallinputs = "tar czvf %s.tar.gz WORLDSWORSESOLUTIONTOAPROBLEM --dereference" % fname
            print tarallinputs
            system(tarallinputs)
            system("rm -r WORLDSWORSESOLUTIONTOAPROBLEM")
    
    
    if not options.dataCollections and not options.dataCollectionslong:
        print "Create tarball ${CMSSW_VERSION}.tar.gz"
        system("tar --exclude-caches-all --exclude-vcs -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp")
    
    # makeExeAndFriendsTarrball() is necessary now to apply WORLDSWORSESOLUTIONTOAPROBLEM 
    
    if options.goMakeEff:
        exeName = "calcEff"
        submitFile = submitFileGME
    elif options.goMakeEffPhoton:
        exeName = "calcEffPhoton"
        submitFile = submitFileGMEP
        makeExeAndFriendsTarrball(filestoTransferGMEP, "gmep")
    elif options.goMakeBeff:
        exeName = "beffCalc"
        submitFile = submitFileGBE
    elif options.goMakeSigEff:
        exeName = "makeSignalHistograms"
        submitFile = submitFileGSE
    elif options.goMakeTopPlots:
        exeName = "makeTopPlots"
        submitFile = submitFileGTP
        makeExeAndFriendsTarrball(filestoTransferGTP, "gtp")
    elif options.goTTPlots:
        exeName = "simpleAnalyzer"
        submitFile = submitFileTT
        makeExeAndFriendsTarrball(filestoTransferTT, "TT")
    else:
        exeName = "makePlots"
        submitFile = submitFileGMP
        makeExeAndFriendsTarrball(filestoTransferGMP, "gmp")
    
    nFilesPerJob = options.numfile
    
    fileParts = [submitFile]
    sc = SampleCollection("../" + sampleSetsFile, "../" + sampleCollectionsFile)
    
    if options.dataCollections or options.dataCollectionslong:
        scl = sc.sampleCollectionList()
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            print sampleCollection
            if options.dataCollectionslong:
                sys.stdout.write("\t")
                for sample in sl:
                    sys.stdout.write("%s  "%sample[1])
                print ""
                print ""
        exit(1)
    
    
    lumis = sc.sampleCollectionLumiList()
    lumi = sc.getFixedLumi()
    if options.refLumi != None:
        lumi = lumis[options.refLumi]
        print "Normalizing to %s pb-1" % (lumi)
    
    for ds in datasets:
        total_files = 0
        total_events = 0
        ds = ds.strip()
        #print ds
        print "# --- {0} --- #".format(ds)
        
        # s: file, n:name, e:nEvts
        for s, n, e in sc.sampleList(ds):
            #print "\t{0} {1} {2}".format(s, n, e)
            n_files = sum(1 for line in open(s))
            total_files += n_files
            total_events += e
            #n_files = subprocess.check_output(["grep", "root " + s +" | wc -l"])  
            #print "\t%s"%n
            print "\t{0}: n_files={1}, n_events={2}".format(n, n_files, e)
            try:
                f = open(s)
            except IOError:
                fShort = s.split("/")[-1]
                if(os.path.isfile(fShort)):
                    os.remove(fShort)
                system("xrdcp root://cmseos.fnal.gov/$(echo %s | sed 's|/eos/uscms||') ."%s)
                print "fShort = {0}".format(fShort)
                f = open(fShort)
            if not f == None:
                count = 0
                for l in f:
                    if '.root' in l and not 'failed' in l:
                        count = count + 1
                for startFileNum in xrange(0, count, nFilesPerJob):
                    fileParts.append("Arguments = %s $ENV(CMSSW_VERSION) %i %i %f %s %s\n"%(n, nFilesPerJob, startFileNum, lumi, s, year))
                    fileParts.append("Output = logs/%s_%s_%i.stdout\n"%(exeName, n, startFileNum))
                    fileParts.append("Error = logs/%s_%s_%i.stderr\n"%(exeName, n, startFileNum))
                    fileParts.append("Log = logs/%s_%s_%i.log\n"%(exeName, n, startFileNum))
                    fileParts.append("Queue\n\n")
    
                f.close()
        print "\t{0}: total_files={1}, total_events={2}".format(ds, total_files, total_events)
        #print "{0}\n{1}".format(ds, total_files)
    
    fout = open("condor_submit.txt", "w")
    fout.write(''.join(fileParts))
    fout.close()
    
    if not options.noSubmit: 
        system('mkdir -p logs')
        system("echo 'condor_submit condor_submit.txt'")
        system('condor_submit condor_submit.txt')

if __name__ == "__main__":
    main()

