#!/usr/bin/env python

import optparse
import subprocess
import datetime
import subprocess
import sys
import os
from os import system, environ
import json

sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path
from samples import SampleCollection

# https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763
import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# submit()
def submit(datasets, era, numfile=5, noSubmit=False, verbose=False, dataCollections=False, dataCollectionslong=False, systematic_string=''):
    print "# ---------- Submitting condor jobs for {0} ---------- #".format(era)
    # split string into list
    if datasets:
        datasets = datasets.split(',')
    else:
        datasets = []
        print "ERROR: datasets not found"
    # sample config files
    eras = ["2016", "2017", "2017_BE", "2017_F", "2018", "2018_PreHEM", "2018_PostHEM"]
    if era not in eras:
        print "Please use -y to enter era (2016, 2017, 2017_BE, 2017_F, 2018, 2018_PreHEM, or 2018_PostHEM)."
        exit(1)
    # era may include period
    # year does not include period
    year = era[0:4]
    sampleSetsFile = "sampleSets_PostProcessed_" + year + ".cfg"
    sampleCollectionsFile = "sampleCollections_" + year + ".cfg"

    filestoTransferGMP = [environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/makePlots",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6/",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6/",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6/",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleSets_PostProcessed_2016.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleSets_PostProcessed_2017.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleSets_PostProcessed_2018.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleCollections_2016.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleCollections_2017.cfg",
                          environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/sampleCollections_2018.cfg",
                         ]

    submitFileGMP = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh,gmp.tar.gz,$ENV(CMSSW_VERSION).tar.gz
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)
+maxWallTime = 2880

"""

    submitFile = ""
    exeName = ""

    # load samples
    sc = SampleCollection("../" + sampleSetsFile, "../" + sampleCollectionsFile)

    # make directory for condor submission
    now = datetime.datetime.now()
    dirName = "DataMC_%s_%s_submission_%s" % (era, systematic_string, now.strftime("%Y-%m-%d_%H-%M-%S"))
    os.system("mkdir %s"%dirName)
    #os.chdir(dirName)
    # use class to cd into directory and then return to previous directory
    with cd(dirName):

        def makeExeAndFriendsTarrball(filestoTransfer, fname):
            if not dataCollections and not dataCollectionslong:
                #WORLDSWORSTSOLUTIONTOAPROBLEM
                system("mkdir -p WORLDSWORSTSOLUTIONTOAPROBLEM")
                for fn in filestoTransfer:
                    system("cd WORLDSWORSTSOLUTIONTOAPROBLEM; ln -s %s" % fn)

                tarallinputs = "tar czf %s.tar.gz WORLDSWORSTSOLUTIONTOAPROBLEM --dereference" % fname
                if verbose:
                    # Use v option if verbose is set
                    tarallinputs = "tar czvf %s.tar.gz WORLDSWORSTSOLUTIONTOAPROBLEM --dereference" % fname
                    print "Create tarball {0}.tag.gz".format(fname)
                    print tarallinputs
                system(tarallinputs)
                system("rm -r WORLDSWORSTSOLUTIONTOAPROBLEM")


        if not dataCollections and not dataCollectionslong:
            if verbose:
                print "Create tarball ${CMSSW_VERSION}.tar.gz"
            system("tar --exclude-caches-all --exclude-vcs -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp")

        # makeExeAndFriendsTarrball() is necessary now to apply WORLDSWORSTSOLUTIONTOAPROBLEM

        exeName = "makePlots"
        submitFile = submitFileGMP
        makeExeAndFriendsTarrball(filestoTransferGMP, "gmp")

        nFilesPerJob = numfile

        fileParts = [submitFile]

        if dataCollections or dataCollectionslong:
            scl = sc.sampleCollectionList()
            for sampleCollection in scl:
                sl = sc.sampleList(sampleCollection)
                print sampleCollection
                if dataCollectionslong:
                    sys.stdout.write("\t")
                    for sample in sl:
                        sys.stdout.write("%s  "%sample[1])
                    print ""
                    print ""
            exit(1)

        debug = False
        total_files = 0
        total_events = 0
        jobMap = {}
        for ds in datasets:
            sample_files = 0
            sample_events = 0
            ds = ds.strip()

            # s: file, n:name, e:nEvts
            for s, n, e in sc.sampleList(ds):
                if debug:
                    print "s={0}, n={1}, e={2}".format(s, n, e)
                n_files = sum(1 for line in open(s))
                sample_files += n_files
                sample_events += e
                total_files += n_files
                total_events += e
                jobMap[n] =  n_files
                if verbose:
                    print "\t{0}: n_files={1}, n_events={2}".format(n, n_files, e)
                try:
                    f = open(s)
                except IOError:
                    fShort = s.split("/")[-1]
                    if(os.path.isfile(fShort)):
                        os.remove(fShort)
                    system("xrdcp root://cmseos.fnal.gov/$(echo %s | sed 's|/eos/uscms||') ."%s)
                    if verbose:
                        print "fShort = {0}".format(fShort)
                    f = open(fShort)
                if not f == None:
                    count = 0
                    for l in f:
                        if '.root' in l and not 'failed' in l:
                            count = count + 1
                    for startFileNum in xrange(0, count, nFilesPerJob):
                        fileParts.append("Arguments = %s $ENV(CMSSW_VERSION) %i %i %s %s %s\n"%(n, nFilesPerJob, startFileNum, systematic_string, s, era))
                        fileParts.append("Output = logs/$(Cluster)_$(Process)_%s_%s_%i.stdout\n"%(exeName, n, startFileNum))
                        fileParts.append("Error = logs/$(Cluster)_$(Process)_%s_%s_%i.stderr\n"%(exeName, n, startFileNum))
                        fileParts.append("Log = logs/$(Cluster)_$(Process)_%s_%s_%i.log\n"%(exeName, n, startFileNum))
                        fileParts.append("Queue\n\n")

                    f.close()
            # number of files and events in sample
            if verbose:
                print "\t{0}: sample_files={1}, sample_events={2}".format(ds, sample_files, sample_events)

        # total number of files and events
        if verbose:
            print "# --- Totals --- #"
            print "total_files={0}, total_events={1}".format(total_files, total_events)

        with open("nJobs.json", "w") as outfile:
            json.dump(jobMap, outfile, sort_keys=True, indent=4)

        fout = open("condor_submit.txt", "w")
        fout.write(''.join(fileParts))
        fout.close()

        if not noSubmit:
            system('mkdir -p logs')
            system("echo 'condor_submit condor_submit.txt'")
            system('condor_submit condor_submit.txt')

    print "Condor submission directory: {0}".format(dirName)
    return dirName

# main()
def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")

    parser.add_option ('-n',  dest='numfile',               type='int',          default = 5,     help="number of files per job (default: 5)")
    parser.add_option ('-d',  dest='datasets',              type='string',       default = '',    help="List of datasets 'ZJetsToNuNu_2016,GJets_2016,DYJetsToLL_2016'")
    parser.add_option ('-l',  dest='dataCollections',       action='store_true', default = False, help="List all datacollections")
    parser.add_option ('-L',  dest='dataCollectionslong',   action='store_true', default = False, help="List all datacollections and sub collections")
    parser.add_option ('-y',  dest='era',                   type='string',       default = None,  help="Year of Data or MC to analyze.")
    parser.add_option ('-c',  dest='noSubmit',              action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('-v',  dest='verbose',               action='store_true', default = False, help="Print more things.")
    parser.add_option ('-S',  dest='systematic',            type='string',       default = 'Base',help="String denoting which systematic to run: {empty string, JESUp, JESDown, METUnClustUp, METUnClustDown}")

    options, args = parser.parse_args()

    numfile = options.numfile
    datasets = options.datasets
    dataCollections = options.dataCollections
    dataCollectionslong = options.dataCollectionslong
    era = options.era
    noSubmit = options.noSubmit
    verbose = options.verbose
    systematic = options.systematic

    submit(datasets, era, numfile, noSubmit, verbose, dataCollections, dataCollectionslong, systematic)

if __name__ == "__main__":
    main()

