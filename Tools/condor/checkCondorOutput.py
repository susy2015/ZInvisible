# checkCondorOutput.py

import argparse
import glob
import os
import re
import json

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory",      "-d", default="",    help="directory for condor jobs")
    parser.add_argument("--files_per_job",  "-n", default=-1,    help="number of files per job")
    
    options = parser.parse_args()
    directory = options.directory
    files_per_job = int(options.files_per_job)
    jobsMap = {}

    if not directory:
        print "Please enter a directory using the -d option."
        return
    output_directory = directory + "/output/"
    json_file = directory + "/nJobs.json"
    if not os.path.isdir(output_directory):
        print "ERROR: The directory {0} does not exist.".format(output_directory)
        return
    if not os.path.isfile(json_file):
        print "ERROR: The file {0} does not exist.".format(json_file)
        return
    if files_per_job < 1:
        print "Please enter number of files per job (greater than 0) using the -n option."
        return
    
    with open(json_file, "r") as infile:
        jobsMap = json.load(infile)
 
    # example file names
    # DataMC_2016_submission_2019-05-16_10-06-59/output/histoutput_ZJetsToNuNu_HT_100to200_2016_0.root
    # DataMC_2016_submission_2019-05-16_10-06-59/output/histoutput_TTbarDiLep_2016_0.root
    # DataMC_2016_submission_2019-05-16_10-06-59/output/histoutput_Data_SingleMuon_2016_PeriodH_118.root
    
    ####################################
    # WARNING: Do not break the regex  #
    ####################################
    # regex to get job number
    regex = re.compile(".*histoutput_(.*)_([0-9]+).root")

    root_files = glob.glob(output_directory + "*.root")
    samples = {}
    n_files = 0
    nSubmittedTotal = 0
    nReturnedTotal = 0
    for f in root_files:
        match = regex.match(f)
        sample = match.group(1)
        job_number = int(match.group(2))
        if sample not in samples:
            samples[sample] = [job_number]
        else:
            samples[sample] += [job_number]
        #print "{0} {1}".format(match.group(1), match.group(2))
        n_files += 1
    
    print "--- Loop over submitted samples ---"
    # loop over submitted samples
    for sample in jobsMap:
        nJobs = jobsMap[sample]["nJobs"]
        print "nJobs={0}".format(nJobs)
        nReturned = -1
        try:
            nReturned = len(samples[sample])
        except:
            print "ERROR: The sample {0} was not found in {1}".format(sample, output_directory)
        if nReturned >= 0:
            # nSubmitted: divide nJobs by files_per_job and round up: e.g. for 2 files per job, number of jobs is 4/2 = 2, 5/2 = 3
            nSubmitted = int(round(nJobs/files_per_job))
            diff = nSubmitted - nReturned 
            nSubmittedTotal += nSubmitted
            nReturnedTotal += nReturned
            print "{0}: {1} - {2} = {3}".format(sample, nSubmitted, nReturned, diff)
        else:
            print "ERROR: No jobs found for for {0}".format(sample)

    
    print "--- Loop over returned samples ---"
    # loop over returned samples
    for sample in samples:
        samples[sample].sort()
        nReturned = len(samples[sample]) 
        nJobs = -1
        try:
            nJobs = jobsMap[sample]["nJobs"]
        except:
            print "ERROR: The sample {0} was not found in {1}".format(sample, json_file)
        
        if nJobs >= 0:
            # nSubmitted: divide nJobs by files_per_job and round up: e.g. for 2 files per job, number of jobs is 4/2 = 2, 5/2 = 3
            nSubmitted = int(round(nJobs/files_per_job))
            diff = nSubmitted - nReturned 
            print "{0}: {1} - {2} = {3}".format(sample, nSubmitted, nReturned, diff)
        else:
            print "ERROR: No jobs found for for {0}".format(sample)
        
        #print "{0} {1}".format(sample, samples[sample])
        #print "{0} {1}".format(sample, nReturned)
        
        expected = (i for i in xrange(0, max(samples[sample]), files_per_job))
        s_expected = set(expected) 
        s_returned = set(samples[sample])
        s_diff = s_expected.difference(s_returned)
        if s_diff:
            print "WARNING: This sample has missing jobs: {0} jobs {1}".format(sample, s_diff)

    print "n_files = {0}".format(n_files)
    print "nSubmittedTotal = {0}".format(nSubmittedTotal)
    print "nReturnedTotal = {0}".format(nReturnedTotal)

if __name__ == "__main__":
    main()


