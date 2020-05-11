# checkCondorOutput.py
import subprocess
import argparse
import glob
import os
import re
import json

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory",         "-d", default="",                                   help="directory for condor jobs")
    parser.add_argument("--check_output",      "-o", default = True,   action = "store_false",     help="check number of output root files, use flag to turn off")
    parser.add_argument("--verbose",           "-v", default = False,  action = "store_true",      help="verbose flag to print more things")
    
    options       = parser.parse_args()
    directory     = options.directory
    check_output  = options.check_output
    verbose       = options.verbose
    jobsMap = {}

    if not directory:
        print "Please enter a directory using the -d option."
        return
    output_directory = directory + "/output/"
    log_directory    = directory + "/logs/"
    json_file        = directory + "/nJobs.json"
    
    if output_directory[-1] != "/":
        output_directory += "/"
    if log_directory[-1] != "/":
        log_directory += "/"
    if not os.path.isdir(output_directory) and check_output:
        print "ERROR: The directory {0} does not exist.".format(output_directory)
        return
    if not os.path.isdir(log_directory):
        print "ERROR: The directory {0} does not exist.".format(log_directory)
        return
    if not os.path.isfile(json_file):
        print "ERROR: The file {0} does not exist.".format(json_file)
        return
    
    
    # check log files
    if check_output:
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
        
        if verbose:
            print "--- Loop over submitted samples ---"
        # loop over submitted samples
        for sample in jobsMap:
            nJobs = jobsMap[sample]["nJobs"]
            nReturned = -1
            try:
                nReturned = len(samples[sample])
            except:
                print "ERROR: The sample {0} was not found in {1}".format(sample, output_directory)
            if nReturned >= 0:
                nSubmitted = nJobs
                diff = nSubmitted - nReturned 
                nSubmittedTotal += nSubmitted
                nReturnedTotal += nReturned
                if verbose:
                    print "{0}: {1} - {2} = {3}".format(sample, nSubmitted, nReturned, diff)
            else:
                print "ERROR: No jobs found for for {0}".format(sample)

        
        if verbose:
            print "--- Loop over returned samples ---"
        # loop over returned samples
        for sample in samples:
            samples[sample].sort()
            nReturned = len(samples[sample]) 
            nFilesPerJob = jobsMap[sample]["nFilesPerJob"]
            nJobs = -1
            try:
                nJobs = jobsMap[sample]["nJobs"]
            except:
                print "ERROR: The sample {0} was not found in {1}".format(sample, json_file)
            
            if nJobs >= 0:
                nSubmitted = nJobs
                diff = nSubmitted - nReturned 
                if verbose:
                    print "{0}: {1} - {2} = {3}".format(sample, nSubmitted, nReturned, diff)
            else:
                print "ERROR: No jobs found for for {0}".format(sample)
            
            # use nFilesPerJob as step size, use max(samples[sample]) as max value
            expected = (i for i in xrange(0, max(samples[sample]) + nFilesPerJob, nFilesPerJob))
            s_expected = set(expected) 
            s_returned = set(samples[sample])
            s_diff = s_expected.difference(s_returned)
            if s_diff:
                print "WARNING: This sample has missing jobs: {0} jobs {1}".format(sample, s_diff)


        print "--- Totals ---"
        #print "nFiles = {0}".format(n_files)
        print "nSubmittedJobs = {0}".format(nSubmittedTotal)
        print "nReturnedJobs = {0}".format(nReturnedTotal)
        print "nMissingJobs = {0}".format(nSubmittedTotal - nReturnedTotal)
    
    # check log files
    # look for these key words (ignoring case)
    print "--- Logs ---"
    key_words = ["error", "warn", "not found"]
    for word in key_words:
        process = subprocess.Popen("grep -i \"{0}\" {1}*".format(word, log_directory),
                                   shell=True,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        out_list = stdout.split("\n")     # get list for each line of output
        out_list = filter(None, out_list) # remove empty strings from list 
        
        #print out_list
        print " > number of occurances of '{0}' in logs: {1}".format(word, len(out_list))

if __name__ == "__main__":
    main()


