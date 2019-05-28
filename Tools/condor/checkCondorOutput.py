# checkCondorOutput.py

import argparse
import glob
import os
import re

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory",      "-d", default="",    help="directory for condor jobs")
    parser.add_argument("--files_per_job",  "-n", default=-1,    help="number of files per job")
    
    options = parser.parse_args()
    directory = options.directory
    files_per_job = int(options.files_per_job)

    if not directory:
        print "Please enter a directory using the -d option."
        return
    output_directory = directory + "/output/"
    if not os.path.isdir(output_directory):
        print "ERROR: The directory {0} does not exist.".format(output_directory)
        return
    if files_per_job < 1:
        print "Please enter number of files per job (greater than 0) using the -n option."
        return
 
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
    for sample in samples:
        samples[sample].sort()
        #print "{0} {1}".format(sample, samples[sample])
        expected = (i for i in xrange(0, max(samples[sample]), files_per_job))
        s_expected = set(expected) 
        s_returned = set(samples[sample])
        diff = s_expected.difference(s_returned)
        if diff:
            print "WARNING: This sample has missing jobs: {0} {1}".format(sample, diff)

    print "n_files = {0}".format(n_files)

if __name__ == "__main__":
    main()


