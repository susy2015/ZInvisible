# removeOutput.py

# Remove output and log files from condor jobs

import os
import shutil

def main():
    tag = "DataMC"
    tuples = os.walk(".")
    for t in tuples:
        d = t[0]
        if d[2:8] == tag:
            #if "output" in d: 
            if "logs" in d or "output" in d: 
                print "Removing {0}".format(d)
                ## Try to remove tree; if failed show an error using try...except on screen
                try:
                        shutil.rmtree(d)
                except OSError as e:
                        print ("Error: %s - %s." % (e.filename, e.strerror))

if __name__ == "__main__":
    main()


