# Generate systematics json file from csv

import csv
import json

with open("systematics.csv", "r") as csvFile:
    csvObject = csv.reader(csvFile, delimiter = ',')

    jsonDict = {} 

    for row in csvObject:
        # keys are in row starting with key
        if (row[0] == "key"):
            keys = row
        jsonDict[row[0]] = {}
        # loop over columns to match keys to values
        for i in xrange(1, len(keys)):
            jsonDict[row[0]][keys[i]] = row[i]

    with open("systematics.json", "w") as jsonFile:
        json.dump(jsonDict,jsonFile, sort_keys=True, indent=4, separators=(',', ': '))

