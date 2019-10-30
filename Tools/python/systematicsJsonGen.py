# Generate systematics json file from csv

import csv
import json

with open("systematics.csv", "r") as csvFile:
    csvObject = csv.reader(csvFile, delimiter = ',')

    jsonDict = {} 

    for row in csvObject:
        if (row[0] == "key"):
            continue
        jsonDict[row[0]]= {
            "name"    : row[1],
            "nominal" : row[2],
            "up"      : row[3],
            "down"    : row[4]
        }

    with open("systematics.json", "w") as jsonFile:
        json.dump(jsonDict,jsonFile, sort_keys=True, indent=4, separators=(',', ': '))

