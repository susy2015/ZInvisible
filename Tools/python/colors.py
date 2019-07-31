# colors.py
import ROOT
import json

# xkcd color survey
# see https://xkcd.com/color/rgb.txt
# see https://xkcd.com/color/rgb/
# download text file using "wget https://xkcd.com/color/rgb.txt"
# generate json file from the text file by running main()

# print color, RGB, and index
def printColors():
    with open("rgb.json", "r") as input_file:
        data = json.load(input_file)
        for color in data:
            index = ROOT.TColor.GetColor(data[color])
            print "{0}: {1}: {2}".format(color, data[color], index) 

# return TColor index given xkcd color using RGB json file 
def getColorIndex(color):
    with open("rgb.json", "r") as input_file:
        data = json.load(input_file)
        index = ROOT.TColor.GetColor(data[color])
        #print "{0}: {1}: {2}".format(color, data[color], index) 
        return index

# make RGB json file from xkcd colors text file
def main():
    input_name = "rgb.txt"
    output_name = "rgb.json"
    color_map = {}
    with open(input_name, "r") as input_file:
        for line in input_file:
            line = line.strip()
            if "#" in line and line[0] != "#":
                #print line
                split = line.split("#")
                key = split[0].strip()
                value = "#" + split[1].strip()
                color_map[key] = value
    with open(output_name, "w") as output_file:
        json.dump(color_map, output_file, indent=4)

if __name__ == "__main__":
    # run main() to create json file
    main()
    
    # example of getting color index
    #getColorIndex("tomato")
    
    # print colors
    #printColors()

