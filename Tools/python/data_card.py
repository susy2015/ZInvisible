# data_card.py

from tools import isclose

def writeLine(f, line):
    f.write(line + "\n")

def makeDataCard(BinObject, directory, era):
    
    print ">>> Making data card for {0}.".format(era)
    
    out_name = "{0}zinv_dataCard_{1}.txt".format(directory, era)
    # open file
    out_file = open(out_name, "w+")
    nChannels = len(BinObject.all_bins) 
    sample = "zinv"
    writeLine(out_file, "luminosity = {0}".format(-999))
    writeLine(out_file, "channels = {0}".format(nChannels))
    writeLine(out_file, "sample = {0}".format(sample))
    
    # search bin number
    channel = "channel = "
    # prediction
    rate = "rate = "
    # effective number of events
    cs_event = "cs_event = "
    # weighted average of the weights
    avg_weight = "avg_weight = "
    # prediction:                   p     = bin value
    # uncertainty:                  sigma = bin error
    # average weight:               avg_w = sigma^2 / p 
    # effective number of events:   N_eff = p / avg_w
    
    for bin_i in BinObject.all_bins:
        # bin_i starts from 0
        # for datacard, use b starting from 1
        b           = str(int(bin_i) + 1)
        pred        = float(BinObject.binValues[era][bin_i]["pred"])
        sigma       = float(BinObject.binValues[era][bin_i]["pred_error_propagated"])
        n_eff       = float(BinObject.binValues[era][bin_i]["n_eff"])
        avg_w       = float(BinObject.binValues[era][bin_i]["avg_w"])
        n_eff_final = float(BinObject.binValues[era][bin_i]["n_eff_final"])
        avg_w_final = float(BinObject.binValues[era][bin_i]["avg_w_final"])
        channel     += "bin{0} ".format(b)
        rate        += "{0} ".format(pred)
        cs_event    += "{0} ".format(n_eff_final)
        avg_weight  += "{0} ".format(avg_w_final)
        x = pred
        y = n_eff * avg_w
        if (not isclose(x, y)):
            print "ERROR: bin {0}, pred = {1} and Neff * avgW = {2} are not equal".format(b, x, y)
   
    writeLine(out_file, channel)
    writeLine(out_file, rate)
    writeLine(out_file, cs_event)
    writeLine(out_file, avg_weight)
    
    # close file
    out_file.close()




