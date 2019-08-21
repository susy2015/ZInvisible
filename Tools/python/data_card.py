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
    
    for bin_i in BinObject.all_bins:
        # bin_i starts from 0
        # for datacard, use b starting from 1
        b = str(int(bin_i) + 1)
        pred    = float(BinObject.binValues[era][bin_i]["pred"])
        sigma   = float(BinObject.binValues[era][bin_i]["pred_error"])
        if pred == 0:
            print "ERROR: bin {0}, pred = {1}; seting avg weight to 1.0".format(b, pred)
            avg_w   = 1.0
        else:
            avg_w   = (sigma ** 2) / pred
        n_eff   = pred / avg_w
        n_eff_final = int(n_eff)
        if n_eff_final == 0:
            print "ERROR: bin {0}, n_eff_final = {1}; leaving avg weight unchanged".format(b, n_eff_final)
            avg_w_final = avg_w
        else:
            avg_w_final = pred / n_eff_final
        channel     += "bin{0} ".format(b)
        rate        += "{0} ".format(pred)
        cs_event    += "{0} ".format(n_eff_final)
        avg_weight  += "{0} ".format(avg_w_final)
        x = pred
        y = n_eff_final * avg_w_final
        if (not isclose(x, y)):
            print "ERROR: bin {0}, pred = {1} and Neff * avgW = {2}".format(b, x, y)
   
    writeLine(out_file, channel)
    writeLine(out_file, rate)
    writeLine(out_file, cs_event)
    writeLine(out_file, avg_weight)
    
    # close file
    out_file.close()




