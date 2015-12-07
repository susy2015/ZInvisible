from ROOT import *
from optparse import OptionParser
import array

def rebin1D(h, bins):
    """Rebin histo h to bins and recompute the errors."""
    new_h = TH1D("%s_rebin"%(h.GetName()), "%s_rebin"%(h.GetName()), 
                 len(bins)-1, array.array('d', bins))
    new_h.Sumw2()
    for i in xrange(h.GetNbinsX()):
        bin_to_fill = new_h.FindBin(h.GetBinCenter(i+1))
        new_h.SetBinContent(bin_to_fill, h.GetBinContent(i+1) + new_h.GetBinContent(bin_to_fill))
        new_h.SetBinError(bin_to_fill, TMath.Sqrt(h.GetBinError(i+1)*h.GetBinError(i+1) + new_h.GetBinError(bin_to_fill)*new_h.GetBinError(bin_to_fill) ) )    
    return new_h

def add(hs):
    """Add all histograms in list hs and return output."""
    for h in hs[1:]:
        hs[0].Add(h)
    return hs[0]

def makeRatio(f, hname1, hnames2, bins, newname):
    """Make the data/stack ratio for the new bins."""
    # Get all histos
    h1 = f.Get(hname1)
    h2s = [f.Get(hname2) for hname2 in hnames2]
    # Add stack together
    h2 = add(h2s)
    # Rebin the histograms
    h1 = rebin1D(h1, bins)
    h2 = rebin1D(h2, bins)
    # Make the ratio
    h1.Divide(h2)
    h1.SetName(newname)
    h1.SetTitle(newname)
    return h1


if __name__ ==  "__main__":
    
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Grab histogram from FILE", metavar="FILE")

    (options, args) = parser.parse_args()

    # Get the file
    f = TFile.Open(options.filename)
    
    # How to rebin
    bins = [0,1,2,3,4,5,6,7,8,20]
    
    # Run over the relevant histograms
    cuts = ["muZinv", "muZinv_0b", "muZinv_g1b"]
    
    # Prepare writing to a file
    fout = TFile.Open("dataMCweights.root","RECREATE")
    fout.cd()

    # Data/MC plots in loose region for njets
    for cut in cuts:
        hname1 = "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDatadata"%(cut)
        hnames2 = ["DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDYstack"%(cut),
                   "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDY HT<100stack"%(cut),
                   "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}stack"%(cut),
                   "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvsingle topstack"%(cut),
                   "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}Zstack"%(cut),
                   "DataMC_SingleMuon_nj_%s_loose0cntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDibosonstack"%(cut)
                   ]
        newname = "DataMC_nj_%s_loose0"%(cut)
        newh = makeRatio(f, hname1, hnames2, bins, newname)

        newh.Write()

    fout.Close()
    f.Close()
