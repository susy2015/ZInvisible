from optparse import OptionParser

import array, sys
# Check if we want the help, if not import ROOT (otherwise ROOT overwrites the help)
if '-h' not in sys.argv and '--help' not in sys.argv:
    from ROOT import TH1D, TMath, TFile, TCanvas, TF1

from math import sqrt
import math

############################
##  Some utilities first  ##
############################

c = TCanvas("c1", "c1", 800, 800)

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

def makeRatio(h1, h2s, bins=None, newname="test"):
    """Make the data/stack ratio for the new bins."""
    # Add stack together
    h2 = add(h2s)
    if bins != None:
        # Rebin the histograms
        h1 = rebin1D(h1, bins)
        h2 = rebin1D(h2, bins)
    # Make the ratio
    h1.Divide(h2)
    h1.SetName(newname)
    h1.SetTitle(newname)
    return h1

def reweight(h, hsf):
    """Reweight histogram h according to histogram hsf: h_i * hsf_i, for each bin i"""
    new_h = h.Clone()
    for i in xrange(new_h.GetNbinsX()):
        # Get the scale factor
        sfbin = hsf.FindBin(new_h.GetBinCenter(i+1))
        sf = hsf.GetBinContent(sfbin)
        sf_e = hsf.GetBinError(sfbin)
        if sf == 0: continue
        # Get the old bin info
        bc_old = new_h.GetBinContent(i+1)
        be_old = new_h.GetBinError(i+1)
        if bc_old == 0: continue
        # Get the new bin content
        bc = bc_old*sf
        # Get the new bin error, just use basic error propagation without correlations: [E(AB)/(AB)]^2 = [E(A)/A]^2 + [E(B)/B]^2
        be = bc * TMath.Sqrt( (be_old/bc_old)*(be_old/bc_old) + (sf_e/sf)*(sf_e/sf) )
        # Write the new info
        new_h.SetBinContent(i+1,bc)
        new_h.SetBinError(i+1,be)
    return new_h

def subtract(h, hlist):
    """Subtract all histograms in hlist from h"""
    new_h = h.Clone()
    new_h.Sumw2()
    for hl in hlist:
        new_h.Add(hl,-1.)
    return new_h

##############################
##  Main scale factor code  ##
##############################

##################
## Njet weights ##
##################
def njetWeights(filename):
    # Get the file
    f = TFile.Open(filename)

    # Prepare writing to a file
    fout = TFile.Open("dataMCweights.root","RECREATE")
    fout.cd()

    # How to rebin
    bins = [0,1,2,3,4,5,6,7,8,20]
    bins_TT = [0,1,2,3,4,5,6,7,20]

    # Run over the relevant histograms
    cuts_DY = ["muZinv", "muZinv_0b", "muZinv_g1b"]
    cuts_TT = ["elmuZinv", "elmuZinv_0b", "elmuZinv_g1b"]
    selection = "ht200_dphi"
    # histo names
    hname1 = "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDatadata"
    hnames2 = ["cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDYstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDY HT<100stack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}stack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvsingle topstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}Zstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDibosonstack"
               ]

    # dictionary to keep track of all the scale factors
    SFs = {}

    # First get TTbar reweighting, region is very pure, just take ratio of data/full stack and apply that to ttbar later
    for cut in cuts_TT:
        hname1_TT = hname1 % {"cut":cut, "selection":selection}
        hnames2_TT = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_TT)
        h2s = [f.Get(hname2_TT) for hname2_TT in hnames2_TT]
        newname = "DataMC_nj_%s_%s"%(cut,selection)
        # Make the ratio
        newh = makeRatio(h1, h2s, bins_TT, newname)
        SFs["TT_%s"%(cut)] = newh
        newh.Write()

    # Now get DY reweighting:
    # Procedure: 1. Apply TTbar reweighting to TTbar MC
    #            2. Subtract non-DY MC from data
    #            3. Make ratio of subtracted data and DY
    for cut in cuts_DY:
        hname1_DY = hname1 % {"cut":cut, "selection":selection}
        hnames2_DY = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_DY)
        h2s = [f.Get(hname2_DY) for hname2_DY in hnames2_DY]

        # apply weights to ttbar
        h2s[2] = reweight(h2s[2], SFs["TT_%s"%(cut.replace("mu","elmu"))])

        # subtract relevant histograms from data
        data_subtracted = subtract(h1, h2s[2:])

        newname = "DataMC_nj_%s_%s"%(cut,selection)
        newh = makeRatio(data_subtracted, h2s[:1], bins, newname)

        newh.Write()

    # Close files
    fout.Close()
    f.Close()

##################
## norm weights ##
##################
def normWeight(filename):
    # Get the file
    f = TFile.Open(filename)
    # Run over the relevant histograms
    cuts_DY = ["muZinv_0b"]
    selection = "blnotag"
    # histo names
    hname1 = "cntCSVSZinv/DataMC_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDatadata"
    hnames2 = ["cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDYstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDY HT<100stack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvt#bar{t}stack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvsingle topstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvt#bar{t}Zstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDibosonstack"
               ]

    # Procedure: 1. Grab the njet reweighted MC
    #            2. Subtract non-DY MC from data
    #            3. Make ratio of subtracted data and DY
    for cut in cuts_DY:
        hname1_DY = hname1 % {"cut":cut, "selection":selection}
        hnames2_DY = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_DY)
        h2s = [f.Get(hname2_DY) for hname2_DY in hnames2_DY]

        # subtract relevant histograms from data
        data_subtracted = subtract(h1, h2s[2:])

        newname = "DataMC_nb_%s_%s"%(cut,selection)
        newh = makeRatio(data_subtracted, h2s[:1], newname=newname)
        #newh = makeRatio(h1, h2s, newname=newname)

        print "Data/MC normalization scale factor in region %s_%s: %.3f +- %.3f" % (cuts_DY[0], selection, newh.GetBinContent(1), newh.GetBinError(1))

    f.Close()

def shapeSyst(filename):
    # Get the file
    f = TFile.Open(filename)
    fout = TFile.Open("syst_shape.root", "RECREATE")
    # Run over the relevant histograms
    # histo names
    hnameData = "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDatadata"
    hnames2 =  ["%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDYstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDY HT<100stack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)st#bar{t}stack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)ssingle topstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)st#bar{t}Zstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDibosonstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sRarestack"]

    varList = [#["met", "cleanMetPt",             [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 1500], "MET" ],
               #["mt2", "best_had_brJet_MT2Zinv", [0, 50, 100, 150, 200, 250, 300, 350, 400, 1500],      "M_{T2}" ],
               ["met", "cleanMetPt",             [0, 100, 200, 300, 600, 1500], "MET" ],
               ["mt2", "best_had_brJet_MT2Zinv", [0, 100, 200, 300, 400, 1500],      "M_{T2}" ],
               ["nt",  "nTopCandSortedCntZinv",  [0, 1, 2, 8 ],                                         "N(t)" ],
               ["nb",  "cntCSVSZinv",            [0, 1, 2, 3, 8 ],                                      "N(b)" ]]

    for var in varList:
        # Procedure: 1. Grab the njet reweighted MC
        #            2. Subtract non-DY MC from data
        #            3. Make ratio of subtracted data and DY
        # Get all histos
        hData = f.Get(hnameData%{"name":var[0], "var":var[1]})
        h2s   = [f.Get(hname2%{"name":var[0], "var":var[1]}) for hname2 in hnames2]
        
        # subtract relevant histograms from data
        data_subtracted = subtract(hData, h2s[2:])

        newname = "ShapeRatio_%s"%var[0]
        newh = makeRatio(data_subtracted, h2s[:1], newname=newname, bins=var[2])

        for i in xrange(1, newh.GetNbinsX() + 1):
            print newh.GetBinContent(i)

        fit = TF1("fit_%s"%var[0], "pol1")
        newh.SetLineWidth(2)
        newh.GetXaxis().SetTitle(var[3])
        newh.GetYaxis().SetTitle("Data/MC")
        newh.GetYaxis().SetTitleOffset(1.15)
        newh.SetStats(0)
        newh.SetTitle("")
        #if not "nt" in var[0] and not "nb" in var[0]:
        #    newh.Fit(fit, "")
        #else:
        #    newh.Draw()
        newh.Draw()
        #
        #    fit.Draw("same")
        #newh.DrawCopy()
        c.Print("shapeSyst_%s.png"%var[0])
        c.Print("shapeSyst_%s.pdf"%var[0])

        fout.cd()
        newh.Write()

    f.Close()
    fout.Close()

def systHarvest():
    # Get the file
    #f = TFile.Open(filename)
    #fout = TFile.Open("syst_shape.root", "RECREATE")
    # Run over the relevant histograms
    # histo names

    # Get shape central value uncertainty
    f = TFile("systematics.root")
    hShape_MET_Nom = f.Get("nSearchBin/systWgtMET_cleanMetPtnSearchBinnSearchBinNominalsingle")
    hShape_MET_Var = f.Get("nSearchBin/systWgtMET_cleanMetPtnSearchBinnSearchBinvariedsingle")
    hShape_MT2_Nom = f.Get("nSearchBin/systWgtMT2_best_had_brJet_MT2ZinvnSearchBinnSearchBinNominalsingle")
    hShape_MT2_Var = f.Get("nSearchBin/systWgtMT2_best_had_brJet_MT2ZinvnSearchBinnSearchBinvariedsingle")
    hShape_NT_Nom  = f.Get("nSearchBin/systWgtNT_nTopCandSortedCntZinvnSearchBinnSearchBinNominalsingle")
    hShape_NT_Var  = f.Get("nSearchBin/systWgtNT_nTopCandSortedCntZinvnSearchBinnSearchBinvariedsingle")
    hShape_NB_Nom  = f.Get("nSearchBin/systWgtNB_cntCSVSZinvnSearchBinnSearchBinNominalsingle")
    hShape_NB_Var  = f.Get("nSearchBin/systWgtNB_cntCSVSZinvnSearchBinnSearchBinvariedsingle")

    hShape_MET_ratio = hShape_MET_Var.Clone(hShape_MET_Nom.GetName()+"_ratio")
    hShape_MET_ratio.Divide(hShape_MET_Nom)

    hShape_MT2_ratio = hShape_MT2_Var.Clone(hShape_MT2_Nom.GetName()+"_ratio")
    hShape_MT2_ratio.Divide(hShape_MT2_Nom)

    hShape_NT_ratio = hShape_NT_Var.Clone(hShape_NT_Nom.GetName()+"_ratio")
    hShape_NT_ratio.Divide(hShape_NT_Nom)

    hShape_NB_ratio = hShape_NB_Var.Clone(hShape_NB_Nom.GetName()+"_ratio")
    hShape_NB_ratio.Divide(hShape_NB_Nom)

    hShape_final = hShape_MET_ratio.Clone("shape_central")

    for i in xrange(1, 46):
        uncertMET = hShape_MET_ratio.GetBinContent(i) - 1 if hShape_MET_ratio.GetBinContent(i) > 0 else 0
        uncertMT2 = hShape_MT2_ratio.GetBinContent(i) - 1 if hShape_MT2_ratio.GetBinContent(i) > 0 else 0
        uncertNT = hShape_NT_ratio.GetBinContent(i) - 1 if hShape_NT_ratio.GetBinContent(i) > 0 else 0
        uncertNB = hShape_NB_ratio.GetBinContent(i) - 1 if hShape_NB_ratio.GetBinContent(i) > 0 else 0
        hShape_final.SetBinContent(i, sqrt(uncertMET**2 + uncertMT2**2 + uncertNT**2 + uncertNB**2))

    fout = TFile("syst_all.root", "RECREATE")
    hShape_final.Write()

    # Get shape stats uncertainty
    f2 = TFile("syst_nJetWgt.root")
    hShapeStat = f2.Get("syst68Max").Clone("shape_stat")

    hAvgWgt = f2.Get("avgWgt").Clone("avgWgt")
    hNEff   = f2.Get("neff").Clone("neff")

    fout.cd()
    hShapeStat.Write()


    # Get PDF and scale uncertainties 
    f3 = TFile("syst_scalePDF.root")
    hScaleUp = f3.Get("nSearchBin_ratio_scale_up").Clone("scale_up")
    hScaleDn = f3.Get("nSearchBin_ratio_scale_down").Clone("scale_down")
    hPDFUp = f3.Get("nSearchBin_ratio_pdf_up").Clone("pdf_up")
    hPDFDn = f3.Get("nSearchBin_ratio_pdf_down").Clone("pdf_down")

    fout.cd()
    hScaleUp.Write()
    hScaleDn.Write()
    hPDFUp.Write()
    hPDFDn.Write()

    # Get central, MC stats, closure, trigger, MEU, and JEU
    f4 = TFile("/uscms/home/nstrobbe/nobackup/HadronicStop/DataTest/CMSSW_7_4_8/src/ZInvisible/Tools/condor/dataplots_muon_Jan24.root")
    hMC = f4.Get("nSearchBin/NJetWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nusingle")
    hMCstats = hMC.Clone("MC_stats")
    for i in xrange(1, 46):
        if hMC.GetBinContent(i) > 0.00000001:
            hMCstats.SetBinContent(i, hMC.GetBinError(i) / hMC.GetBinContent(i))
        else:
            hMCstats.SetBinContent(i, 0.0)
    fout.cd()
    hMCstats.Write()

    hClosureZ = f4.Get("nSearchBin/nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nusingle")
    hClosureDY = f4.Get("nSearchBin/nSearchBinnSearchBinnSearchBinDY#rightarrow#mu#mu no #mu, Z eff+accsingle")
    hClosureRatio = hClosureZ.Clone("MC_closure")
    hClosureRatio.Add(hClosureDY, -1)
    hClosureRatio.Divide(hClosureZ)
    for i in xrange(1, 46):
        hClosureRatio.SetBinContent(i, abs(hClosureRatio.GetBinContent(i)))
    fout.cd()
    hClosureRatio.Write()

    hJECNom = f4.Get("nSearchBin/syst_JESUncert_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Njet+norm weightsingle")
    hJECUp = f4.Get("nSearchBinJEUUp/syst_JESUncert_nSearchBinnSearchBinJEUUpnSearchBinJEUUpZ#rightarrow#nu#nu JEC Upsingle")
    hJECDn = f4.Get("nSearchBinJEUDn/syst_JESUncert_nSearchBinnSearchBinJEUDnnSearchBinJEUDnZ#rightarrow#nu#nu JEC Downsingle")

    hMECNom = f4.Get("nSearchBin/syst_MESUncert_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Njet+norm weightsingle")
    hMECUp = f4.Get("nSearchBinMEUUp/syst_MESUncert_nSearchBinnSearchBinMEUUpnSearchBinMEUUpZ#rightarrow#nu#nu MEC Upsingle")
    hMECDn = f4.Get("nSearchBinMEUDn/syst_MESUncert_nSearchBinnSearchBinMEUDnnSearchBinMEUDnZ#rightarrow#nu#nu MEC Downsingle")

    hJECUp_ratio = hJECUp.Clone("jeu_up")
    hJECUp_ratio.Add(hJECNom, -1.0)
    hJECUp_ratio.Divide(hJECNom)

    hJECDn_ratio = hJECDn.Clone("jeu_down")
    hJECDn_ratio.Add(hJECNom, -1.0)
    hJECDn_ratio.Divide(hJECNom)

    hMECUp_ratio = hMECUp.Clone("meu_up")
    hMECUp_ratio.Add(hMECNom, -1.0)
    hMECUp_ratio.Divide(hMECNom)

    hMECDn_ratio = hMECDn.Clone("meu_down")
    hMECDn_ratio.Add(hMECNom, -1.0)
    hMECDn_ratio.Divide(hMECNom)

    fout.cd()
    hJECUp_ratio.Write()
    hJECDn_ratio.Write()
    hMECUp_ratio.Write()
    hMECDn_ratio.Write()

    hTrigNom = f4.Get("nSearchBin/TriggerWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Trigger weight Centralsingle")
    hTrigUp =  f4.Get("nSearchBin/TriggerWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Trigger weight Upsingle")
    hTrigDn =  f4.Get("nSearchBin/TriggerWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Trigger weight Downsingle")

    hTrigUp_ratio = hTrigUp.Clone("trig_up")
    hTrigUp_ratio.Add(hTrigNom, -1.0)
    hTrigUp_ratio.Divide(hTrigNom)

    hTrigDn_ratio = hTrigDn.Clone("trig_dn")
    hTrigDn_ratio.Add(hTrigNom, -1.0)
    hTrigDn_ratio.Divide(hTrigNom)
    
    fout.cd()
    hTrigUp_ratio.Write()
    hTrigDn_ratio.Write()

    hPrediction = f4.Get("nSearchBin/TriggerWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Njet+norm weightsingle").Clone("central_prediction")
    fout.cd()
    hPrediction.Write()

    # make proto data card 

    hJEC_ratio_sym = hJECUp_ratio.Clone("hJEC_ratio_sym")
    hMEC_ratio_sym = hMECUp_ratio.Clone("hMEC_ratio_sym")
    hScale_sym     = hScaleUp.Clone("hScale_sym")
    hPDF_sym       = hPDFUp.Clone("hPDF_sym")
    hTrig_sym      = hTrigUp_ratio.Clone("hTrig_sym")

    for i in xrange(1, 46):
        hJEC_ratio_sym.SetBinContent(i, max(abs(hJECUp_ratio.GetBinContent(i)),  abs(hJECDn_ratio.GetBinContent(i))))
        hMEC_ratio_sym.SetBinContent(i, max(abs(hMECUp_ratio.GetBinContent(i)),  abs(hMECDn_ratio.GetBinContent(i))))
        hScale_sym    .SetBinContent(i, max(abs(hScaleUp.GetBinContent(i)),      abs(hScaleDn.GetBinContent(i))))
        hPDF_sym      .SetBinContent(i, max(abs(hPDFUp.GetBinContent(i)),        abs(hPDFDn.GetBinContent(i))))
        hTrig_sym     .SetBinContent(i, max(abs(hTrigUp_ratio.GetBinContent(i)), abs(hTrigDn_ratio.GetBinContent(i))))

    fout.cd()
    hJEC_ratio_sym.Write()
    hMEC_ratio_sym.Write()
    hScale_sym    .Write()
    hPDF_sym      .Write()
    hTrig_sym     .Write()

    hOther = hJEC_ratio_sym.Clone("hOther")
    for i in xrange(1, 46):
        hOther.SetBinContent(i, sqrt(hJEC_ratio_sym.GetBinContent(i)**2 + hMEC_ratio_sym.GetBinContent(i)**2 + hScale_sym.GetBinContent(i)**2 + hPDF_sym.GetBinContent(i)**2 + hTrig_sym.GetBinContent(i)**2))
    hOther.Write()

    hists = [("syst_unc_shape_central_up",   hShape_final), 
             ("syst_unc_shape_central_dn",   hShape_final), 
             ("syst_unc_shape_stat_up",      hShapeStat), 
             ("syst_unc_shape_stat_dn",      hShapeStat), 
             #("stat_unc_up",                 hMCstats), 
             #("stat_unc_dn",                 hMCstats), 
             ("syst_unc_jeu_up",             hJEC_ratio_sym),
             ("syst_unc_jeu_dn",             hJEC_ratio_sym), 
             ("syst_unc_meu_up",             hMEC_ratio_sym), 
             ("syst_unc_meu_dn",             hMEC_ratio_sym), 
             ("syst_unc_scale_up",           hScale_sym),
             ("syst_unc_scale_dn",           hScale_sym), 
             ("syst_unc_pdf_up",             hPDF_sym),
             ("syst_unc_pdf_dn",             hPDF_sym),
             ("syst_unc_trig_up",            hTrig_sym),
             ("syst_unc_trig_dn",            hTrig_sym),
             ]
    
    print "luminosity = 2153.74"
    print "channels = 45"
    print "sample = zinv"
    print ""

    print "%-25s = %s"%("channel", ' '.join(["%8s" % ("bin%i" % i) for i in xrange(1, 46)]))
    print ""

    print "%-25s = %s"%("rate", ' '.join(["%8.5f" % hPrediction.GetBinContent(i) for i in xrange(1, 46)]))
    print ""

    print "%-25s = %s"%("cs_event", ' '.join(["%8.0f" % math.floor(hNEff.GetBinContent(i)) for i in xrange(1, 46)]))

    data = []
    for i in xrange(1, 46):
        if hNEff.GetBinContent(i) > 0:
            data.append("%8.5f" % (hPrediction.GetBinContent(i)/math.floor(hNEff.GetBinContent(i))))
        else:
            data.append("%8.5f" % 0.00)
    print "%-25s = %s"%("avg_weight", ' '.join(data))

    print ""
    print "stat_unc_up = xxx yy zz"
    print "stat_unc_dn = xxx yy zz"
    print ""

    for (name, h) in hists:
        print "%-25s = %s"%(name, ' '.join(["%8.5f" % (h.GetBinContent(i)) for i in xrange(1, 46)]))

    fout.Close()



def systScalePDF(filename):
    # Open file
    f = TFile.Open(filename)

    # Get scale and PDF histograms
    hScale     = f.Get("nSearchBin/syst_ScaleWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Njet+norm weightsingle")
    hScaleUp   = f.Get("nSearchBin/syst_ScaleWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Scale weight Upsingle")
    hScaleDown = f.Get("nSearchBin/syst_ScaleWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Scale weight Downsingle")
    hPDF     = f.Get("nSearchBin/syst_PDFWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu Njet+norm weightsingle")
    hPDFUp   = f.Get("nSearchBin/syst_PDFWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu PDF weight Upsingle")
    hPDFDown = f.Get("nSearchBin/syst_PDFWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nu PDF weight Downsingle")

    # Get normalization histograms
    hnScale     = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwscale_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu Njet+norm weightsingle")
    hnScaleUp   = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwscale_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu Scale weight Upsingle")
    hnScaleDown = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwscale_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu Scale weight Downsingle")
    hnPDF     = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwpdf_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu Njet+norm weightsingle")
    hnPDFUp   = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwpdf_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu PDF weight Upsingle")
    hnPDFDown = f.Get("cntNJetsPt30Eta24Zinv/DataMCwwpdf_SingleMuon_nj_0b_blnotagcntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvZ#rightarrow#nu#nu PDF weight Downsingle")

    # Determine scale factors
    sf_Scale     = hnScale.Integral()
    sf_ScaleUp   = hnScaleUp.Integral()
    sf_ScaleDown = hnScaleDown.Integral()
    sf_PDF     = hnPDF.Integral()
    sf_PDFUp   = hnPDFUp.Integral()
    sf_PDFDown = hnPDFDown.Integral()

    # Normalize histograms
    hScale.Scale(1./sf_Scale)
    hScaleUp.Scale(1./sf_ScaleUp)
    hScaleDown.Scale(1./sf_ScaleDown)
    hPDF.Scale(1./sf_PDF)
    hPDFUp.Scale(1./sf_PDFUp)
    hPDFDown.Scale(1./sf_PDFDown)

    # Get the ratios (Up - nominal)/nominal etc
    hratio_ScaleUp = hScaleUp.Clone("nSearchBin_ratio_scale_up")
    hratio_ScaleUp.Add(hScale,-1.0)
    hratio_ScaleUp.Divide(hScale)
    hratio_ScaleDown = hScale.Clone("nSearchBin_ratio_scale_down")
    hratio_ScaleDown.Add(hScaleDown,-1.0)
    hratio_ScaleDown.Divide(hScale)
    hratio_PDFUp = hPDFUp.Clone("nSearchBin_ratio_pdf_up")
    hratio_PDFUp.Add(hPDF,-1.0)
    hratio_PDFUp.Divide(hPDF)
    hratio_PDFDown = hPDF.Clone("nSearchBin_ratio_pdf_down")
    hratio_PDFDown.Add(hPDFDown,-1.0)
    hratio_PDFDown.Divide(hPDF)

    # Write the histograms
    fout = TFile.Open("syst_scalePDF.root","RECREATE")
    fout.cd()

    hScale.Write("nSearchBin_scale_nominal")
    hScaleUp.Write("nSearchBin_scale_up")
    hScaleDown.Write("nSearchBin_scale_down")
    hratio_ScaleUp.Write()
    hratio_ScaleDown.Write()
    hPDF.Write("nSearchBin_pdf_nominal")
    hPDFUp.Write("nSearchBin_pdf_up")
    hPDFDown.Write("nSearchBin_pdf_down")
    hratio_PDFUp.Write()
    hratio_PDFDown.Write()

    # Print some info
    print "Scale uncertainty: "
    for i in range(hratio_ScaleUp.GetNbinsX()):
        print "Bin %s: %.3f, %.3f" % (i, hratio_ScaleUp.GetBinContent(i+1), hratio_ScaleDown.GetBinContent(i+1))
    print "PDF uncertainty: "
    for i in range(hratio_PDFUp.GetNbinsX()):
        print "Bin %s: %.3f, %.3f" % (i, hratio_PDFUp.GetBinContent(i+1), hratio_PDFDown.GetBinContent(i+1))

    fout.Close()
    f.Close()

def extrapolationSyst():
    print "hello"


if __name__ ==  "__main__":

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Grab histogram from FILE", metavar="FILE")
    parser.add_option("--norm", dest="normweight", default=False,
                      help="Compute the normalization weight", action='store_true')
    parser.add_option("--njet", dest="njetweight", default=False,
                      help="Compute the njet weights", action='store_true')
    parser.add_option("--shape", dest="shape", default=False,
                      help="Compute the shape systematics", action='store_true')
    parser.add_option("--systHarvest", dest="systHarvest", default=False,
                      help="Harvest systematics", action='store_true')
    parser.add_option("--systScale", dest="systScale", default=False,
                      help="Grab information for scale and PDF systematics", action='store_true')

    (options, args) = parser.parse_args()

    if options.njetweight:
        njetWeights(options.filename)
    if options.normweight:
        normWeight(options.filename)
    if options.shape:
        shapeSyst(options.filename)
    if options.systHarvest:
        systHarvest()
    if options.systScale:
        systScalePDF(options.filename)        
