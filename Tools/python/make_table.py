# make_table.py
import ROOT
import os

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Reference: 
# https://github.com/mkilpatr/EstToolsSUSY/blob/SBv4/SUSYNano19/getUncertainty.py

all_values=("norm_tex", "shape_tex", "mc_tex", "pred_tex")
table_header='Search region & \\met [GeV]  &  $R_Z$  &  $S_{\gamma}$  & $N_{\\rm MC}$ & $N_{\\rm pred}$ \\\\ \n'

# ordered bin list
binlist = (
    'bin_lm_nb0_nivf0_highptisr_nj2to5_MET_pt450to550', 
    'bin_lm_nb0_nivf0_highptisr_nj2to5_MET_pt550to650', 
    'bin_lm_nb0_nivf0_highptisr_nj2to5_MET_pt650to750', 
    'bin_lm_nb0_nivf0_highptisr_nj2to5_MET_pt750toinf', 
    'bin_lm_nb0_nivf0_highptisr_nj6_MET_pt450to550', 
    'bin_lm_nb0_nivf0_highptisr_nj6_MET_pt550to650', 
    'bin_lm_nb0_nivf0_highptisr_nj6_MET_pt650to750', 
    'bin_lm_nb0_nivf0_highptisr_nj6_MET_pt750toinf', 
    'bin_lm_nb0_nivf1_highptisr_nj2to5_MET_pt450to550', 
    'bin_lm_nb0_nivf1_highptisr_nj2to5_MET_pt550to650', 
    'bin_lm_nb0_nivf1_highptisr_nj2to5_MET_pt650to750', 
    'bin_lm_nb0_nivf1_highptisr_nj2to5_MET_pt750toinf', 
    'bin_lm_nb0_nivf1_highptisr_nj6_MET_pt450to550', 
    'bin_lm_nb0_nivf1_highptisr_nj6_MET_pt550to650', 
    'bin_lm_nb0_nivf1_highptisr_nj6_MET_pt650to750', 
    'bin_lm_nb0_nivf1_highptisr_nj6_MET_pt750toinf', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_lowptb_MET_pt300to400', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_lowptb_MET_pt400to500', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_lowptb_MET_pt500to600', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_lowptb_MET_pt600toinf', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_medptb_MET_pt300to400', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_medptb_MET_pt400to500', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_medptb_MET_pt500to600', 
    'bin_lm_nb1_nivf0_lowmtb_lowptisr_medptb_MET_pt600toinf', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_lowptb_MET_pt450to550', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_lowptb_MET_pt550to650', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_lowptb_MET_pt650to750', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_lowptb_MET_pt750toinf', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_medptb_MET_pt450to550', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_medptb_MET_pt550to650', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_medptb_MET_pt650to750', 
    'bin_lm_nb1_nivf0_lowmtb_highptisr_medptb_MET_pt750toinf', 
    'bin_lm_nb1_nivf1_lowmtb_medptisr_lowptb_MET_pt300to400', 
    'bin_lm_nb1_nivf1_lowmtb_medptisr_lowptb_MET_pt400to500', 
    'bin_lm_nb1_nivf1_lowmtb_medptisr_lowptb_MET_pt500toinf', 
    'bin_lm_nb2_lowmtb_lowptisr_lowptb12_MET_pt300to400', 
    'bin_lm_nb2_lowmtb_lowptisr_lowptb12_MET_pt400to500', 
    'bin_lm_nb2_lowmtb_lowptisr_lowptb12_MET_pt500toinf', 
    'bin_lm_nb2_lowmtb_lowptisr_medptb12_MET_pt300to400', 
    'bin_lm_nb2_lowmtb_lowptisr_medptb12_MET_pt400to500', 
    'bin_lm_nb2_lowmtb_lowptisr_medptb12_MET_pt500toinf', 
    'bin_lm_nb2_lowmtb_lowptisr_highptb12_nj7_MET_pt300to400', 
    'bin_lm_nb2_lowmtb_lowptisr_highptb12_nj7_MET_pt400to500', 
    'bin_lm_nb2_lowmtb_lowptisr_highptb12_nj7_MET_pt500toinf', 
    'bin_lm_nb2_lowmtb_highptisr_lowptb12_MET_pt450to550', 
    'bin_lm_nb2_lowmtb_highptisr_lowptb12_MET_pt550to650', 
    'bin_lm_nb2_lowmtb_highptisr_lowptb12_MET_pt650toinf', 
    'bin_lm_nb2_lowmtb_highptisr_medptb12_MET_pt450to550', 
    'bin_lm_nb2_lowmtb_highptisr_medptb12_MET_pt550to650', 
    'bin_lm_nb2_lowmtb_highptisr_medptb12_MET_pt650toinf', 
    'bin_lm_nb2_lowmtb_highptisr_highptb12_nj7_MET_pt450to550', 
    'bin_lm_nb2_lowmtb_highptisr_highptb12_nj7_MET_pt550to650', 
    'bin_lm_nb2_lowmtb_highptisr_highptb12_nj7_MET_pt650toinf', 
    'bin_hm_nb1_lowmtb_nj7_nrtgeq1_MET_pt250to300', 
    'bin_hm_nb1_lowmtb_nj7_nrtgeq1_MET_pt300to400', 
    'bin_hm_nb1_lowmtb_nj7_nrtgeq1_MET_pt400to500', 
    'bin_hm_nb1_lowmtb_nj7_nrtgeq1_MET_pt500toinf', 
    'bin_hm_nb2_lowmtb_nj7_nrtgeq1_MET_pt250to300', 
    'bin_hm_nb2_lowmtb_nj7_nrtgeq1_MET_pt300to400', 
    'bin_hm_nb2_lowmtb_nj7_nrtgeq1_MET_pt400to500', 
    'bin_hm_nb2_lowmtb_nj7_nrtgeq1_MET_pt500toinf', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt450to550', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt550toinf', 
    'bin_hm_nb2_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt250to350', 
    'bin_hm_nb2_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt350to450', 
    'bin_hm_nb2_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt450to550', 
    'bin_hm_nb2_highmtb_nt0_nrt0_nw0_htgt1000_MET_pt550toinf', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_MET_pt250to550', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_MET_pt550to650', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htlt1000_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_MET_pt250to550', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_MET_pt550to650', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_ht1000to1500_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_MET_pt250to550', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_MET_pt550to650', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nw0_htgt1500_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htlt1300_MET_pt450toinf', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrt0_nwgeq1_htgt1300_MET_pt450toinf', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_MET_pt450to550', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_MET_pt550to650', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htlt1000_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_MET_pt450to550', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_MET_pt550to650', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_ht1000to1500_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_MET_pt250to350', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_MET_pt350to450', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_MET_pt450to550', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_MET_pt550to650', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nw0_htgt1500_MET_pt650toinf', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1_MET_pt250to550', 
    'bin_hm_nb1_highmtb_ntgeq1_nrt0_nwgeq1_MET_pt550toinf', 
    'bin_hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0_MET_pt250to550', 
    'bin_hm_nb1_highmtb_ntgeq1_nrtgeq1_nw0_MET_pt550toinf', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1_MET_pt250to550', 
    'bin_hm_nb1_highmtb_nt0_nrtgeq1_nwgeq1_MET_pt550toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt250to550', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt250to550', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt250to550', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htlt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw1_htgt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt450to550', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt450to550', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt450to550', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt550to650', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt650toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw1_MET_pt250to550', 
    'bin_hm_nbeq2_highmtb_nt1_nrt0_nw1_MET_pt550toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htlt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_MET_pt250to350', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_MET_pt350to450', 
    'bin_hm_nbeq2_highmtb_nt1_nrt1_nw0_htgt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw1_MET_pt250to550', 
    'bin_hm_nbeq2_highmtb_nt0_nrt1_nw1_MET_pt550toinf', 
    'bin_hm_nbeq2_highmtb_nt2_nrt0_nw0_MET_pt250to450', 
    'bin_hm_nbeq2_highmtb_nt2_nrt0_nw0_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt0_nw2_MET_pt250toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300_MET_pt250to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt2_nw0_htlt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300_MET_pt250to450', 
    'bin_hm_nbeq2_highmtb_nt0_nrt2_nw0_htgt1300_MET_pt450toinf', 
    'bin_hm_nbeq2_highmtb_nrtntnwgeq3_MET_pt250toinf', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htlt1000_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_ht1000to1500_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw0_htgt1500_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt0_nw1_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt0_nrt0_nw1_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt0_nrt0_nw1_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htlt1000_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_ht1000to1500_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt350to550', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw0_htgt1500_MET_pt550toinf', 
    'bin_hm_nb3_highmtb_nt1_nrt0_nw1_MET_pt250toinf', 
    'bin_hm_nb3_highmtb_nt1_nrt1_nw0_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt1_nrt1_nw0_MET_pt350toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt1_nw1_MET_pt250toinf', 
    'bin_hm_nb3_highmtb_nt2_nrt0_nw0_MET_pt250toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt0_nw2_MET_pt250toinf', 
    'bin_hm_nb3_highmtb_nt0_nrt2_nw0_MET_pt250to350', 
    'bin_hm_nb3_highmtb_nt0_nrt2_nw0_MET_pt350toinf', 
    'bin_hm_nb3_highmtb_nrtntnwgeq3_MET_pt250toinf'
)

# lookup map for latex labels
labelMap = {
    'lowptisr': r'$300\leq\ptisr<500$\,GeV',
    'ntgeq1': r'$\nt\geq1$',
    'nt2': r'$\nt=2$',
    'nivf0': r'$\nsv=0$',
    'nivf1': r'$\nsv\geq1$',
    'nw2': r'$\nw=2$',
    'nj2to5': r'$2\leq\nj\leq5$',
    'nb2': r'$\nb\geq2$',
    'nbeq2': r'$\nb=2$',
    'nb3': r'$\nb\geq3$',
    'nb1': r'$\nb=1$',
    'nbgeq1': r'$\nb\geq1$',
    'nb0': r'$\nb=0$',
    'nrt2': r'$\nrt=2$',
    'medptisr': r'$\ptisr\geq300$\,GeV',
    'highptisr': r'$\ptisr\geq500$\,GeV',
    'nj7': r'$\nj\geq7$',
    'highptb': r'$\ptb\geq70$\,GeV',
    'hm': r'high \dm',
    'nw0': r'$\nw=0$',
    'nwgeq1': r'$\nw\geq1$',
    'nw1': r'$\nw=1$',
    'nrt0': r'$\nrt=0$',
    'nrt1': r'$\nrt=1$',
    'lowptb': r'$\ptb<40$\,GeV',
    'medptb': r'$40<\ptb<70$\,GeV',
    'nt0': r'$\nt=0$',
    'lm': r'low \dm',
    'lowptb12': r'$\ptbonetwo<80$\,GeV',
    'highptb12': r'$\ptbonetwo\geq140$\,GeV',
    'lowmtb': r'$\mtb<175$~\GeV',
    'highmtb': r'$\mtb\geq175$~\GeV',
    'nt1': r'$\nt=1$',
    'medptb12': r'$80<\ptbonetwo<140$\,GeV',
    'nrtgeq1': r'$\nrt\geq1$',
    'nj6': r'$\nj\geq6$',
    'nrtntnwgeq2': r'$\nt+\nrt+\nw\geq2$',
    'nrtntnwgeq3': r'$\nt+\nrt+\nw\geq3$',
    'htlt1000': r'$\Ht<1000$',
    'htgt1000': r'$\Ht\geq1000$',
    'ht1000to1500': r'$1000\leq\Ht<1500$',
    'htgt1500': r'$\Ht\geq1500$',
    'htlt1300': r'$\Ht<1300$',
    'htgt1300': r'$\Ht\geq1300$',
    'lmNoDPhi': r'low $\Delta m$',
    'hmNoDPhi': r'high $\Delta m$',
    'MET': r'',
}

class Table:
    def __init__(self):
        pass
        
    def makeYieldTable(self, BinObject, total_era, output="pred_sr.tex", makeDoc=False, size=0.6):
        ''' Make a Latex-formatted table with each bkg plus unc, total bkg plus unc, and observed data for every bin. '''
        self.size = size
        s  = ""
        if makeDoc:
            s += self.beginDocument()
        # TODO: remove
        #s += self.beginTable()
        #s += table_header
        s += self.makeTable(BinObject, total_era)
        if makeDoc:
            s += self.endDocument()
        with open(output, 'w') as f:
            print >> f, s

    def beginDocument(self):
        '''begin latex document'''
        s  = '\\documentclass{article}\n'
        s += '\\usepackage[utf8]{inputenc}\n'
        s += '\\usepackage{xspace}\n'
        s += '\\usepackage{graphicx}\n'
        s += '\\usepackage{geometry}\n'
        s += '\\usepackage{longtable}\n'
        s += '\\usepackage{cancel}\n'
        s += '\\usepackage{amsmath}\n'
        s += '\\input{VariableNames.tex}\n'
        s += '\\geometry{margin=0.1cm}\n'
        s += '\\begin{document}\n'
        s += '\\footnotesize\n'
        s += '\\tabcolsep=0.2cm\n'
        s += '\\centering\n'
        return s

    def endDocument(self):
        '''end latex document'''
        s = '\\end{document}'
        return s
       
    def beginTable(self, caption="", label=""):
        '''Add a break between the bins to fit on each page'''
        # WARNING: label must go after caption for table reference to work
        s  = '\\begin{table}[!h]\n'
        s += '\\begin{center}\n'
        s += '\\caption{%s}\n' % caption
        s += '\\label{%s}\n' % label
        s += '\\resizebox*{%.2f\\textwidth}{!}{\n' % self.size
        s += '\\begin{tabular}{|c||c||c|c|c|c|}\n'
        s += '\\hline\n'
        return s
    
    def endTable(self):
        '''Add a break between the bins to fit on each page'''
        '''Include label and caption'''
        s  = '\\hline\n'
        s += '\\end{tabular}\n'
        s += '}\n'
        s += '\\end{center}\n'
        s += '\\end{table}\n'
        return s
    
    def makeTable(self, BinObject, total_era):
        ''' Put together the table chunk for the given nj,nb,mtb,nt,nw,ht mega-bin. '''
        binRanges = {0:52, 53:93, 94:134, 135:182}
        sections=[]
        s  = ""
        ibin = 0
        # WARNING: binlist contains string bin names
        for bin in binlist: 
            # put caption before table
            if ibin in binRanges:
                # bin range
                firstBin = ibin
                lastBin  = binRanges[ibin]
                # low/high dm
                region = "low \dm"
                if lastBin >= 53:
                    region = "high \dm"
                caption  = "Prediction for the \zinv background $\\left(\Np\\right)$ in {0} search bins {1}--{2}.".format(region, firstBin, lastBin)
                caption += " The normalization factor $\\left(\Rz\\right)$, shape factor $\\left(\Sg\\right)$, and number of \znunu MC events $\\left(\Nmc\\right)$ are also shown for each search bin including their statistical uncertainties."
                caption += " The uncertainty for the prediction $\\left(\Np\\right)$ is calculated by propagating the statistical uncertainties of \Rz, \Sg, and \Nmc."
                caption += " See Eq.~\\ref{eq:zinv_pred}."
                label    = "tab:zinvPredToBin{0}".format(lastBin)
                s += self.beginTable(caption, label)
                s += table_header
            sec, met = bin.lstrip('bin_').rsplit('_', 1)
            met = met.lstrip("pt")
            if sec not in sections:
                sections.append(sec)
                s += self.chunkHeader(sec)
            xlow, xhigh = met.lstrip('met').split('to')
            metlabel = r'$>%s$'%xlow if xhigh=='inf' else '$-$'.join([xlow, xhigh])
            s += '%d & ' % ibin
            s += metlabel
            for value in list(all_values):
                s += " & {0} ".format(BinObject.binValues[total_era][str(ibin)][value])
            s += ' \\\\ \n'
            # first increment ibin
            ibin += 1
            # now these are bin numbers to end table
            if ibin == 53 or ibin == 94 or ibin == 135 or ibin == 183:
                # TODO: remove
                # # last bin for previous table
                # lastBin = ibin - 1
                # # low/high dm
                # region = "low \dm"
                # if lastBin >= 53:
                #     region = "high \dm"
                # caption  = "Prediction for the \zinv background $\\left(\Np\\right)$ in {0} search bins {1}--{2}.".format(region, firstBin, lastBin)
                # caption += " The normalization factor $\\left(\Rz\\right)$, shape factor $\\left(\Sg\\right)$, and number of \znunu MC events $\\left(\Nmc\\right)$ are also shown for each search bin including their statistical uncertainties."
                # caption += " The uncertainty for the prediction $\\left(\Np\\right)$ is calculated by propagating the statistical uncertainties of \Rz, \Sg, and \Nmc."
                # caption += " See Eq.~\\ref{eq:zinv_pred}."
                # label    = "tab:zinvPredToBin{0}".format(lastBin)
                # TODO: remove
                s += self.endTable()
                # if ibin < 183:
                #     # first bin for next table 
                #     firstBin = ibin
                #     s += self.beginTable()
                #     s += table_header
        return s
    
    # formats the prediction nEvents +/- error
    def formatPrediction(self,n,e_low,e_up):
        if n>=10:
            n = str(int(round(n,0)))
            e_low = str(int(round(e_low,0)))
            e_up  = str(int(round(e_up,0)))
        elif n>=1:
            n = str(round(n,1))
            e_low = str(round(e_low,1))
            e_up  = str(round(e_up,1))
        else:
            n = str(round(n,2))
            e_low = str(round(e_low,2))
            e_up  = str(round(e_up,2))
        if n=='0.0':
            if e_up=='0.0':
                return ' & $<$0.01'
            return ' & $<$' + e_up
        if e_low==e_up:
            return ' & $ %s\\pm%s $ ' %(n, e_up)
        else:
            return ' & $ %s\,^{+%s}_{-%s} $ ' %(n, e_up, e_low)
    
    # puts together the bin header for bins of nJets, mtb, nTop (no selection on nB)
    def chunkHeader(self,sec):
        ''' Put together the mega-bin chunk header. '''
        # felines
        cats = sec.split('_')
        # canines
        # only include non-empty labels, e.g. MET has not label
        labs = [labelMap[c] for c in cats if labelMap[c]]
        ncolumn = len(all_values)+2
        s  = '\\hline\n'
        s += '\\multicolumn{'+str(ncolumn)+'}{c}{'
        s += ', '.join(labs)
        s += '} \\\\ \n'
        s += '\\hline\n' 
        return s

def main():
    print "Running \"python make_table.py\" directly is not supported."
    print "Plesae use run_modules.py"  

if __name__ == "__main__":
    main()

