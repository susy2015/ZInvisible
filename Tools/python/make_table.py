# make_table.py
import ROOT
import os

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)



# Reference: 
# https://github.com/mkilpatr/EstToolsSUSY/blob/SBv4/SUSYNano19/getUncertainty.py

pred_total_name = 'hpred'
all_samples=('ttbarplusw', 'znunu', 'ttZ', 'diboson', 'qcd')
graph_names=('httbar_stack_5', 'hznunu_stack_4', 'httz_stack_2', 'hdiboson_stack_1', 'hqcd_stack_3')
table_header='Search region & \\met [GeV]  &  Lost lepton  &  \\znunu  & rare & QCD  &  total SM  &  $N_{\\rm data}$  \\\\ \n'

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
    'lowptisr': r'($300\leq\ptisr<500$\,GeV)',
    'ntgeq1': r'($\nt\geq1$)',
    'nt2': r'($\nt=2$)',
    'nivf0': r'($\nsv=0$)',
    'nivf1': r'($\nsv\geq1$)',
    'nw2': r'($\nw=2$)',
    'nj2to5': r'($2\leq\nj\leq5$)',
    'nb2': r'($\nb\geq2$)',
    'nbeq2': r'($\nb=2$)',
    'nb3': r'($\nb\geq3$)',
    'nb1': r'($\nb=1$)',
    'nbgeq1': r'($\nb\geq1$)',
    'nb0': r'($\nb=0$)',
    'nrt2': r'($\nrt=2$)',
    'medptisr': r'($\ptisr\geq300$\,GeV)',
    'highptisr': r'($\ptisr\geq500$\,GeV)',
    'nj7': r'($\nj\geq7$)',
    'highptb': r'($\ptb\geq70$\,GeV)',
    'hm': r'(high \dm)',
    'nw0': r'($\nw=0$)',
    'nwgeq1': r'($\nw\geq1$)',
    'nw1': r'($\nw=1$)',
    'nrt0': r'($\nrt=0$)',
    'nrt1': r'($\nrt=1$)',
    'lowptb': r'($\ptb<40$\,GeV)',
    'medptb': r'($40<\ptb<70$\,GeV)',
    'nt0': r'($\nt=0$)',
    'lm': r'(low \dm)',
    'lowptb12': r'($\ptbonetwo<80$\,GeV)',
    'highptb12': r'($\ptbonetwo\geq140$\,GeV)',
    'lowmtb': r'($\mtb<175$~\GeV)',
    'highmtb': r'($\mtb\geq175$~\GeV)',
    'nt1': r'($\nt=1$)',
    'medptb12': r'($80<\ptbonetwo<140$\,GeV)',
    'nrtgeq1': r'($\nrt\geq1$)',
    'nj6': r'($\nj\geq6$)',
    'nrtntnwgeq2': r'($(\nt+\nrt+\nw)\geq2$)',
    'nrtntnwgeq3': r'($(\nt+\nrt+\nw)\geq3$)',
    'htlt1000': r'($\Ht<1000$)',
    'htgt1000': r'($\Ht\geq1000$)',
    'ht1000to1500': r'($1000\leq\Ht<1500$)',
    'htgt1500': r'($\Ht\geq1500$)',
    'htlt1300': r'($\Ht<1300$)',
    'htgt1300': r'($\Ht\geq1300$)',
    'lmNoDPhi': r'(low $\Delta m$)',
    'hmNoDPhi': r'(high $\Delta m$)',
    'MET': r'',
    }


class Table:
    def __init__(self):
        self.yields_data = {}
        self.yields      = {}
        self.allVals     = {}
        
    # read in yields
    def readYields(self, pred_file):
        ''' Read in predicted bkg yields and the stat. unc. 
        
        pred_file -- input root file 
        '''
        if pred_file:
            f = ROOT.TFile(pred_file)
            for hname, sample in zip(graph_names, all_samples):
                h = f.Get(hname)
                for ibin in xrange(0, h.GetNbinsX()):
                    bin = binlist[ibin]
                    if bin not in self.yields:
                        self.yields[bin] = {}
                        statUnc_pieces[bin] = {}
                    y = h.GetBinContent(ibin)
                    e_up = h.GetBinError(ibin)
                    e_low = h.GetBinError(ibin)
                    self.yields[bin][sample] = y
                    if sample == 'rare': statUnc_pieces[bin][sample] = (min(e_up,y), min(e_up,y))  # don't want MC stat unc > 100%
                    else :               statUnc_pieces[bin][sample] = (e_low, e_up)
            #h = f.Get('data')
            #for ibin in xrange(0, h.GetNbinsX()):
                    bin = binlist[ibin]
                    #self.yields_data[bin] = (h.GetBinContent(ibin+1), h.GetBinError(ibin+1))
                    self.yields_data[bin] = (1, 1)
            # get total pred (w/ asymmetric errors)
            h = f.Get(pred_total_name)
            for ibin in xrange(0, h.GetNbinsX()):
                bin = binlist[ibin]
                e_up = h.GetBinError(ibin)
                e_low = h.GetBinError(ibin)
                statUnc[bin] = (e_low, e_up)
            f.Close()
        else:
            for ibin in xrange(0, len(binlist)):
                bin = binlist[ibin]
                self.yields[bin] = {}
                for hname, sample in zip(graph_names, all_samples):
                    #self.yields[bin][sample] = y
                    #self.yields_data[bin] = (h.GetBinContent(ibin+1), h.GetBinError(ibin+1))
                    self.yields[bin][sample] = 1
                    self.yields_data[bin] = (1, 1)



    def writeFullUnc(self, pred_file):
        ''' Update the input root file, add a hist with total prediction and full uncertainty. '''
        if pred_file:
            f = ROOT.TFile(pred_file, 'UPDATE')
            h = TGraphAsymmErrors(f.Get(pred_total_name).Clone('bkgtotal_unc_sr'))
            h_pieces = {}
            for hname, sample in zip(graph_names, all_samples):
                h_pieces[sample] = TGraphAsymmErrors(f.Get(hname).Clone(sample+'_unc_sr'))
            print "%30s %10s %16s" % ('bin', 'total pred', 'total unc.')
            for ibin in xrange(0, h.GetN()):
                bin = binlist[ibin]
                val = h.GetY()[ibin]
                e_low, e_up = fullUnc[bin]
                h.SetPointEYlow(ibin, e_low)
                h.SetPointEYhigh(ibin, e_up)
                print "%30s %10.2f +%8.2f -%8.2f" % (bin, val, e_up, e_low)
                self.allVals[bin] = {'bkg':(val,e_low,e_up)}
                for sample in all_samples:
                    val = self.yields[bin][sample]
                    e_low, e_up = fullUnc_pieces[sample][bin]
                    h_pieces[sample].SetPointEYlow(ibin, e_low)
                    h_pieces[sample].SetPointEYhigh(ibin, e_up)
                    self.allVals[bin][sample] = (val,e_low,e_up)  
            h.Write('bkgtotal_unc_sr', ROOT.TObject.kOverwrite)
            for sample in all_samples : h_pieces[sample].Write(sample+'_unc_sr', ROOT.TObject.kOverwrite)
            f.Close()
        else:
            for ibin in xrange(0, len(binlist)):
                bin = binlist[ibin]
                self.allVals[bin] = {}
                #self.allVals[bin]['bkg'] = (val,e_low,e_up)
                self.allVals[bin]['bkg'] = (1,1,1)
                for sample in all_samples:
                    #self.allVals[bin][sample] = (val,e_low,e_up)  
                    self.allVals[bin][sample] = (1,1,1)  



    def makeYieldTable(self, output="pred_sr.tex"):
        ''' Make a Latex-formatted table with each bkg plus unc, total bkg plus unc, and observed data for every bin. '''
        s  = self.beginDocument()
        s += self.beginTable()
        s += table_header
        s += '\\hline\n'
        s += self.makeTable()
        s += self.endTable()
        s += self.endDocument()
        print '\nprinting yield table...\n'
        print s
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
        s += '\\tabcolsep=0.01cm\n'
        s += '\\centering\n'
        return s

    def endDocument(self):
        '''end latex document'''
        s = '\\end{document}'
        return s
       
    def beginTable(self):
        '''Add a break between the bins to fit on each page'''
        s  = '\\begin{table}[!h]\n'
        s += '\\begin{center}\n'
        #s += '\\resizebox*{0.6\\textwidth}{!}{\n'
        s += '\\begin{tabular}{|c||c||c|c|c|c|c|c|}\n'
        s += '\\hline\n'
        return s
    
    def endTable(self):
        '''Add a break between the bins to fit on each page'''
        s  = '\\hline\n'
        s += '\\end{tabular}\n'
        s += '\\caption[]{}\n'
        s += '\\end{center}\n'
        s += '\\end{table}\n'
        return s
    
    
    def makeTable(self):
        ''' Put together the table chunk for the given nj,nb,mtb,nt mega-bin. '''
        sections=[]
        s=''
        ibin=0
        for bin in binlist: 
            sec, met = bin.lstrip('bin_').rsplit('_', 1)
            met = met.lstrip("pt")
            #print "sec = {0}, met = {1}".format(sec, met)
            if sec not in sections:
                sections.append(sec)
                s += self.chunkHeader(sec)
    #         metbins = binMap[sec]['bin']
    #         print metbins
    #         idx = metbins.index(int(met))
            xlow, xhigh = met.lstrip('met').split('to')
            metlabel = r'$>%s$'%xlow if xhigh=='inf' else '$-$'.join([xlow, xhigh])
            s += '%d & '%ibin
            ibin = ibin+1
            s += metlabel
            for bkg in list(all_samples)+['bkg']:
                if bkg == 'diboson': continue
                n, e_low, e_up = self.allVals[bin][bkg]
                if bkg == 'ttZ':
                    n1, e1_low, e1_up = self.allVals[bin]["diboson"]
                    n += n1
                    #e_low = sumUnc([e_low, e1_low])
                    #e_up  = sumUnc([e_up, e1_up])
                    e_low = 1 
                    e_up  = 1
                s += self.formatPrediction(n,e_low,e_up)
            n, e = self.yields_data[bin]
            s += ' & ' + str(int(n))
            s += ' \\\\ \n'
            if ibin == 53 or ibin == 94 or ibin == 135:
                s += self.endTable()
                s += self.beginTable()
                s += table_header
                s += self.endTable()
                s += self.beginTable()
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
        cats = sec.split('_')
        labs = [labelMap[c] for c in cats]
        ncolumn = len(all_samples)+3
        s  = '\\hline\n'
        s += '\\multicolumn{'+str(ncolumn)+'}{c}{'
        s += ', '.join(labs)
        s += '} \\\\ \n'
        s += '\\hline\n' 
        return s




def main():
    T = Table()
    T.readYields("")
    T.writeFullUnc("")
    T.makeYieldTable("latex_files/zinv_pred_sr.tex")

if __name__ == "__main__":
    main()








