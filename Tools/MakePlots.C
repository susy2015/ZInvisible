#include "SusyAnaTools/Tools/Plotter.h"
#include "SusyAnaTools/Tools/samples.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
    using namespace std;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"plot",            no_argument, 0, 'p'},
        {"savehist",        no_argument, 0, 's'},
        {"fromFile",        no_argument, 0, 'f'},
        {"condor",          no_argument, 0, 'c'},
        {"histFile",  required_argument, 0, 'H'},
        {"dataSets",  required_argument, 0, 'D'},
        {"numFiles",  required_argument, 0, 'N'}
    };

    bool doPlots = true, doSave = true, fromTuple = true, runOnCondor = false;
    string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir;
    int nFiles = -1;

    while((opt = getopt_long(argc, argv, "psfcH:D:N:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'p':
            if(doPlots) doSave  = false;
            else        doPlots = true;
            break;

        case 's':
            if(doSave) doPlots = false;
            else       doSave  = true;
            break;

        case 'f':
            fromTuple = false;
            break;

        case 'c':
            runOnCondor = true;
            break;

        case 'H':
            histFile = optarg;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;
        }
    }

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        histFile = "histoutput_" + dataSets + ".root";
        doSave = true;
        doPlots = false;
        fromTuple = true;
        sampleloc = "";
    }

    AnaSamples::SampleSet        ss(sampleloc);
    AnaSamples::SampleCollection sc(ss);

    const double zAcc = 1.0;
//    const double zAcc = 0.5954;
//    const double zAcc = 0.855;
    const double znunu_mumu_ratio = 5.942;

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_400to600"]};
        fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
    }
    else if(dataSets.size() == 0)
    {
        fileMap["IncDY"] = sc["IncDY"];
        fileMap["DYJetsToLL"]  = sc["DYJetsToLL"];
        fileMap["ZJetsToNuNu"] = sc["ZJetsToNuNu"];
        fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["DYJetsToLL_HT_200to400"] = {ss["DYJetsToLL_HT_200to400"]};
        fileMap["DYJetsToLL_HT_100to200"] = {ss["DYJetsToLL_HT_100to200"]};
        fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
        fileMap["ZJetsToNuNu_HT_200to400"] = {ss["ZJetsToNuNu_HT_200to400"]};
        fileMap["ZJetsToNuNu_HT_100to200"] = {ss["ZJetsToNuNu_HT_100to200"]};
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = sc[dataSets];
        if(dataSets.compare("DYJetsToLL") == 0)
        {
            fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
            fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
            fileMap["DYJetsToLL_HT_200to400"] = {ss["DYJetsToLL_HT_200to400"]};
            fileMap["DYJetsToLL_HT_100to200"] = {ss["DYJetsToLL_HT_100to200"]};
        }
        else if(dataSets.compare("ZJetsToNuNu"))
        {
            fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
            fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
            fileMap["ZJetsToNuNu_HT_200to400"] = {ss["ZJetsToNuNu_HT_200to400"]};
            fileMap["ZJetsToNuNu_HT_100to200"] = {ss["ZJetsToNuNu_HT_100to200"]};
        }
    }
    else if(ss[dataSets] != ss.null())
    {
        for(const auto& fsnp : sc)
        {
            if(dataSets.find(fsnp.first) !=std::string::npos)
            {
                fileMap[fsnp.first] = {ss[dataSets]};
                fileMap[dataSets] = {ss[dataSets]};
                break;
            }
        }
    }

    Plotter::DatasetSummary dsDY_ll_inc(    "DY#rightarrow#mu#mu Inc",               fileMap["IncDY"], "pdgIdZDec=13;passMuZinvSel", "");

    Plotter::DatasetSummary dsDY_ll(        "DY#rightarrow#mu#mu",                   fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_ll_base(   "DY#rightarrow#mu#mu baseline",          fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passZinvBaseline", "");
    Plotter::DatasetSummary dsDY_ll_baseNT( "DY#rightarrow#mu#mu baselineNoTag",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passZinvBaselineNoTag", "");
    Plotter::DatasetSummary dsDY_ll_zWeight("DY#rightarrow#mu#mu no #mu, Z eff",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt");
    Plotter::DatasetSummary dsDY_ll_zAcc(   "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_gen(    "DY#rightarrow#mu#mu Gen",               fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_llclean(   "DY#rightarrow#mu#mu no #mu",            fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_ll_mu30_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 30, 20 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>30;cutMuPt2>20", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_mu30_mu30( "DY#rightarrow#mu#mu #mu p_{T} > 30, 30 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>30;cutMuPt2>30", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_mu40_mu30( "DY#rightarrow#mu#mu #mu p_{T} > 40, 30 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>40;cutMuPt2>30", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_mu35_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 35, 20 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>35;cutMuPt2>20", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_mu40_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 40, 20 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>40;cutMuPt2>20", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_mu45_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 45, 20 GeV",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;cutMuPt1>45;cutMuPt2>20", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_singleMu45("DY#rightarrow#mu#mu SingleMu45eta2p1",           fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passSingleMu45",          "zEffWgt;zAccWgt");

    Plotter::DatasetSummary dsDY_ll_scaled(          "DY#rightarrow#mu#mu",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel",                         "",   znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_scaled_baseNT( "DY#rightarrow#mu#mu baselineNoTag",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passZinvBaselineNoTag", "", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_zWeight_scaled("DY#rightarrow#mu#mu no #mu, Z eff",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt",         znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_zAcc_scaled(   "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_trig_scaled(   "DY#rightarrow#mu#mu no #mu, Z eff+acc+trig", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passDiMuIsoTrig", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_llclean_scaled(   "DY#rightarrow#mu#mu no #mu",            fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "",                znunu_mumu_ratio / zAcc);

    Plotter::DatasetSummary dsDY_ll_forAcc(     "DY#rightarrow#mu#mu no Sel forAcc",                      fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMu");
    Plotter::DatasetSummary dsDY_ll_forEff(     "DY#rightarrow#mu#mu no Sel forEff",                      fileMap["DYJetsToLL"],  "pdgIdZDec=13", "");
    Plotter::DatasetSummary dsDY_ll_forEff_num( "DY#rightarrow#mu#mu no Sel N(#mu) match wgt",     fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMatchMuInAcc");
    Plotter::DatasetSummary dsDY_ll_forEff_den( "DY#rightarrow#mu#mu no Sel N(#mu) wgt",           fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMuInAcc");

    Plotter::DatasetSummary dsDY_nunu(          "Z#rightarrow#nu#nu",                fileMap["ZJetsToNuNu"], "passStandLeptVeto", "");
    Plotter::DatasetSummary dsDY_test(          "test",                              fileMap["ZJetsToNuNu"], "passLeptVeto", "");

    Plotter::DatasetSummary dsDY_ll_0t(        "DY#rightarrow#mu#mu N(t) = 0",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=0", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_1t(        "DY#rightarrow#mu#mu N(t) = 1",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=1", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_2t(        "DY#rightarrow#mu#mu N(t) = 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=2", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_3t(        "DY#rightarrow#mu#mu N(t) > 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv>2", "zEffWgt;zAccWgt");

    Plotter::DatasetSummary dsDY_ll_0b(        "DY#rightarrow#mu#mu N(b) = 0",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=0", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_1b(        "DY#rightarrow#mu#mu N(b) = 1",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=1", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_2b(        "DY#rightarrow#mu#mu N(b) = 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=2", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_3b(        "DY#rightarrow#mu#mu N(b) > 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv>2", "zEffWgt;zAccWgt");

    Plotter::DatasetSummary dsDY_nunu_0t(          "Z#rightarrow#nu#nu N(t) = 0",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;nTopCandSortedCntZinv=0", "");
    Plotter::DatasetSummary dsDY_nunu_1t(          "Z#rightarrow#nu#nu N(t) = 1",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;nTopCandSortedCntZinv=1", "");
    Plotter::DatasetSummary dsDY_nunu_2t(          "Z#rightarrow#nu#nu N(t) = 2",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;nTopCandSortedCntZinv=2", "");
    Plotter::DatasetSummary dsDY_nunu_3t(          "Z#rightarrow#nu#nu N(t) > 2",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;nTopCandSortedCntZinv>2", "");
                                        
    Plotter::DatasetSummary dsDY_nunu_0b(          "Z#rightarrow#nu#nu N(b) = 0",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_nunu_1b(          "Z#rightarrow#nu#nu N(b) = 1",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;cntCSVSZinv=1", "");
    Plotter::DatasetSummary dsDY_nunu_2b(          "Z#rightarrow#nu#nu N(b) = 2",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;cntCSVSZinv=2", "");
    Plotter::DatasetSummary dsDY_nunu_3b(          "Z#rightarrow#nu#nu N(b) > 2",    fileMap["ZJetsToNuNu"], "passStandLeptVeto;cntCSVSZinv>2", "");

    Plotter::DatasetSummary dsDY_ll_HT_100to200(  "DY#rightarrow#mu#mu 100 < H_{T} < 200 GeV",    fileMap["DYJetsToLL_HT_100to200"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_200to400(  "DY#rightarrow#mu#mu 200 < H_{T} < 400 GeV",    fileMap["DYJetsToLL_HT_200to400"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_400to600(  "DY#rightarrow#mu#mu 400 < H_{T} < 600 GeV",    fileMap["DYJetsToLL_HT_400to600"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_600toInf(  "DY#rightarrow#mu#mu H_{T} > 600 GeV",          fileMap["DYJetsToLL_HT_600toInf"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);

    Plotter::DatasetSummary dsDY_nunu_HT_100to200("DY#rightarrow#nu#nu 100 < H_{T} < 200 GeV",    fileMap["ZJetsToNuNu_HT_100to200"],  "passStandLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_200to400("DY#rightarrow#nu#nu 200 < H_{T} < 400 GeV",    fileMap["ZJetsToNuNu_HT_200to400"],  "passStandLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_400to600("DY#rightarrow#nu#nu 400 < H_{T} < 600 GeV",    fileMap["ZJetsToNuNu_HT_400to600"],  "passStandLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_600toInf("DY#rightarrow#nu#nu H_{T} > 600 GeV",          fileMap["ZJetsToNuNu_HT_600toInf"],  "passStandLeptVeto", "");



    vector<Plotter::HistSummary> vh;

    //Basic Raw compairson plots  
    vh.push_back(PHS("inc_mht",     {PDC("single", "mht",                   {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "M(H_{t}) [GeV]",          "Events"));
    vh.push_back(PHS("inc_met",     {PDC("single", "met",                   {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "MET [GeV]",               "Events"));
    vh.push_back(PHS("inc_ht",      {PDC("single", "ht",                    {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "H_{T} [GeV]",             "Events"));
    vh.push_back(PHS("inc_genht",   {PDC("single", "genHt",                 {dsDY_ll_inc, dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "gen H_{T} [GeV]",             "Norm Events"));

    vh.push_back(PHS("minljdR",   {PDC("single", "minljdR",                 {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2.0,  true,  false,  "min #DeltaR(l,j)",             "Events"));

    vh.push_back(PHS("mht",         {PDC("single", "mht",                   {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "M(H_{t}) [GeV]",          "Norm Events"));
    vh.push_back(PHS("mt2",         {PDC("single", "mt2",                   {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MT2 [GeV]",               "Norm Events"));
    vh.push_back(PHS("met",         {PDC("single", "met",                   {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MET [GeV]",               "Norm Events"));
    vh.push_back(PHS("ht",          {PDC("single", "ht",                    {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",             "Norm Events"));
    vh.push_back(PHS("nMuons",      {PDC("single", "cutMuVec(size)",        {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(#mu)",                  "Events"));
    vh.push_back(PHS("jetPt",       {PDC("single", "jetsLVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Jet p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("nJets",       {PDC("single", {{"cntNJetsPt30Eta24", dsDY_ll}, {"cntNJetsPt30Eta24", dsDY_nunu}, {"cntNJetsPt30Eta24Zinv", dsDY_ll}})}, {1, 2}, "", 20,  0, 20,    true,  false, "N(jets)",  "Events"));
    vh.push_back(PHS("nJetRemoved", {PDC("single", "nJetsRemoved", {{dsDY_ll}, {dsDY_nunu}})}, {1, 2}, "", 20,  0, 20,    true,  false, "N(jets)",  "Events"));
    vh.push_back(PHS("muPt",        {PDC("single", "cutMuVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("cleanJetPt",  {PDC("single", "cleanJetVec(pt)",       {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("muMtw",       {PDC("single", "muonsMtw",              {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 600,   true,  false, "#mu M_{T}(W) [GeV]",      "Events"));
    vh.push_back(PHS("nMuGenMatch", {PDC("single", "genMatchMuInAcc(size)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(#mu) gen Matched",      "Events"));
    vh.push_back(PHS("genZPt",      {PDC("single", "genZPt",                {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})}, {4, 2}, "", 100, 0, 1500,  true,  true,  "gen Z p_{T} [GeV]",       "Norm Events"));
    vh.push_back(PHS("genZPt_gr",   {PDC("single", "genZPt",                {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})}, {4, 2}, "", 200, 0, 3000,  true,  true,  "gen Z p_{T} [GeV]",       "Norm Events"));
    //vh.push_back(PHS("genMuPt",     {PDC("single", "genMu(pt)",             {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "gen #mu p_{T} [GeV]",     "Norm Events"));
    vh.push_back(PHS("genMuEta",    {PDC("single", "genMu(eta)",            {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3,    true,  true,  "gen #mu #eta",            "Norm Events"));
    vh.push_back(PHS("genZEta",     {PDC("single", "genZEta",               {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3,    false, true,  "gen Z #eta",              "Norm Events"));
    vh.push_back(PHS("genZmass",    {PDC("single", "genZmass",              {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1000,  true,  true,  "gen M(Z) [GeV]",          "Norm Events"));

    // Simple Reco vs. Gen plots
    vh.push_back(PHS("bestRecoZPt"  , {PDC("single",  {{"bestRecoZPt", dsDY_ll}, {"genZPt", dsDY_ll_gen}})},    {1, 2}, "", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]",           "Norm Events"));    
    vh.push_back(PHS("bestRecoZMass", {PDC("single",  {{"bestRecoZM", dsDY_ll},  {"genZmass", dsDY_ll_gen}})},  {1, 2}, "", 100, 0, 1000,  true,  true,  "M(Z) [GeV]",              "Norm Events"));    

    vh.push_back(PHS("zPtRes",             {PDC("single", "ZPtRes",       {dsDY_ll, dsDY_ll_base})}, {1, 1}, "",           100,  -1,  1,  false,  true,  "Z p_{T} (reco - gen)/gen",   "Norm Events"));
    vh.push_back(PHS("zPtRes_Zpt_lt_200",  {PDC("single", "ZPtRes",       {dsDY_ll, dsDY_ll_base})}, {1, 1}, "genZPt<200", 100,  -1,  1,  false,  true,  "Z p_{T} (reco - gen)/gen",   "Norm Events"));
    vh.push_back(PHS("zPtRes_Zpt_gt_200",  {PDC("single", "ZPtRes",       {dsDY_ll, dsDY_ll_base})}, {1, 1}, "genZPt>200", 100,  -1,  1,  false,  true,  "Z p_{T} (reco - gen)/gen",   "Norm Events"));
    vh.push_back(PHS("zEtaRes",            {PDC("single", "ZEtaRes",      {dsDY_ll, dsDY_ll_base})}, {1, 1}, "",           100,  -5,  5,  false,  true,  "Z #eta (reco - gen)",        "Norm Events"));
    vh.push_back(PHS("zPhiRes",            {PDC("single", "ZPhiRes",      {dsDY_ll, dsDY_ll_base})}, {1, 1}, "",           100,  -5,  5,  false,  true,  "Z #phi (reco - gen)",        "Norm Events"));
    vh.push_back(PHS("zMRes",              {PDC("single", "ZMRes",        {dsDY_ll, dsDY_ll_base})}, {1, 1}, "",           100,  -1,  1,  false,  true,  "M(Z) (reco - gen)/gen",      "Norm Events"));

    vh.push_back(PHS("mu1dRMin",              {PDC("single", "mu1dRMin",        {dsDY_ll,        dsDY_ll_base})}, {1, 1}, "",           100,   0,  5,  false,  true,  "min dR (#mu1,jet)",        "Norm Events"));
    vh.push_back(PHS("mu2dRMin",              {PDC("single", "mu2dRMin",        {dsDY_ll,        dsDY_ll_base})}, {1, 1}, "",           100,   0,  5,  false,  true,  "min dR (#mu2,jet)",        "Norm Events"));
    vh.push_back(PHS("mudR",                  {PDC("single", "mudR",            {dsDY_ll,        dsDY_ll_base})}, {1, 1}, "",           100,   0,  5,  false,  true,  "min dR (#mu1,#mu2)",       "Norm Events"));
    vh.push_back(PHS("genMudR",               {PDC("single", "genMudR",         {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "",           100,   0,  5,  false,  true,  "gen min dR (#mu1,#mu2)",   "Norm Events"));
    vh.push_back(PHS("genMudR_Zpt_gt_800",    {PDC("single", "genMudR",         {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "genZPt>800", 100,   0,  5,  false,  true,  "gen min dR (#mu1,#mu2)",   "Norm Events"));
    vh.push_back(PHS("genMudPhi",             {PDC("single", "genMudPhi",       {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "",           100,  -4,  4,  false,  true,  "gen min dPhi (#mu1,#mu2)", "Norm Events"));
    vh.push_back(PHS("genMudEta",             {PDC("single", "genMudEta",       {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "",           100,  -5,  5,  false,  true,  "gen min dEta (#mu1,#mu2)", "Norm Events"));

    // Here are the interesting plots for the closure test 

    // DataCollections which are reused multiple times are defined "upfront" to save time 
    Plotter::DataCollection cleanht(     "single", {{"cleanHt", dsDY_nunu},            {"cleanHt", dsDY_ll_zAcc},               {"cleanHt",               dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanmht(    "single", {{"cleanMHt", dsDY_nunu},           {"cleanMHt", dsDY_ll_zAcc},              {"cleanMHt",              dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu},        {"cleanMHtPhi", dsDY_ll_zAcc},           {"cleanMHtPhi",           dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanMet(    "single", {{"cleanMetPt", dsDY_nunu},         {"cleanMetPt", dsDY_ll_zAcc},            {"cleanMetPt",            dsDY_ll_trig_scaled} });
    Plotter::DataCollection nCleanJet(   "single", {{"cntNJetsPt30Eta24", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_trig_scaled} });
    Plotter::DataCollection nTop(        "single", {{"nTopCandSortedCnt", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc}, {"nTopCandSortedCntZinv", dsDY_ll_trig_scaled} });
    Plotter::DataCollection nBottom(     "single", {{"cntCSVS", dsDY_nunu},            {"cntCSVSZinv", dsDY_ll_zAcc},           {"cntCSVSZinv",           dsDY_ll_trig_scaled} });

    Plotter::DataCollection scaled_ht(   "single", "ht",  {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_mht(  "single", "mht", {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_met(  "single", "met", {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });

    Plotter::DataCollection scaled_cleanht(     "single", {{"cleanHt", dsDY_nunu},               {"cleanHt", dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_cleanmht(    "single", {{"cleanMHt", dsDY_nunu},              {"cleanMHt", dsDY_ll_zAcc_scaled}              });
    Plotter::DataCollection scaled_cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu},           {"cleanMHtPhi", dsDY_ll_zAcc_scaled}           });
    Plotter::DataCollection scaled_cleanMet(    "single", {{"cleanMetPt", dsDY_nunu},            {"cleanMetPt", dsDY_ll_zAcc_scaled}            });
    Plotter::DataCollection scaled_nCleanJet(   "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc_scaled} });
    Plotter::DataCollection scaled_nTop(        "single", {{"nTopCandSortedCntZinv", dsDY_nunu}, {"nTopCandSortedCntZinv", dsDY_ll_zAcc_scaled} });
    Plotter::DataCollection scaled_nBottom(     "single", {{"cntCSVSZinv", dsDY_nunu},           {"cntCSVSZinv", dsDY_ll_zAcc_scaled}           });

    PDC trg_cleanht(  "single","cleanHt",              {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_cleanmht( "single","cleanMHt",             {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_cleanMet( "single","cleanMetPt",           {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nCleanJet("single","cntNJetsPt30Eta24Zinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nTop(     "single","nTopCandSortedCntZinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nBottom(  "single","cntCSVSZinv",          {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});

    Plotter::DataCollection scaled_stacked_DYtoll_genht(  "stack", "genHt",    {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

    Plotter::DataCollection scaled_stacked_DYtoll_cleanht(  "stack", "cleanHt",    {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_ht(       "stack", "ht",         {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_cleanmht( "stack", "cleanMHt",   {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_cleanMet( "stack", "cleanMetPt", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

    Plotter::DataCollection cleanht_znunu(     "single", {{"cleanHt",    dsDY_nunu}});
    Plotter::DataCollection ht_znunu(          "single", {{"ht",         dsDY_nunu}});
    Plotter::DataCollection cleanmht_znunu(    "single", {{"cleanMHt",   dsDY_nunu}});
    Plotter::DataCollection mht_znunu(         "single", {{"mht",        dsDY_nunu}});
    Plotter::DataCollection cleanMet_znunu(    "single", {{"cleanMetPt", dsDY_nunu}});
    Plotter::DataCollection met_znunu(         "single", {{"met",        dsDY_nunu}});

    Plotter::DataCollection scaled_stacked_DYtonunu_cleanht(  "stack", "cleanHt",    {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtonunu_cleanmht( "stack", "cleanMHt",   {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtonunu_cleanMet( "stack", "cleanMetPt", {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});

    vh.push_back(PHS("trigStudy_mu1Pt", {PDC("single", "cutMuPt1",  {dsDY_ll_baseNT, dsDY_ll})}, {2, 1}, "", 100, 0, 1500,  true,  false,  "#mu p_{T} [GeV]", "Events"));
    vh.push_back(PHS("trigStudy_mu2Pt", {PDC("single", "cutMuPt2",  {dsDY_ll_baseNT, dsDY_ll})}, {2, 1}, "", 100, 0, 1500,  true,  false,  "#mu p_{T} [GeV]", "Events"));

    vh.push_back(PHS("trigStudy_cleanht",   {trg_cleanht  },   {1, 1}, "",     50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigStudy_cleanmht",  {trg_cleanmht },   {1, 1}, "",     50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigStudy_cleanmet",  {trg_cleanMet },   {1, 1}, "",     50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigStudy_nTop",      {trg_nTop     },   {1, 1}, "",     10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigStudy_cleannJet", {trg_nCleanJet},   {1, 1}, "",     20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigStudy_nBottom",   {trg_nBottom  },   {1, 1}, "",     10,  0,       10, true,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigStudy_baseline_cleanht",   {trg_cleanht  },   {1, 1}, "passZinvBaseline",     50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleanmht",  {trg_cleanmht },   {1, 1}, "passZinvBaseline",     50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleanmet",  {trg_cleanMet },   {1, 1}, "passZinvBaseline",     50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigStudy_baseline_nTop",      {trg_nTop     },   {1, 1}, "passZinvBaseline",     10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleannJet", {trg_nCleanJet},   {1, 1}, "passZinvBaseline",     20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigStudy_baseline_nBottom",   {trg_nBottom  },   {1, 1}, "passZinvBaseline",     10,  0,       10, true,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigStudy_baselineNoTag_cleanht",   {trg_cleanht  },   {1, 1}, "passZinvBaselineNoTag",     50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleanmht",  {trg_cleanmht },   {1, 1}, "passZinvBaselineNoTag",     50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleanmet",  {trg_cleanMet },   {1, 1}, "passZinvBaselineNoTag",     50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_nTop",      {trg_nTop     },   {1, 1}, "passZinvBaselineNoTag",     10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleannJet", {trg_nCleanJet},   {1, 1}, "passZinvBaselineNoTag",     20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_nBottom",   {trg_nBottom  },   {1, 1}, "passZinvBaselineNoTag",     10,  0,       10, true,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_singleMu45_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_singleMu45_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_singleMu45_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_singleMu45_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_singleMu45_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_singleMu45_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_30_20_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_30_20_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_30_20_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_30_20_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_30_20_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_30_20_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu30_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_30_30_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_30_30_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_30_30_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_30_30_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_30_30_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_30_30_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu30_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_40_30_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_40_30_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_40_30_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_40_30_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_40_30_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_40_30_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu40_mu30, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_35_20_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_35_20_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_35_20_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_35_20_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_35_20_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_35_20_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu35_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_40_20_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_40_20_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_40_20_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_40_20_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_40_20_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_40_20_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu40_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_45_20_cleanht",   {PDC("ratio", "cleanHt",               {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigRatio_45_20_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigRatio_45_20_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  50,  0, 1500, false,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigRatio_45_20_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigRatio_45_20_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  20,  0,   20, false,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigRatio_45_20_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_mu45_mu20, dsDY_ll_zAcc})},   {1, 1}, "passZinvBaselineNoTag",  10,  0,   10, false,  false,  "N(b)",           "Events"));


    vh.push_back(PHS("genht_DY_mumu_stack",           {scaled_stacked_DYtoll_genht},    {1, 1}, "",                                    100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_nr0",     {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "nJetsRemoved=0",   100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_nr1",     {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "nJetsRemoved=1",   100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_nr2",     {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "nJetsRemoved=2",   100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_nr3",     {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "nJetsRemoved>2",   100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_er",      {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "",                 150, 0,     3000, true,  false,  "H_{T} [GeV]",    "Events"));

    vh.push_back(PHS("cleanht_DY_mumu_stack",           {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("ht_DY_mumu_stack",                {scaled_stacked_DYtoll_ht,       ht_znunu      },    {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht_DY_mumu_stack",          {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "",                                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet_DY_mumu_stack",          {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "",                                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("baseline_cleanht_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "passZinvBaseline",                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanmht_DY_mumu_stack", {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "passZinvBaseline",                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmet_DY_mumu_stack", {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "passZinvBaseline",                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
                                                                                                       
    vh.push_back(PHS("nb_0_cleanht_DY_mumu_stack",      {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanmht_DY_mumu_stack",     {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmet_DY_mumu_stack",     {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("cleanht_DY_nunu_stack",           {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht_DY_nunu_stack",          {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "",                                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet_DY_nunu_stack",          {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "",                                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("baseline_cleanht_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "passZinvBaseline",                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanmht_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passZinvBaseline",                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmet_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passZinvBaseline",                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("nb_0_cleanht_DY_nunu_stack",      {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanmht_DY_nunu_stack",     {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmet_DY_nunu_stack",     {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("norm_cleanht",                 {cleanht},      {3, 1}, "",                                                             100, 0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_baseline_cleanht",        {cleanht},      {3, 1}, "passZinvBaseline",                                             50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_0_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_1_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_2_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_3_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nt_0_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nt_1_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nt_2_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nt_3_cleanht",            {cleanht},      {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_0_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_0_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_0_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_0_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_1_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_1_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_1_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_1_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_2_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_2_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_2_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_2_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_3_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_3_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_3_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_3_cleanht",       {cleanht},      {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  true,   "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("norm_cleanmht",                {cleanmht},     {3, 1}, "",                                                             100, 0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_baseline_cleanmht",       {cleanmht},     {3, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_0_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_1_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_2_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_3_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nt_0_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nt_1_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nt_2_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nt_3_cleanmht",           {cleanmht},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_0_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_0_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_0_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_0_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_1_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_1_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_1_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_1_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_2_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_2_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_2_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_2_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_3_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_3_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_3_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_3_cleanmht",      {cleanmht},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  true,   "MH_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("norm_cleanmet",                {cleanMet},     {3, 1}, "",                                                             100, 0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_baseline_cleanmet",       {cleanMet},     {3, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_0_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_1_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_2_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_3_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nt_0_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nt_1_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nt_2_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nt_3_cleanmet",           {cleanMet},     {3, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_0_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_0_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_0_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_0_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_1_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_1_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_1_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_1_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_2_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_2_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_2_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_2_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_0_nt_3_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_1_nt_3_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_2_nt_3_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nb_3_nt_3_cleanmet",      {cleanMet},     {3, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  true,   "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("norm_nTop",                    {nTop},         {3, 1}, "",                                                             10,  0,       10, true,  true,   "N(t)",           "Norm Events"));
    vh.push_back(PHS("norm_baseline_nTop",           {nTop},         {3, 1}, "passZinvBaseline",                                             10,  0,       10, true,  true,   "N(t)",           "Norm Events"));
    vh.push_back(PHS("norm_cleannJet",               {nCleanJet},    {3, 1}, "",                                                             20,  0,       20, true,  true,   "N(jet)",         "Norm Events"));
    vh.push_back(PHS("norm_baseline_cleannJet",      {nCleanJet},    {3, 1}, "passZinvBaseline",                                             20,  0,       20, true,  true,   "N(jet)",         "Norm Events"));
    vh.push_back(PHS("norm_nBottom",                 {nBottom},      {3, 1}, "",                                                             10,  0,       10, true,  true,   "N(b)",           "Norm Events"));
    vh.push_back(PHS("norm_baseline_nBottom",        {nBottom},      {3, 1}, "passZinvBaseline",                                             10,  0,       10, true,  true,   "N(b)",           "Norm Events"));
    vh.push_back(PHS("norm_cleanmhtphi",             {cleanmhtphi},  {3, 1}, "",                                                             100, -3.14, 3.14, false, true,   "#phi(MH_{T})",   "Norm Events"));

    vh.push_back(PHS("cleanht",                {scaled_cleanht},     {2, 1}, "",                                                             100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanht",       {scaled_cleanht},     {2, 1}, "passZinvBaseline",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baselineNoTag_cleanht",  {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag",                                        50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_ht",            {scaled_ht},          {2, 1}, "passZinvBaseline",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_0_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_1_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_2_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_3_cleanht",           {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanht",      {scaled_cleanht},     {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht",               {scaled_cleanmht},    {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmht",      {scaled_cleanmht},    {2, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baselineNoTag_cleanmht", {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag",                                        50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_mht",           {scaled_mht},         {2, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_0_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_1_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_2_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_3_cleanmht",          {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanmht",     {scaled_cleanmht},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet",               {scaled_cleanMet},    {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("cleanmet_gr",            {scaled_cleanMet},    {2, 1}, "",                                                             200, 0,     3000, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("baseline_cleanmet",      {scaled_cleanMet},    {2, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("baselineNoTag_cleanmet", {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag",                                        50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("baseline_met",           {scaled_met},         {2, 1}, "passZinvBaseline",                                             50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_0_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=0",                100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_1_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_2_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_3_cleanmet",          {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=0",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=0",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=0",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=0",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanmet",     {scaled_cleanMet},    {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nTop",                   {scaled_nTop},        {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("baseline_nTop",          {scaled_nTop},        {2, 1}, "passZinvBaseline",                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_0_nTop",              {scaled_nTop},        {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=0",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_1_nTop",              {scaled_nTop},        {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=1",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_2_nTop",              {scaled_nTop},        {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv=2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_3_nTop",              {scaled_nTop},        {2, 1}, "passZinvBaselineNoTag;cntCSVSZinv>2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("cleannJet",              {scaled_nCleanJet},   {2, 1}, "",                                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("baseline_cleannJet",     {scaled_nCleanJet},   {2, 1}, "passZinvBaseline",                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("nBottom",                {scaled_nBottom},     {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("baseline_nBottom",       {scaled_nBottom},     {2, 1}, "passZinvBaseline",                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("cleanmhtphi",            {scaled_cleanmhtphi}, {2, 1}, "",                                                             100, -3.14, 3.14, false, false,  "#phi(MH_{T})",   "Events"));

    vh.push_back(PHS("muPt_met_lt_150", {PDC("single", "cutMuVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "cleanMetPt<150", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("cleanht_met_lt_150",   {scaled_cleanht}, {1, 1}, "cleanMetPt<150", 50, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmet_met_lt_150", {scaled_cleanMet}, {1, 1}, "cleanMetPt<150", 50, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("cleanmht_met_lt_150", {scaled_cleanmht}, {1, 1}, "cleanMetPt<150", 50, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));

    Plotter::DataCollection  cleanht_DY_ll_nb("single", "cleanHt",    {dsDY_ll_0b, dsDY_ll_1b, dsDY_ll_2b, dsDY_ll_3b});
    Plotter::DataCollection cleanmht_DY_ll_nb("single", "cleanMHt",   {dsDY_ll_0b, dsDY_ll_1b, dsDY_ll_2b, dsDY_ll_3b});
    Plotter::DataCollection cleanMet_DY_ll_nb("single", "cleanMetPt", {dsDY_ll_0b, dsDY_ll_1b, dsDY_ll_2b, dsDY_ll_3b});

    Plotter::DataCollection  cleanht_DY_ll_nt("single", "cleanHt",    {dsDY_ll_0t, dsDY_ll_1t, dsDY_ll_2t, dsDY_ll_3t});
    Plotter::DataCollection cleanmht_DY_ll_nt("single", "cleanMHt",   {dsDY_ll_0t, dsDY_ll_1t, dsDY_ll_2t, dsDY_ll_3t});
    Plotter::DataCollection cleanMet_DY_ll_nt("single", "cleanMetPt", {dsDY_ll_0t, dsDY_ll_1t, dsDY_ll_2t, dsDY_ll_3t});

    Plotter::DataCollection  cleanht_DY_nunu_nb("single", "cleanHt",    {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b});
    Plotter::DataCollection cleanmht_DY_nunu_nb("single", "cleanMHt",   {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b});
    Plotter::DataCollection cleanMet_DY_nunu_nb("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b});

    Plotter::DataCollection  cleanht_DY_nunu_nt("single", "cleanHt",    {dsDY_nunu_0t, dsDY_nunu_1t, dsDY_nunu_2t, dsDY_nunu_3t});
    Plotter::DataCollection cleanmht_DY_nunu_nt("single", "cleanMHt",   {dsDY_nunu_0t, dsDY_nunu_1t, dsDY_nunu_2t, dsDY_nunu_3t});
    Plotter::DataCollection cleanMet_DY_nunu_nt("single", "cleanMetPt", {dsDY_nunu_0t, dsDY_nunu_1t, dsDY_nunu_2t, dsDY_nunu_3t});

    PDC cleanMet_DY_ll_nb_1_0("ratio", "cleanMetPt", {dsDY_nunu_1b, dsDY_nunu_0b});
    PDC cleanMet_DY_ll_nb_2_0("ratio", "cleanMetPt", {dsDY_nunu_2b, dsDY_nunu_0b});
    PDC cleanMet_DY_ll_nb_3_0("ratio", "cleanMetPt", {dsDY_nunu_3b, dsDY_nunu_0b});

    vh.push_back(PHS("DY_mumu_nb_Ratio_cleanmet",               {cleanMet_DY_ll_nb_1_0, cleanMet_DY_ll_nb_2_0, cleanMet_DY_ll_nb_3_0}, {1, 1}, "",                      50, 0,     1500, false,  true,  "MET [GeV]", "Norm Ratio"));
    vh.push_back(PHS("DY_mumu_nb_baseline_Ratio_cleanmet",      {cleanMet_DY_ll_nb_1_0, cleanMet_DY_ll_nb_2_0, cleanMet_DY_ll_nb_3_0}, {1, 1}, "passZinvBaseline",      50, 0,     1500, false,  true,  "MET [GeV]", "Norm Ratio"));
    vh.push_back(PHS("DY_mumu_nb_baselineNoTag_Ratio_cleanmet", {cleanMet_DY_ll_nb_1_0, cleanMet_DY_ll_nb_2_0, cleanMet_DY_ll_nb_3_0}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, false,  true,  "MET [GeV]", "Norm Ratio"));

    vh.push_back(PHS("DY_mumu_nb_cleanht",   {cleanht_DY_ll_nb}, {1, 1}, "", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_mumu_nb_cleanmet", {cleanMet_DY_ll_nb}, {1, 1}, "", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_mumu_nb_cleanmht", {cleanmht_DY_ll_nb}, {1, 1}, "", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_mumu_nt_cleanht",   {cleanht_DY_ll_nt}, {1, 1}, "", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_mumu_nt_cleanmet", {cleanMet_DY_ll_nt}, {1, 1}, "", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_mumu_nt_cleanmht", {cleanmht_DY_ll_nt}, {1, 1}, "", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_nunu_nb_cleanht",   {cleanht_DY_nunu_nb}, {1, 1}, "", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_nunu_nb_cleanmet", {cleanMet_DY_nunu_nb}, {1, 1}, "", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_nunu_nb_cleanmht", {cleanmht_DY_nunu_nb}, {1, 1}, "", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_nunu_nt_cleanht",   {cleanht_DY_nunu_nt}, {1, 1}, "", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_nunu_nt_cleanmet", {cleanMet_DY_nunu_nt}, {1, 1}, "", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_nunu_nt_cleanmht", {cleanmht_DY_nunu_nt}, {1, 1}, "", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_mumu_nb_baseline_cleanht",   {cleanht_DY_ll_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_mumu_nb_baseline_cleanmet", {cleanMet_DY_ll_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_mumu_nb_baseline_cleanmht", {cleanmht_DY_ll_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_mumu_nt_baseline_cleanht",   {cleanht_DY_ll_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_mumu_nt_baseline_cleanmet", {cleanMet_DY_ll_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_mumu_nt_baseline_cleanmht", {cleanmht_DY_ll_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_nunu_nb_baseline_cleanht",   {cleanht_DY_nunu_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_nunu_nb_baseline_cleanmet", {cleanMet_DY_nunu_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_nunu_nb_baseline_cleanmht", {cleanmht_DY_nunu_nb}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));

    vh.push_back(PHS("DY_nunu_nt_baseline_cleanht",   {cleanht_DY_nunu_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     2000, true,  true,  "H_{T} [GeV]",    "Norm Events"));
    vh.push_back(PHS("DY_nunu_nt_baseline_cleanmet", {cleanMet_DY_nunu_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MET [GeV]",      "Norm Events"));
    vh.push_back(PHS("DY_nunu_nt_baseline_cleanmht", {cleanmht_DY_nunu_nt}, {1, 1}, "passZinvBaselineNoTag", 50, 0,     1500, true,  true,  "MH_{T} [GeV]",   "Norm Events"));


    // Plots the information needed to calculate acceptance corrections, plot itself is meaningless
    vh.push_back(PHS("accInfo",  {PDC("single", {{"genMuInAcc(size)", dsDY_ll_forEff}, {"genMu(size)", dsDY_ll_forEff}})},    {1, 2}, "",  20,  0,       20, true, false,  "gen N(#mu)",     "Events"));

    vh.push_back(PHS("genMuPt",  {PDC("single", {{"genMuInAcc(pt)", dsDY_ll_forEff},            {"genMu(pt)", dsDY_ll_forEff}})},             {1, 1}, "",  150,  0,  1500, true,  false,  "gen #mu p_{T} [GeV]", "Events"));
    vh.push_back(PHS("accHt",    {PDC("ratio",  {{"cleanHt", dsDY_ll_forEff_den},               {"cleanHt", dsDY_ll_forAcc}})},               {1, 1}, "",  300,  0,  3000, false, false,  "H_{T} [GeV]",         "Acc"));
    vh.push_back(PHS("accMt",    {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_den}, {"nTopCandSortedCntZinv", dsDY_ll_forAcc}})}, {1, 1}, "",  10,   0,  10,   false, false,  "N(t)",                "Acc"));
    vh.push_back(PHS("accNb",    {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_den},           {"cntCSVSZinv", dsDY_ll_forAcc}})},           {1, 1}, "",  10,   0,  10,   false, false,  "N(b)",                "Acc"));
    vh.push_back(PHS("accNjet",  {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forAcc}})}, {1, 1}, "",  20,   0,  20,   false, false,  "N(j)",                "Acc"));

    // Plots for calculation of muon efficiency 
    vh.push_back(PHS("effMuPt",      {PDC("ratio",  {{"genMatchMuInAcc(pt)", dsDY_ll_forEff},       {"genMuInAcc(pt)", dsDY_ll_forEff}})},             {1, 1}, "", 200, 0,     2000, false, false,  "#mu p_{T} [GeV]",  "Efficiency"));
    vh.push_back(PHS("effMuEta",     {PDC("ratio",  {{"genMatchMuInAcc(eta)", dsDY_ll_forEff},      {"genMuInAcc(eta)", dsDY_ll_forEff}})},            {1, 1}, "", 100, -3.0,  3.0,  false, false,  "#mu #eta",         "Efficiency"));
    vh.push_back(PHS("effMucleanHt", {PDC("ratio",  {{"cleanHt", dsDY_ll_forEff_num},               {"cleanHt", dsDY_ll_forEff_den}})},                {1, 1}, "", 300, 0,     3000, false, false,  "H_{T} [GeV]",      "Efficiency"));
    vh.push_back(PHS("effMuNb",      {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_num},           {"cntCSVSZinv", dsDY_ll_forEff_den}})},            {1, 1}, "", 10, 0,      10.0, false, false,  "N(b)",             "Efficiency"));
    vh.push_back(PHS("effMuNt",      {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_num}, {"nTopCandSortedCntZinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 10, 0,      10.0, false, false,  "N(t)",             "Efficiency"));
    vh.push_back(PHS("effMuNjet",    {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_num}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 20, 0,      20.0, false, false,  "N(j)",             "Efficiency"));

    // Plots of weight and efficiencies used as sanity checks
    vh.push_back(PHS("zEffWgt", {PDC("single", "zEffWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Eff Weight",         ""));
    vh.push_back(PHS("zEff",    {PDC("single", "zEff",     {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0,   false, false,  "Z Eff Weight",         ""));
    vh.push_back(PHS("zAccWgt", {PDC("single", "zAccWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Acc Weight",         ""));
    vh.push_back(PHS("zAcc",    {PDC("single", "zAcc",     {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0,   false, false,  "Z Acc Weight",         ""));

    // Test plots, these are simply meaningless demonstrations of features
    Plotter::DataCollection dcDY_test(   "single", {{"mht", dsDY_ll}, {"jetsLVec(pt)", dsDY_nunu}, {"genZmass", dsDY_test}, {"genZPt", dsDY_test}});
    Plotter::DataCollection dcDY_ratio(  "ratio",  {{"cleanJetVec(pt)", dsDY_ll}, {"jetsLVec(pt)", dsDY_ll}});
    Plotter::DataCollection dcDY_stack(  "stack",  "mht", {dsDY_ll, dsDY_nunu});

    vh.push_back(PHS("test",  {dcDY_stack, dcDY_test},  {1, 2}, "genZmass>80", 100, 0, 1000,   true, true,  "???",         "Norm Events"));
    vh.push_back(PHS("test2", {dcDY_ratio},  {1, 1}, "", 100, 0, 500, true, false, "jpt???", "Ratio"));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles);
    if(doSave)  plotter.saveHists();
    if(doPlots) plotter.plot();
}
