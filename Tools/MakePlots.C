#include "Plotter.h"
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
        {"numFiles",  required_argument, 0, 'N'},
        {"startFile", required_argument, 0, 'M'},
        {"numEvts",   required_argument, 0, 'E'},
        {"plotDir",   required_argument, 0, 'P'},
	{"luminosity",      required_argument, 0, 'L'}
    };

    bool doPlots = true, doSave = true, fromTuple = true, runOnCondor = false;
    string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi = AnaSamples::luminosity;

    while((opt = getopt_long(argc, argv, "psfcH:D:N:M:E:P:L:", long_options, &option_index)) != -1)
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

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            nEvts = int(atoi(optarg));
            break;

        case 'P':
            plotDir = optarg;
            break;

	case 'L':
            lumi = atof(optarg); 
            break;
        }
    }

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "histoutput_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        doSave = true;
        doPlots = false;
        fromTuple = true;
        //sampleloc = "";
    }

    AnaSamples::SampleSet        ss(sampleloc, lumi);
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
        else if(dataSets.compare("ZJetsToNuNu") == 0)
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

    Plotter::DatasetSummary dsDY_ll_noSel(      "DY#rightarrow#mu#mu no sel",            fileMap["DYJetsToLL"],   "pdgIdZDec=13", "");
    Plotter::DatasetSummary dsDY_ll_gen_noSel(  "DY#rightarrow#mu#mu Gen no Sel",        fileMap["DYJetsToLL"],   "pdgIdZDec=13", "");

    Plotter::DatasetSummary dsDY_ll(        "DY#rightarrow#mu#mu",                   fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_ll_base(   "DY#rightarrow#mu#mu baseline",          fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passBaselineZinv", "");
    Plotter::DatasetSummary dsDY_ll_baseNT( "DY#rightarrow#mu#mu baselineNoTag",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsDY_ll_zWeight("DY#rightarrow#mu#mu no #mu, Z eff",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt");
    Plotter::DatasetSummary dsDY_ll_zAcc(   "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_zAccAcc("DY#rightarrow#mu#mu no #mu, Z acc",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;ngenMuInAcc>1.5", "zAccWgt");
    Plotter::DatasetSummary dsDY_ll_zAccAcc_nw("DY#rightarrow#mu#mu no #mu, Z acc NW",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;ngenMuInAcc>1.5", "");
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
    Plotter::DatasetSummary dsDY_ll_scaled_baseNT( "DY#rightarrow#mu#mu baselineNoTag",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passBaselineNoTagZinv", "", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_zWeight_scaled("DY#rightarrow#mu#mu no #mu, Z eff",     fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt",         znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_zAcc_scaled(   "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_trig_scaled(   "DY#rightarrow#mu#mu no #mu, Z eff+acc+trig", fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel;passDiMuIsoTrig", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_llclean_scaled(   "DY#rightarrow#mu#mu no #mu",            fileMap["DYJetsToLL"],   "pdgIdZDec=13;passMuZinvSel", "",                znunu_mumu_ratio / zAcc);

    Plotter::DatasetSummary dsDY_ll_forAcc(     "DY#rightarrow#mu#mu no Sel forAcc",                      fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMu");
    Plotter::DatasetSummary dsDY_ll_forEff(     "DY#rightarrow#mu#mu no Sel forEff",                      fileMap["DYJetsToLL"],  "pdgIdZDec=13", "");
    Plotter::DatasetSummary dsDY_ll_forEff_num( "DY#rightarrow#mu#mu no Sel N(#mu) match wgt",     fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMatchMuInAcc");
    Plotter::DatasetSummary dsDY_ll_forEff_den( "DY#rightarrow#mu#mu no Sel N(#mu) wgt",           fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMuInAcc");

    Plotter::DatasetSummary dsDY_nunu(          "Z#rightarrow#nu#nu",                fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_SB0b(   "Z#rightarrow#nu#nu, N(b) = 0",        fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_SBMb(   "Z#rightarrow#nu#nu, Direct MC",       fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_test(          "test",                              fileMap["ZJetsToNuNu"], "passLeptVeto", "");

    Plotter::DatasetSummary dsDY_ll_0t(        "DY#rightarrow#mu#mu N(t) = 0",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=0", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_1t(        "DY#rightarrow#mu#mu N(t) = 1",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=1", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_2t(        "DY#rightarrow#mu#mu N(t) = 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=2", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_3t(        "DY#rightarrow#mu#mu N(t) > 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv>2", "zEffWgt;zAccWgt");

    Plotter::DatasetSummary dsDY_ll_0b(        "DY#rightarrow#mu#mu N(b) = 0",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=0", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_1b(        "DY#rightarrow#mu#mu N(b) = 1",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=1", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_2b(        "DY#rightarrow#mu#mu N(b) = 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=2", "zEffWgt;zAccWgt");
    Plotter::DatasetSummary dsDY_ll_3b(        "DY#rightarrow#mu#mu N(b) > 2",       fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv>2", "zEffWgt;zAccWgt");

    Plotter::DatasetSummary dsDY_nunu_0t(          "Z#rightarrow#nu#nu N(t) = 0",    fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=0", "");
    Plotter::DatasetSummary dsDY_nunu_1t(          "Z#rightarrow#nu#nu N(t) = 1",    fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=1", "");
    Plotter::DatasetSummary dsDY_nunu_2t(          "Z#rightarrow#nu#nu N(t) = 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=2", "");
    Plotter::DatasetSummary dsDY_nunu_3t(          "Z#rightarrow#nu#nu N(t) > 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv>2", "");
                                        
    Plotter::DatasetSummary dsDY_nunu_0b(          "Z#rightarrow#nu#nu N(b) = 0",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_nunu_1b(          "Z#rightarrow#nu#nu N(b) = 1",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1", "");
    Plotter::DatasetSummary dsDY_nunu_2b(          "Z#rightarrow#nu#nu N(b) = 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=2", "");
    Plotter::DatasetSummary dsDY_nunu_3b(          "Z#rightarrow#nu#nu N(b) > 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv>2", "");

    Plotter::DatasetSummary dsDY_nunu_0b_blNoTag(          "Z#rightarrow#nu#nu N(b) = 0",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsDY_nunu_1b_blNoTag(          "Z#rightarrow#nu#nu N(b) = 1",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsDY_nunu_2b_blNoTag(          "Z#rightarrow#nu#nu N(b) = 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=2;passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsDY_nunu_3b_blNoTag(          "Z#rightarrow#nu#nu N(b) > 2",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv>2;passBaselineNoTagZinv", "");

    Plotter::DatasetSummary dsDY_nunu_0b_1fake(  "Z#rightarrow#nu#nu N(b) = 0, 1 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight1fakeb");
    Plotter::DatasetSummary dsDY_nunu_0b_2fake(  "Z#rightarrow#nu#nu N(b) = 0, 2 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight2fakeb");
    Plotter::DatasetSummary dsDY_nunu_0b_3fake(  "Z#rightarrow#nu#nu N(b) = 0, 3 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight3fakeb");

    Plotter::DatasetSummary dsDY_nunu_0b_1fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 1 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv1b", "weight1fakeb");
    Plotter::DatasetSummary dsDY_nunu_0b_2fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 2 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv2b", "weight2fakeb");
    Plotter::DatasetSummary dsDY_nunu_0b_3fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 3 fake b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv3b", "weight3fakeb");

    Plotter::DatasetSummary dsDY_nunu_highCSV(          "Z#rightarrow#nu#nu csvMax > 0.5",    fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV>0.5", "");
    Plotter::DatasetSummary dsDY_nunu_veryhighCSV1b(          "Z#rightarrow#nu#nu csvMax > 0.8 1b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;maxCSV>0.8", "");
    Plotter::DatasetSummary dsDY_nunu_highCSV1b(          "Z#rightarrow#nu#nu csvMax > 0.5 1b",    fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;maxCSV>0.5", "");
    Plotter::DatasetSummary dsDY_nunu_lowCSV(          "Z#rightarrow#nu#nu csvMax < 0.5",    fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV<0.5", "");
    Plotter::DatasetSummary dsDY_nunu_verylowCSV(          "Z#rightarrow#nu#nu csvMax < 0.4",    fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV<0.4", "");

    Plotter::DatasetSummary dsDY_ll_HT_100to200(  "DY#rightarrow#mu#mu 100 < H_{T} < 200 GeV",    fileMap["DYJetsToLL_HT_100to200"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_200to400(  "DY#rightarrow#mu#mu 200 < H_{T} < 400 GeV",    fileMap["DYJetsToLL_HT_200to400"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_400to600(  "DY#rightarrow#mu#mu 400 < H_{T} < 600 GeV",    fileMap["DYJetsToLL_HT_400to600"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ll_HT_600toInf(  "DY#rightarrow#mu#mu H_{T} > 600 GeV",          fileMap["DYJetsToLL_HT_600toInf"],   "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);

    Plotter::DatasetSummary dsDY_nunu_HT_100to200("DY#rightarrow#nu#nu 100 < H_{T} < 200 GeV",    fileMap["ZJetsToNuNu_HT_100to200"],  "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_200to400("DY#rightarrow#nu#nu 200 < H_{T} < 400 GeV",    fileMap["ZJetsToNuNu_HT_200to400"],  "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_400to600("DY#rightarrow#nu#nu 400 < H_{T} < 600 GeV",    fileMap["ZJetsToNuNu_HT_400to600"],  "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_HT_600toInf("DY#rightarrow#nu#nu H_{T} > 600 GeV",          fileMap["ZJetsToNuNu_HT_600toInf"],  "passLeptVeto", "");

    Plotter::DatasetSummary dsDY_ll_passMuonVeto(     "passMuonVeto",         fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuonVeto", "");
    Plotter::DatasetSummary dsDY_ll_passEleVeto(      "passEleVeto",          fileMap["DYJetsToLL"], "pdgIdZDec=13;passEleVeto", "");
    Plotter::DatasetSummary dsDY_ll_passIsoTrkVeto(   "passIsoTrkVeto",       fileMap["DYJetsToLL"], "pdgIdZDec=13;passIsoTrkVeto", "");
    Plotter::DatasetSummary dsDY_ll_passnJets(        "passnJets",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passnJets", "");
    Plotter::DatasetSummary dsDY_ll_passdPhis(        "passdPhis",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passdPhis", "");
    Plotter::DatasetSummary dsDY_ll_passBJets(        "passBJets",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passBJets", "");
    Plotter::DatasetSummary dsDY_ll_passMET(          "passMET",              fileMap["DYJetsToLL"], "pdgIdZDec=13;passMET", "");
    Plotter::DatasetSummary dsDY_ll_passTagger(       "passTagger",           fileMap["DYJetsToLL"], "pdgIdZDec=13;passTagger", "");
    //Plotter::DatasetSummary dsDY_ll_passNewCuts(      "passNewCuts",          fileMap["DYJetsToLL"], "pdgIdZDec=13;passNewCuts", "");
    Plotter::DatasetSummary dsDY_ll_passBaseline(     "passBaseline",         fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaseline", "");
    Plotter::DatasetSummary dsDY_ll_passBaselineNoTag("passBaselineNoTag",    fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineNoTag", "");

    Plotter::DatasetSummary dsDY_ll_passMuonVetoZinv(     "passMuonVetoZinv",         fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuonVetoZinv", "");
    Plotter::DatasetSummary dsDY_ll_passEleVetoZinv(      "passEleVetoZinv",          fileMap["DYJetsToLL"], "pdgIdZDec=13;passEleVetoZinv", "");
    Plotter::DatasetSummary dsDY_ll_passIsoTrkVetoZinv(   "passIsoTrkVetoZinv",       fileMap["DYJetsToLL"], "pdgIdZDec=13;passIsoTrkVetoZinv", "");
    Plotter::DatasetSummary dsDY_ll_passnJetsZinv(        "passnJetsZinv",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passnJetsZinv", "");
    Plotter::DatasetSummary dsDY_ll_passdPhisZinv(        "passdPhisZinv",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passdPhisZinv", "");
    Plotter::DatasetSummary dsDY_ll_passBJetsZinv(        "passBJetsZinv",            fileMap["DYJetsToLL"], "pdgIdZDec=13;passBJetsZinv", "");
    Plotter::DatasetSummary dsDY_ll_passMETZinv(          "passMETZinv",              fileMap["DYJetsToLL"], "pdgIdZDec=13;passMETZinv", "");
    Plotter::DatasetSummary dsDY_ll_passTaggerZinv(       "passTaggerZinv",           fileMap["DYJetsToLL"], "pdgIdZDec=13;passTaggerZinv", "");
    //Plotter::DatasetSummary dsDY_ll_passNewCutsZinv(      "passNewCutsZinv",          fileMap["DYJetsToLL"], "pdgIdZDec=13;passNewCutsZinv", "");
    Plotter::DatasetSummary dsDY_ll_passBaselineZinv(     "passBaselineZinv",         fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineZinv", "");
    Plotter::DatasetSummary dsDY_ll_passBaselineNoTagZinv("passBaselineNoTagZinv",    fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineNoTagZinv", "");

    Plotter::DataCollection scaled_chargedEMFrac( "single", "cleanChargedHadEFrac", {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_neutralEMFrac( "single", "cleanNeutralEMEFrac",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_chargedHadFrac("single", "cleanChargedEMEFrac",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

    Plotter::DataCollection scaled_chargedEMFrac_j1( "single", "cleanChargedHadEFrac[0]", {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_neutralEMFrac_j1( "single", "cleanNeutralEMEFrac[0]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_chargedHadFrac_j1("single", "cleanChargedEMEFrac[0]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

    Plotter::DataCollection scaled_chargedEMFrac_j2( "single", "cleanChargedHadEFrac[1]", {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_neutralEMFrac_J2( "single", "cleanNeutralEMEFrac[1]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
    Plotter::DataCollection scaled_chargedHadFrac_j2("single", "cleanChargedEMEFrac[1]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

    vector<Plotter::HistSummary> vh;

    //Basic Raw compairson plots  
    vh.push_back(PHS("inc_mht",     {PDC("single", "mht",                   {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "M(H_{t}) [GeV]",          "Events"));
    vh.push_back(PHS("inc_met",     {PDC("single", "met",                   {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "MET [GeV]",               "Events"));
    vh.push_back(PHS("inc_ht",      {PDC("single", "ht",                    {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  "H_{T} [GeV]",             "Events"));
    vh.push_back(PHS("inc_genht",   {PDC("single", "genHt",                 {dsDY_ll_inc, dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "gen H_{T} [GeV]",             "Norm Events"));

    vh.push_back(PHS("minljdR",                {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                      100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));
    vh.push_back(PHS("minljdR_baseline",       {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",      100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));
    vh.push_back(PHS("minljdR_baselineNoTag",  {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));

    vh.push_back(PHS("dPhi1",                   {PDC("single", "dPhiVecZinv[0]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                             100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j1)",   "Norm Events"));
    vh.push_back(PHS("dPhi2",                   {PDC("single", "dPhiVecZinv[1]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                             100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j2)",   "Norm Events"));
    vh.push_back(PHS("dPhi3",                   {PDC("single", "dPhiVecZinv[2]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                             100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j3)",   "Norm Events"));
    vh.push_back(PHS("dPhi1_nJet",              {PDC("single", "dPhiVecZinv[0]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv",                100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j1)",   "Norm Events"));
    vh.push_back(PHS("dPhi2_nJet",              {PDC("single", "dPhiVecZinv[1]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv",                100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j2)",   "Norm Events"));
    vh.push_back(PHS("dPhi3_nJet",              {PDC("single", "dPhiVecZinv[2]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv",                100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j3)",   "Norm Events"));
    vh.push_back(PHS("dPhi1_nJetMET",           {PDC("single", "dPhiVecZinv[0]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;passMETZinv",    100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j1)",   "Norm Events"));
    vh.push_back(PHS("dPhi2_nJetMET",           {PDC("single", "dPhiVecZinv[1]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;passMETZinv",    100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j2)",   "Norm Events"));
    vh.push_back(PHS("dPhi3_nJetMET",           {PDC("single", "dPhiVecZinv[2]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;passMETZinv",    100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j3)",   "Norm Events"));
    vh.push_back(PHS("dPhi1_nJet_met_gt_600",   {PDC("single", "dPhiVecZinv[0]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j1)",   "Norm Events"));
    vh.push_back(PHS("dPhi2_nJet_met_gt_600",   {PDC("single", "dPhiVecZinv[1]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j2)",   "Norm Events"));
    vh.push_back(PHS("dPhi3_nJet_met_gt_600",   {PDC("single", "dPhiVecZinv[2]",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(MET,j3)",   "Norm Events"));
    vh.push_back(PHS("cleanJetPt",                 {PDC("single", "cleanJetVec(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                      100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJetPt_baseline",        {PDC("single", "cleanJetVec(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",      100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJetPt_baselineNoTag",   {PDC("single", "cleanJetVec(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ",               {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                      100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ_baseline",      {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",      100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ_baselineNoTag", {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt",               {PDC("single", "cleanJetVec[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                                    100, 0,     1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt_baseline",      {PDC("single", "cleanJetVec[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",                    100, 0,     1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt_baseline_dPhi", {PDC("single", "cleanJetVec[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv;dPhiVecZinv[0]<0.5", 100, 0,     1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt_baselineNoTag", {PDC("single", "cleanJetVec[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",               100, 0,     1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt",               {PDC("single", "cleanJetVec[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                                    100, 0,     1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt_baseline",      {PDC("single", "cleanJetVec[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",                    100, 0,     1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt_baseline_dPhi", {PDC("single", "cleanJetVec[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv;dPhiVecZinv[0]<0.5", 100, 0,     1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt_baselineNoTag", {PDC("single", "cleanJetVec[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",               100, 0,     1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt",               {PDC("single", "cleanJetVec[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                                    100, 0,     1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt_baseline",      {PDC("single", "cleanJetVec[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",                    100, 0,     1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt_baseline_dPhi", {PDC("single", "cleanJetVec[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv;dPhiVecZinv[0]<0.5", 100, 0,     1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt_baselineNoTag", {PDC("single", "cleanJetVec[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",               100, 0,     1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));

    vh.push_back(PHS("cleanJetChargedEMFrac",              {PDC("single", "cleanChargedEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFrac_j1",           {PDC("single", "cleanChargedEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J1 Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFrac_j2",           {PDC("single", "cleanChargedEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J2 Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFrac_baseline",     {PDC("single", "cleanChargedEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFracJ1_baseline",   {PDC("single", "cleanChargedEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J1 Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFracJ2_baseline",   {PDC("single", "cleanChargedEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J2 Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFrac_baselineNoTag",   {PDC("single", "cleanChargedEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFracJ1_baselineNoTag", {PDC("single", "cleanChargedEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J1 Charged EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedEMFracJ2_baselineNoTag", {PDC("single", "cleanChargedEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J2 Charged EM Frac", "Norm Events"));

    vh.push_back(PHS("cleanJetNeutralEMFrac",              {PDC("single", "cleanNeutralEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFrac_j1",           {PDC("single", "cleanNeutralEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J1 Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFrac_j2",           {PDC("single", "cleanNeutralEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J2 Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFrac_baseline",     {PDC("single", "cleanNeutralEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFracJ1_baseline",   {PDC("single", "cleanNeutralEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J1 Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFracJ2_baseline",   {PDC("single", "cleanNeutralEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J2 Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFrac_baselineNoTag",   {PDC("single", "cleanNeutralEMEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFracJ1_baselineNoTag", {PDC("single", "cleanNeutralEMEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J1 Neutral EM Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetNeutralEMFracJ2_baselineNoTag", {PDC("single", "cleanNeutralEMEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J2 Neutral EM Frac", "Norm Events"));

    vh.push_back(PHS("cleanJetChargedHadFrac",            {PDC("single", "cleanChargedHadEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFrac_j1",         {PDC("single", "cleanChargedHadEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J1 Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFrac_j2",         {PDC("single", "cleanChargedHadEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                 100, 0, 1.0, true, true, "Clean J2 Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFrac_baseline",   {PDC("single", "cleanChargedHadEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFracJ1_baseline", {PDC("single", "cleanChargedHadEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J1 Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFracJ2_baseline", {PDC("single", "cleanChargedHadEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1.0, true, true, "Clean J2 Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFrac_baselineNoTag",   {PDC("single", "cleanChargedHadEFrac",    {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFracJ1_baselineNoTag", {PDC("single", "cleanChargedHadEFrac[0]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J1 Charged Had Frac", "Norm Events"));
    vh.push_back(PHS("cleanJetChargedHadFracJ2_baselineNoTag", {PDC("single", "cleanChargedHadEFrac[1]", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1.0, true, true, "Clean J2 Charged Had Frac", "Norm Events"));

    vh.push_back(PHS("mht",         {PDC("single", "mht",                   {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "M(H_{t}) [GeV]",          "Norm Events"));
    vh.push_back(PHS("mt2",         {PDC("single", "best_had_brJet_MT2Zinv", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MT2 [GeV]",               "Norm Events"));
    vh.push_back(PHS("met",         {PDC("single", "met",                   {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MET [GeV]",               "Norm Events"));
    vh.push_back(PHS("ht",          {PDC("single", "HT",                    {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",             "Norm Events"));
    vh.push_back(PHS("nMuons",      {PDC("single", "cutMuVec(size)",        {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(#mu)",                  "Events"));
    vh.push_back(PHS("jetPt",       {PDC("single", "jetsLVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Jet p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("nJets",       {PDC("single", {{"cntNJetsPt30Eta24", dsDY_ll}, {"cntNJetsPt30Eta24", dsDY_nunu}, {"cntNJetsPt30Eta24Zinv", dsDY_ll}})}, {1, 2}, "", 20,  0, 20,    true,  false, "N(jets)",  "Events"));
    vh.push_back(PHS("muPt",        {PDC("single", "cutMuVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("muPt_zpt_gt_600", {PDC("single", "cutMuVec(pt)",          {dsDY_ll, dsDY_nunu})}, {1, 2}, "genZPt>800", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("muPt_res",    {PDC("sinlgle", "genMatchMuInAccRes",   {dsDY_ll, dsDY_nunu})}, {1, 1}, "", 100, 0, 1.0,   true,  false, "p_{T}(#mu) (gen - reco)/gen",         "Events"));
    vh.push_back(PHS("muMtw",       {PDC("single", "muonsMtw",              {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 600,   true,  false, "#mu M_{T}(W) [GeV]",      "Events"));
    vh.push_back(PHS("nMuGenMatch", {PDC("single", "genMatchMuInAcc(size)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(#mu) gen Matched",      "Events"));
    vh.push_back(PHS("genZPt",      {PDC("single", "genZPt",                {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})},    {4, 3}, "", 100, 0, 1500,  true,  false, "gen Z p_{T} [GeV]",       "Events"));
    vh.push_back(PHS("genZPt2",     {PDC("single", "genZPt",                {dsDY_ll, dsDY_ll_forEff, dsDY_ll_zAccAcc, dsDY_ll_zAccAcc_nw})}, {3, 2}, "", 100, 0, 1500,  true,  true,  "gen Z p_{T} [GeV]",       "Norm Events"));
    vh.push_back(PHS("genZPt3",     {PDC("single", "genZPt",                {dsDY_ll_forEff, dsDY_ll_zAccAcc_nw})},                  {2, 1}, "", 100, 0, 1500,  true,  false, "gen Z p_{T} [GeV]",       "Events"));
    vh.push_back(PHS("genZPt_gr",   {PDC("single", "genZPt",                {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})},    {3, 2}, "", 200, 0, 3000,  true,  true,  "gen Z p_{T} [GeV]",       "Norm Events"));
    //vh.push_back(PHS("genMuPt",     {PDC("single", "genMu(pt)",             {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "gen #mu p_{T} [GeV]",     "Norm Events"));
    vh.push_back(PHS("genMuEta",    {PDC("single", "genMu(eta)",            {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3,    true,  true,  "gen #mu #eta",            "Norm Events"));
    vh.push_back(PHS("genZEta",     {PDC("single", "genZEta",               {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3,    false, true,  "gen Z #eta",              "Norm Events"));
    vh.push_back(PHS("genZmass",    {PDC("single", "genZmass",              {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1000,  true,  true,  "gen M(Z) [GeV]",          "Norm Events"));

    // Simple Reco vs. Gen plots
    vh.push_back(PHS("bestRecoZPt"  ,         {PDC("single",  {{"bestRecoZPt", dsDY_ll},       {"genZPt", dsDY_ll_gen}})},          {1, 2}, "",              100, 0, 1500,  true,  true,  "Z p_{T} [GeV]",           "Norm Events"));
    vh.push_back(PHS("bestRecoZMass_nosel",   {PDC("single",  {{"bestRecoZM", dsDY_ll_noSel},  {"genZmass", dsDY_ll_gen_noSel}})},  {1, 2}, "",              100, 0, 1000,  true,  true,  "M(Z) [GeV]",              "Norm Events"));
    vh.push_back(PHS("bestRecoZMass",         {PDC("single",  {{"bestRecoZM", dsDY_ll},        {"genZmass", dsDY_ll_gen}})},        {1, 2}, "",              100, 0, 1000,  true,  true,  "M(Z) [GeV]",              "Norm Events"));
    vh.push_back(PHS("bestRecoZMass_dPhiCut", {PDC("single",  {{"bestRecoZM", dsDY_ll},        {"genZmass", dsDY_ll_gen}})},        {1, 2}, "dPhi0_CUT<0.5", 100, 0, 1000,  true,  true,  "M(Z) [GeV]",              "Norm Events"));
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
    Plotter::DataCollection cleanht(     "single", {{"HTZinv", dsDY_nunu},             {"HTZinv", dsDY_ll_zAcc},                {"HTZinv",                dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanmht(    "single", {{"cleanMHt", dsDY_nunu},           {"cleanMHt", dsDY_ll_zAcc},              {"cleanMHt",              dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu},        {"cleanMHtPhi", dsDY_ll_zAcc},           {"cleanMHtPhi",           dsDY_ll_trig_scaled} });
    Plotter::DataCollection cleanMet(    "single", {{"cleanMetPt", dsDY_nunu},         {"cleanMetPt", dsDY_ll_zAcc},            {"cleanMetPt",            dsDY_ll_trig_scaled} });
    Plotter::DataCollection nCleanJet(   "single", {{"cntNJetsPt30Eta24", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_trig_scaled} });
    Plotter::DataCollection nTop(        "single", {{"nTopCandSortedCnt", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc}, {"nTopCandSortedCntZinv", dsDY_ll_trig_scaled} });
    Plotter::DataCollection nBottom(     "single", {{"cntCSVS", dsDY_nunu},            {"cntCSVSZinv", dsDY_ll_zAcc},           {"cntCSVSZinv",           dsDY_ll_trig_scaled} });

    Plotter::DataCollection scaled_ht(   "single", "ht",  {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_mht(  "single", "mht", {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_met(  "single", "met", {{dsDY_nunu},            {dsDY_ll_zWeight_scaled},               {dsDY_ll_zAcc_scaled}               });

    Plotter::DataCollection scaled_cleanht(     "single", {{"HTZinv", dsDY_nunu},                 {"HTZinv", dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_cleanmt2(    "single", {{"best_had_brJet_MT2Zinv", dsDY_nunu}, {"best_had_brJet_MT2Zinv", dsDY_ll_zAcc_scaled}});
    Plotter::DataCollection scaled_cleanmht(    "single", {{"cleanMHt", dsDY_nunu},               {"cleanMHt", dsDY_ll_zAcc_scaled}              });
    Plotter::DataCollection scaled_cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu},            {"cleanMHtPhi", dsDY_ll_zAcc_scaled}           });
    Plotter::DataCollection scaled_cleanMet(    "single", {{"cleanMetPt", dsDY_nunu},             {"cleanMetPt", dsDY_ll_zAcc_scaled}            });
    Plotter::DataCollection scaled_nCleanJet(   "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc_scaled} });
    Plotter::DataCollection scaled_nTop(        "single", {{"nTopCandSortedCntZinv", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc_scaled} });
    Plotter::DataCollection scaled_nBottom(     "single", {{"cntCSVSZinv", dsDY_nunu},            {"cntCSVSZinv", dsDY_ll_zAcc_scaled}           });

    Plotter::DataCollection scaled_nSearchBin(     "single", {{"nSearchBin", dsDY_nunu},             {"nSearchBin", dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nSearchBinNb0(  "single", {{"nSearchBin", dsDY_nunu},             {"nb0Bins",    dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nSearchBinNb0ZZ("single", {{"nSearchBin", dsDY_nunu_SBMb},        {"nb0Bins",    dsDY_nunu_SB0b}       });

    PDC trg_cleanht(  "single","HTZinv",               {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_cleanmht( "single","cleanMHt",             {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_cleanMet( "single","cleanMetPt",           {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nCleanJet("single","cntNJetsPt30Eta24Zinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nTop(     "single","nTopCandSortedCntZinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
    PDC trg_nBottom(  "single","cntCSVSZinv",          {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});

    Plotter::DataCollection scaled_stacked_DYtoll_genht(  "stack", "genHt",    {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

    Plotter::DataCollection scaled_stacked_DYtoll_cleanht(  "stack", "HTZinv",     {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_ht(       "stack", "ht",         {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_cleanmht( "stack", "cleanMHt",   {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
    Plotter::DataCollection scaled_stacked_DYtoll_cleanMet( "stack", "cleanMetPt", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

    Plotter::DataCollection cleanht_znunu(     "single", {{"HTZinv",    dsDY_nunu}});
    Plotter::DataCollection ht_znunu(          "single", {{"ht",         dsDY_nunu}});
    Plotter::DataCollection cleanmht_znunu(    "single", {{"cleanMHt",   dsDY_nunu}});
    Plotter::DataCollection mht_znunu(         "single", {{"mht",        dsDY_nunu}});
    Plotter::DataCollection cleanMet_znunu(    "single", {{"cleanMetPt", dsDY_nunu}});
    Plotter::DataCollection met_znunu(         "single", {{"met",        dsDY_nunu}});

    Plotter::DataCollection scaled_stacked_DYtonunu_cleanht(  "stack", "HTZinv",     {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});
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

    vh.push_back(PHS("trigStudy_baseline_cleanht",   {trg_cleanht  },   {1, 1}, "passBaselineZinv",     50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleanmht",  {trg_cleanmht },   {1, 1}, "passBaselineZinv",     50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleanmet",  {trg_cleanMet },   {1, 1}, "passBaselineZinv",     50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigStudy_baseline_nTop",      {trg_nTop     },   {1, 1}, "passBaselineZinv",     10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigStudy_baseline_cleannJet", {trg_nCleanJet},   {1, 1}, "passBaselineZinv",     20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigStudy_baseline_nBottom",   {trg_nBottom  },   {1, 1}, "passBaselineZinv",     10,  0,       10, true,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigStudy_baselineNoTag_cleanht",   {trg_cleanht  },   {1, 1}, "passBaselineNoTagZinv",     50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleanmht",  {trg_cleanmht },   {1, 1}, "passBaselineNoTagZinv",     50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleanmet",  {trg_cleanMet },   {1, 1}, "passBaselineNoTagZinv",     50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_nTop",      {trg_nTop     },   {1, 1}, "passBaselineNoTagZinv",     10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_cleannJet", {trg_nCleanJet},   {1, 1}, "passBaselineNoTagZinv",     20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("trigStudy_baselineNoTag_nBottom",   {trg_nBottom  },   {1, 1}, "passBaselineNoTagZinv",     10,  0,       10, true,  false,  "N(b)",           "Events"));

    vh.push_back(PHS("trigRatio_singleMu45_cleanht",   {PDC("ratio", "HTZinv",                {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  50,  0, 2000, false,  false,  "H_{T} [GeV]",    "Ratio"));
    vh.push_back(PHS("trigRatio_singleMu45_cleanmht",  {PDC("ratio", "cleanMHt",              {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  50,  0, 1500, false,  false,  "MH_{T} [GeV]",   "Ratio"));
    vh.push_back(PHS("trigRatio_singleMu45_cleanmet",  {PDC("ratio", "cleanMetPt",            {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  50,  0, 1500, false,  false,  "MET [GeV]",      "Ratio"));
    vh.push_back(PHS("trigRatio_singleMu45_nTop",      {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  10,  0,   10, false,  false,  "N(t)",           "Ratio"));
    vh.push_back(PHS("trigRatio_singleMu45_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  20,  0,   20, false,  false,  "N(jet)",         "Ratio"));
    vh.push_back(PHS("trigRatio_singleMu45_nBottom",   {PDC("ratio", "cntCSVSZinv",           {dsDY_ll_singleMu45, dsDY_ll_zAcc})},   {1, 1}, "passBaselineNoTagZinv",  10,  0,   10, false,  false,  "N(b)",           "Ratio"));

    vh.push_back(PHS("genht_DY_mumu_stack",           {scaled_stacked_DYtoll_genht},    {1, 1}, "",                                    100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanht_DY_mumu_stack_er",      {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "",                 150, 0,     3000, true,  false,  "H_{T} [GeV]",    "Events"));

    vh.push_back(PHS("cleanht_DY_mumu_stack",           {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("ht_DY_mumu_stack",                {scaled_stacked_DYtoll_ht,       ht_znunu      },    {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht_DY_mumu_stack",          {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "",                                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet_DY_mumu_stack",          {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "",                                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("baseline_cleanht_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "passBaselineZinv",                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanmht_DY_mumu_stack", {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "passBaselineZinv",                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmet_DY_mumu_stack", {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "passBaselineZinv",                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
                                                                                                       
    vh.push_back(PHS("nb_0_cleanht_DY_mumu_stack",      {scaled_stacked_DYtoll_cleanht,  cleanht_znunu },    {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanmht_DY_mumu_stack",     {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu},    {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmet_DY_mumu_stack",     {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu},    {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("cleanht_DY_nunu_stack",           {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "",                                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht_DY_nunu_stack",          {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "",                                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet_DY_nunu_stack",          {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "",                                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("baseline_cleanht_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "passBaselineZinv",                     100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanmht_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passBaselineZinv",                     100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmet_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passBaselineZinv",                     100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("nb_0_cleanht_DY_nunu_stack",      {scaled_stacked_DYtonunu_cleanht},   {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanmht_DY_nunu_stack",     {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmet_DY_nunu_stack",     {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));

    vh.push_back(PHS("nSearchBin",               {scaled_nSearchBin},      {2, 1}, "passBaselineZinv",                                             45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_log",           {scaled_nSearchBin},      {2, 1}, "passBaselineZinv",                                             45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBinnb0",            {scaled_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBinnb0_log",        {scaled_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBinnb0ZZ",          {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBinnb0ZZ_log",      {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_pull",          {scaled_nSearchBin},      {2, 1}, "passBaselineZinv",                                             45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("nSearchBin_pull_log",      {scaled_nSearchBin},      {2, 1}, "passBaselineZinv",                                             45,  0,     45,   true,  false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("nSearchBinnb0_pull",       {scaled_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("nSearchBinnb0_pull_log",   {scaled_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   true,  false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("nSearchBinnb0ZZ_pull",     {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("nSearchBinnb0ZZ_pull_log", {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   true,  false,  "Search Bin",     "Events", false));

    vh.push_back(PHS("mT2Zinv",                {scaled_cleanmt2},      {2, 1}, "",                                                             100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_mT2Zinv",       {scaled_cleanmt2},      {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("baselineNoTag_mT2Zinv",  {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nt_0_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nt_1_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nt_2_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nt_3_mT2Zinv",           {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_0_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_0_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_0_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_0_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_1_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_1_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_1_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_1_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_2_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_2_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_2_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_2_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_3_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_3_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_3_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_3_mT2Zinv",      {scaled_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("cleanht",                {scaled_cleanht},       {2, 1}, "",                                                             100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_cleanht",       {scaled_cleanht},       {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baselineNoTag_cleanht",  {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("baseline_ht",            {scaled_ht},            {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_0_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_1_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_2_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nt_3_cleanht",           {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanht",      {scaled_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmht",               {scaled_cleanmht},      {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_cleanmht",      {scaled_cleanmht},      {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baselineNoTag_cleanmht", {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("baseline_mht",           {scaled_mht},           {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_0_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_1_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_2_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nt_3_cleanmht",          {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanmht",     {scaled_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("cleanmet",               {scaled_cleanMet},      {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("cleanmet_gr",            {scaled_cleanMet},      {2, 1}, "",                                                             200, 0,     3000, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("baseline_cleanmet",      {scaled_cleanMet},      {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("baselineNoTag_cleanmet", {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_0_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_1_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_2_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nt_3_cleanmet",          {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_0_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_0_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_0_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_0_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_1_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_1_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_1_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_1_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_2_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_2_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_2_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_2_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_0_nt_3_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_1_nt_3_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_2_nt_3_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nb_3_nt_3_cleanmet",     {scaled_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("nTop",                   {scaled_nTop},          {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("baseline_nTop",          {scaled_nTop},          {2, 1}, "passBaselineZinv",                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("baselineNoTag_nTop",     {scaled_nTop},          {2, 1}, "passBaselineNoTagZinv",                                        10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_0_nTop",              {scaled_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_1_nTop",              {scaled_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_2_nTop",              {scaled_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("nb_3_nTop",              {scaled_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("cleannJet",              {scaled_nCleanJet},     {2, 1}, "",                                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("baseline_cleannJet",     {scaled_nCleanJet},     {2, 1}, "passBaselineZinv",                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("baselineNoTag_cleannJet",{scaled_nCleanJet},     {2, 1}, "passBaselineNoTagZinv",                                        20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("nBottom",                {scaled_nBottom},       {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("baseline_nBottom",       {scaled_nBottom},       {2, 1}, "passBaselineZinv",                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("baselineNoTAg_nBottom",  {scaled_nBottom},       {2, 1}, "passBaselineNoTagZinv",                                        10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("cleanmhtphi",            {scaled_cleanmhtphi},   {2, 1}, "",                                                             100, -3.14, 3.14, false, false,  "#phi(MH_{T})",   "Events"));
    vh.push_back(PHS("baseline_cleanmhtphi",   {scaled_cleanmhtphi},   {2, 1}, "passBaselineZinv",                                             100, -3.14, 3.14, false, false,  "#phi(MH_{T})",   "Events"));

    vh.push_back(PHS("fake1b_nTop",               {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_1b},         {"nTopCandSortedCntZinv1b", dsDY_nunu_0b_1fake}})},         {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));
    vh.push_back(PHS("fake1b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_1b_blNoTag}, {"nTopCandSortedCntZinv1b", dsDY_nunu_0b_1fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));

    vh.push_back(PHS("fake2b_nTop",               {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_2b},         {"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})},         {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));
    vh.push_back(PHS("fake2b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_2b_blNoTag}, {"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));

    vh.push_back(PHS("fake3b_nTop",               {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_3b},         {"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})},         {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));
    vh.push_back(PHS("fake3b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_3b_blNoTag}, {"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, "N(t)", "Norm Events"));

    vh.push_back(PHS("fake2bvs0b_nTop",               {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})}, {2, 1}, "",                      10, 0, 10, true, true, "N(t)", "Norm Events"));
    vh.push_back(PHS("fake2bvs0b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, "N(t)", "Norm Events"));

    vh.push_back(PHS("fake3bvs0b_nTop",               {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})}, {2, 1}, "",                      10, 0, 10, true, true, "N(t)", "Norm Events"));
    vh.push_back(PHS("fake3bvs0b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, "N(t)", "Norm Events"));

    vh.push_back(PHS("fake1b_MT2",               {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_1b},         {"best_had_brJet_MT2Zinv1b", dsDY_nunu_0b_1fake}})},         {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));
    vh.push_back(PHS("fake1b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_1b_blNoTag}, {"best_had_brJet_MT2Zinv1b", dsDY_nunu_0b_1fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));

    vh.push_back(PHS("fake2b_MT2",               {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_2b},         {"best_had_brJet_MT2Zinv2b", dsDY_nunu_0b_2fake}})},         {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));
    vh.push_back(PHS("fake2b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_2b_blNoTag}, {"best_had_brJet_MT2Zinv2b", dsDY_nunu_0b_2fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));

    vh.push_back(PHS("fake3b_MT2",               {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_3b},         {"best_had_brJet_MT2Zinv3b", dsDY_nunu_0b_3fake}})},         {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));
    vh.push_back(PHS("fake3b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_3b_blNoTag}, {"best_had_brJet_MT2Zinv3b", dsDY_nunu_0b_3fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));

    vh.push_back(PHS("fake2bvs0b_MT2",               {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv2b",dsDY_nunu_0b_2fake}})}, {2, 1}, "",                      50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));
    vh.push_back(PHS("fake2bvs0b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv2b",dsDY_nunu_0b_2fake}})}, {2, 1}, "passBaselineNoTagZinv", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));

    vh.push_back(PHS("fake3bvs0b_MT2",               {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv3b",dsDY_nunu_0b_3fake}})}, {2, 1}, "",                      50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));
    vh.push_back(PHS("fake3bvs0b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv3b",dsDY_nunu_0b_3fake}})}, {2, 1}, "passBaselineNoTagZinv", 50, 0, 2000, true, true, "M_{T2} [GeV]", "Norm Events"));

    vh.push_back(PHS("fakedCSVByOrder",                 {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "",                      100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVByOrder0b",               {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "cntCSVSZinv=0",         100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVByOrder_baseline",        {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineZinv",      100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVByOrder_baselineNoTag",   {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineNoTagZinv", 100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVByOrder0b_baselineNoTag", {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0", 100, 0, 1, false, true, "CSV value", "Norm Events"));

    vh.push_back(PHS("fakedCSVAll",                 {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "cntCSVSZinv=0",                       100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVAll0b",               {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "",                                    100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVAll_baseline",        {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineZinv",                    100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVAll_baselineNoTag",   {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineNoTagZinv",               100, 0, 1, false, true, "CSV value", "Norm Events"));
    vh.push_back(PHS("fakedCSVAll0b_baselineNoTag", {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0", 100, 0, 1, false, true, "CSV value", "Norm Events"));

    vh.push_back(PHS("bybTag_nJets",               {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "",                      30, 0, 30, true, true, "N jets", "Norm Events"));
    vh.push_back(PHS("bybTag_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, "N jets", "Norm Events"));

    vh.push_back(PHS("bybTag_nJets_1f",               {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_1fake, dsDY_nunu_1b})}, {2, 1}, "",                      30, 0, 30, true, true, "N jets", "Norm Events"));
    vh.push_back(PHS("bybTag_nJets_1f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_1fake, dsDY_nunu_1b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, "N jets", "Norm Events"));

    vh.push_back(PHS("bybTag_nJets_2f",               {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_2fake, dsDY_nunu_2b})}, {2, 1}, "",                      30, 0, 30, true, true, "N jets", "Norm Events"));
    vh.push_back(PHS("bybTag_nJets_2f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_2fake, dsDY_nunu_2b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, "N jets", "Norm Events"));

    vh.push_back(PHS("bybTag_nJets_3f",               {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_3fake, dsDY_nunu_3b})}, {2, 1}, "",                      30, 0, 30, true, true, "N jets", "Norm Events"));
    vh.push_back(PHS("bybTag_nJets_3f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_3fake, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, "N jets", "Norm Events"));

    vh.push_back(PHS("bybTag_cleanMET",               {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "",                      100, 0, 1500, true, true, "MET [GeV]", "Norm Events"));
    vh.push_back(PHS("bybTag_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "MET [GeV]", "Norm Events"));

    vh.push_back(PHS("bybTag_cleanjPt",               {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "",                      100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("bybTag_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));

    vh.push_back(PHS("bybTag_MT2",               {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "",                      100, 0, 2000, true, true, "MT2 [GeV]", "Norm Events"));
    vh.push_back(PHS("bybTag_MT2_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, "MT2 [GeV]", "Norm Events"));

    vh.push_back(PHS("bybTag_Ntops",               {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "",                      10, 0, 10, true, true, "N tops", "Norm Events"));
    vh.push_back(PHS("bybTag_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, "N tops", "Norm Events"));


    vh.push_back(PHS("byCSV_nJets",               {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "",                      30, 0, 30, true, true, "N jets", "Norm Events"));
    vh.push_back(PHS("byCSV_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, "N jets", "Norm Events"));

    vh.push_back(PHS("byCSV_cleanMET",               {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "",                      100, 0, 1500, true, true, "MET [GeV]", "Norm Events"));
    vh.push_back(PHS("byCSV_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "MET [GeV]", "Norm Events"));

    vh.push_back(PHS("byCSV_cleanjPt",               {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "",                      100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("byCSV_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));

    vh.push_back(PHS("byCSV_MT2",               {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "",                      100, 0, 2000, true, true, "MT2 [GeV]", "Norm Events"));
    vh.push_back(PHS("byCSV_MT2_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, "MT2 [GeV]", "Norm Events"));

    vh.push_back(PHS("byCSV_Ntops",               {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "",                      10, 0, 10, true, true, "N tops", "Norm Events"));
    vh.push_back(PHS("byCSV_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, "N tops", "Norm Events"));


    vh.push_back(PHS("byCSV1b_nJets_baselineNoTag",    {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv",  30, 0,   30, true, true, "N jets",        "Norm Events"));
    vh.push_back(PHS("byCSV1b_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt",             {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "MET [GeV]",     "Norm Events"));
    vh.push_back(PHS("byCSV1b_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("byCSV1b_MT2_baselineNoTag",      {PDC("single", "cleanMetPt",             {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, "MT2 [GeV]",     "Norm Events"));
    vh.push_back(PHS("byCSV1b_Ntops_baselineNoTag",    {PDC("single", "nTopCandSortedCntZinv",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv",  10, 0,   10, true, true, "N tops",        "Norm Events"));

    vh.push_back(PHS("byCSVvh1b_nJets_baselineNoTag",    {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv",  30, 0,   30, true, true, "N jets",        "Norm Events"));
    vh.push_back(PHS("byCSVvh1b_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt",             {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "MET [GeV]",     "Norm Events"));
    vh.push_back(PHS("byCSVvh1b_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("byCSVvh1b_MT2_baselineNoTag",      {PDC("single", "cleanMetPt",             {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, "MT2 [GeV]",     "Norm Events"));
    vh.push_back(PHS("byCSVvh1b_Ntops_baselineNoTag",    {PDC("single", "nTopCandSortedCntZinv",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv",  10, 0,   10, true, true, "N tops",        "Norm Events"));


    vh.push_back(PHS("muPt_met_lt_150", {PDC("single", "cutMuVec(pt)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "cleanMetPt<150", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("cleanht_met_lt_150",   {scaled_cleanht}, {1, 1}, "cleanMetPt<150", 50, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("cleanmet_met_lt_150", {scaled_cleanMet}, {1, 1}, "cleanMetPt<150", 50, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("cleanmht_met_lt_150", {scaled_cleanmht}, {1, 1}, "cleanMetPt<150", 50, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));


    // Plots the information needed to calculate acceptance corrections, plot itself is meaningless
    vh.push_back(PHS("accInfo",  {PDC("single", {{"genMuInAcc(size)", dsDY_ll_forEff}, {"genMu(size)", dsDY_ll_forEff}})},    {1, 2}, "",  20,  0,       20, true, false,  "gen N(#mu)",     "Events"));

    vh.push_back(PHS("genMuPt",  {PDC("single", {{"genMuInAcc(pt)", dsDY_ll_forEff},            {"genMu(pt)", dsDY_ll_forEff}})},             {1, 1}, "",  150,  0,  1500, true,  false,  "gen #mu p_{T} [GeV]", "Events"));
    vh.push_back(PHS("accHt",    {PDC("ratio",  {{"HTZinv", dsDY_ll_forEff_den},                {"HTZinv", dsDY_ll_forAcc}})},                {1, 1}, "",  300,  0,  3000, false, false,  "H_{T} [GeV]",         "Acc"));
    vh.push_back(PHS("accMt",    {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_den}, {"nTopCandSortedCntZinv", dsDY_ll_forAcc}})}, {1, 1}, "",  10,   0,  10,   false, false,  "N(t)",                "Acc"));
    vh.push_back(PHS("accNb",    {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_den},           {"cntCSVSZinv", dsDY_ll_forAcc}})},           {1, 1}, "",  10,   0,  10,   false, false,  "N(b)",                "Acc"));
    vh.push_back(PHS("accNjet",  {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forAcc}})}, {1, 1}, "",  20,   0,  20,   false, false,  "N(j)",                "Acc"));

    // Plots for calculation of muon efficiency 
    vh.push_back(PHS("effMuPt",      {PDC("ratio",  {{"genMatchMuInAcc(pt)", dsDY_ll_forEff},       {"genMuInAcc(pt)", dsDY_ll_forEff}})},             {1, 1}, "", 200, 0,     2000, false, false,  "#mu p_{T} [GeV]",  "Efficiency"));
    vh.push_back(PHS("effMuEta",     {PDC("ratio",  {{"genMatchMuInAcc(eta)", dsDY_ll_forEff},      {"genMuInAcc(eta)", dsDY_ll_forEff}})},            {1, 1}, "", 100, -3.0,  3.0,  false, false,  "#mu #eta",         "Efficiency"));
    vh.push_back(PHS("effMucleanHt", {PDC("ratio",  {{"HTZinv", dsDY_ll_forEff_num},                {"HTZinv", dsDY_ll_forEff_den}})},                 {1, 1}, "", 300, 0,     3000, false, false,  "H_{T} [GeV]",      "Efficiency"));
    vh.push_back(PHS("effMuNb",      {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_num},           {"cntCSVSZinv", dsDY_ll_forEff_den}})},            {1, 1}, "", 10, 0,      10.0, false, false,  "N(b)",             "Efficiency"));
    vh.push_back(PHS("effMuNt",      {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_num}, {"nTopCandSortedCntZinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 10, 0,      10.0, false, false,  "N(t)",             "Efficiency"));
    vh.push_back(PHS("effMuNjet",    {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_num}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 20, 0,      20.0, false, false,  "N(j)",             "Efficiency"));

    // Plots of weight and efficiencies used as sanity checks
    vh.push_back(PHS("zEffWgt", {PDC("single", "zEffWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Eff Weight",         ""));
    vh.push_back(PHS("zEff",    {PDC("single", "zEff",     {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0,   false, false,  "Z Eff Weight",         ""));
    vh.push_back(PHS("zAccWgt", {PDC("single", "zAccWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Acc Weight",         ""));
    vh.push_back(PHS("zAcc",    {PDC("single", "zAcc",     {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0,   false, false,  "Z Acc Weight",         ""));

    PDC dcDY_cutFlow("single", "cleanMetPt", {dsDY_ll_passMuonVeto, dsDY_ll_passEleVeto, dsDY_ll_passIsoTrkVeto, dsDY_ll_passnJets, dsDY_ll_passdPhis, dsDY_ll_passMET, dsDY_ll_passBJets, dsDY_ll_passTagger, dsDY_ll_passBaseline, dsDY_ll_passBaselineNoTag});
    vh.push_back(PHS("gcf_met",    {dcDY_cutFlow},  {1, 1}, "", 100, 0.0,  3000.0,   true, false,   "MET [GeV]",         "Events"));

    PDC dcDY_cutFlowZinv("single", "cleanMetPt", {dsDY_ll_passMuonVetoZinv, dsDY_ll_passEleVetoZinv, dsDY_ll_passIsoTrkVetoZinv, dsDY_ll_passnJetsZinv, dsDY_ll_passdPhisZinv, dsDY_ll_passMETZinv, dsDY_ll_passBJetsZinv, dsDY_ll_passTaggerZinv, dsDY_ll_passBaselineZinv, dsDY_ll_passBaselineNoTagZinv});
    vh.push_back(PHS("gcfZinv_met",    {dcDY_cutFlowZinv},  {1, 1}, "", 100, 0.0,  3000.0,   true, false,   "MET [GeV]",         "Events"));

    // Test plots, these are simply meaningless demonstrations of features
    Plotter::DataCollection dcDY_test(   "single", {{"mht", dsDY_ll}, {"jetsLVec(pt)", dsDY_nunu}, {"genZmass", dsDY_test}, {"genZPt", dsDY_test}});
    Plotter::DataCollection dcDY_ratio(  "ratio",  {{"cleanJetVec(pt)", dsDY_ll}, {"jetsLVec(pt)", dsDY_ll}});
    Plotter::DataCollection dcDY_stack(  "stack",  "mht", {dsDY_ll, dsDY_nunu});
    Plotter::DataCollection dcDY_dataTest(  "data",  "met", {dsDY_nunu});

    vh.push_back(PHS("test",  {dcDY_stack, dcDY_test, dcDY_dataTest},  {3, 2}, "genZmass>80", 100, 0, 1000,   true, true,  "???",         "Norm Events"));
    vh.push_back(PHS("test2", {dcDY_ratio},  {1, 1}, "", 100, 0, 500, true, false, "jpt???", "Ratio"));

    // ------------------------
    // - Data/MC plots
    // ------------------------

    Plotter::DatasetSummary dsData_2015B_nosel("Data Run 2015B", fileMap["SingleMuon_2015B"], "", "");
    Plotter::DatasetSummary dsData_2015C_nosel("Data Run 2015C", fileMap["SingleMuon_2015C"], "", "");
    Plotter::DatasetSummary dsData_2015D_nosel("Data Run 2015D", fileMap["SingleMuon_2015D"], "", "");
    Plotter::DatasetSummary dsDY_nosel("DY", fileMap["DYJetsToLL"], "", "");
    Plotter::DatasetSummary dstt2l_nosel("t#bar{t} dilepton", fileMap["TTbarDiLep"], "", "");

    Plotter::DatasetSummary dsData_2015B_blnotag("Data Run 2015B", fileMap["SingleMuon_2015B"], "passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsData_2015C_blnotag("Data Run 2015C", fileMap["SingleMuon_2015C"], "passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsData_2015D_blnotag("Data Run 2015D", fileMap["SingleMuon_2015D"], "passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dsDY_blnotag("DY", fileMap["DYJetsToLL"], "passBaselineNoTagZinv", "");
    Plotter::DatasetSummary dstt2l_blnotag("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagZinv", "");

    Plotter::DatasetSummary dsData_2015B_bl("Data Run 2015B", fileMap["SingleMuon_2015B"], "passBaselineZinv", "");
    Plotter::DatasetSummary dsData_2015C_bl("Data Run 2015C", fileMap["SingleMuon_2015C"], "passBaselineZinv", "");
    Plotter::DatasetSummary dsData_2015D_bl("Data Run 2015D", fileMap["SingleMuon_2015D"], "passBaselineZinv", "");
    Plotter::DatasetSummary dsDY_bl("DY", fileMap["DYJetsToLL"], "passBaselineZinv", "");
    Plotter::DatasetSummary dstt2l_bl("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineZinv", "");


    Plotter::DataCollection dcData_2015B_met_nosel("data",  "cleanMetPt", {dsData_2015B_nosel});
    Plotter::DataCollection dcData_2015CD_met_nosel("data",  "cleanMetPt", {dsData_2015C_nosel, dsData_2015D_nosel});
    Plotter::DataCollection dcDY_met_nosel("stack",  "cleanMetPt", {dsDY_nosel, dstt2l_nosel});

    Plotter::DataCollection dcData_2015B_met_blnotag("data",  "cleanMetPt", {dsData_2015B_blnotag});
    Plotter::DataCollection dcData_2015CD_met_blnotag("data",  "cleanMetPt", {dsData_2015C_blnotag, dsData_2015D_blnotag});
    Plotter::DataCollection dcDY_met_blnotag("stack",  "cleanMetPt", {dsDY_blnotag, dstt2l_blnotag});

    Plotter::DataCollection dcData_2015B_met_bl("data",  "cleanMetPt", {dsData_2015B_bl});
    Plotter::DataCollection dcData_2015CD_met_bl("data",  "cleanMetPt", {dsData_2015C_bl, dsData_2015D_bl});
    Plotter::DataCollection dcDY_met_bl("stack",  "cleanMetPt", {dsDY_bl, dstt2l_bl});


    Plotter::DataCollection dcData_2015B_nt_nosel("data",  "nTopCandSortedCntZinv", {dsData_2015B_nosel});
    Plotter::DataCollection dcData_2015CD_nt_nosel("data",  "nTopCandSortedCntZinv", {dsData_2015C_nosel, dsData_2015D_nosel});
    Plotter::DataCollection dcDY_nt_nosel("stack",  "nTopCandSortedCntZinv", {dsDY_nosel, dstt2l_nosel});

    Plotter::DataCollection dcData_2015B_nt_blnotag("data",  "nTopCandSortedCntZinv", {dsData_2015B_blnotag});
    Plotter::DataCollection dcData_2015CD_nt_blnotag("data",  "nTopCandSortedCntZinv", {dsData_2015C_blnotag, dsData_2015D_blnotag});
    Plotter::DataCollection dcDY_nt_blnotag("stack",  "nTopCandSortedCntZinv", {dsDY_blnotag, dstt2l_blnotag});

    Plotter::DataCollection dcData_2015B_nt_bl("data",  "nTopCandSortedCntZinv", {dsData_2015B_bl});
    Plotter::DataCollection dcData_2015CD_nt_bl("data",  "nTopCandSortedCntZinv", {dsData_2015C_bl, dsData_2015D_bl});
    Plotter::DataCollection dcDY_nt_bl("stack",  "nTopCandSortedCntZinv", {dsDY_bl, dstt2l_bl});


    Plotter::DataCollection dcData_2015B_nb_nosel("data",  "cntCSVSZinv", {dsData_2015B_nosel});
    Plotter::DataCollection dcData_2015CD_nb_nosel("data",  "cntCSVSZinv", {dsData_2015C_nosel, dsData_2015D_nosel});
    Plotter::DataCollection dcDY_nb_nosel("stack",  "cntCSVSZinv", {dsDY_nosel, dstt2l_nosel});

    Plotter::DataCollection dcData_2015B_nb_blnotag("data",  "cntCSVSZinv", {dsData_2015B_blnotag});
    Plotter::DataCollection dcData_2015CD_nb_blnotag("data",  "cntCSVSZinv", {dsData_2015C_blnotag, dsData_2015D_blnotag});
    Plotter::DataCollection dcDY_nb_blnotag("stack",  "cntCSVSZinv", {dsDY_blnotag, dstt2l_blnotag});

    Plotter::DataCollection dcData_2015B_nb_bl("data",  "cntCSVSZinv", {dsData_2015B_bl});
    Plotter::DataCollection dcData_2015CD_nb_bl("data",  "cntCSVSZinv", {dsData_2015C_bl, dsData_2015D_bl});
    Plotter::DataCollection dcDY_nb_bl("stack",  "cntCSVSZinv", {dsDY_bl, dstt2l_bl});


    Plotter::DataCollection dcData_2015B_nj_nosel("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015B_nosel});
    Plotter::DataCollection dcData_2015CD_nj_nosel("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_nosel, dsData_2015D_nosel});
    Plotter::DataCollection dcDY_nj_nosel("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_nosel, dstt2l_nosel});

    Plotter::DataCollection dcData_2015B_nj_blnotag("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015B_blnotag});
    Plotter::DataCollection dcData_2015CD_nj_blnotag("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_blnotag, dsData_2015D_blnotag});
    Plotter::DataCollection dcDY_nj_blnotag("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_blnotag, dstt2l_blnotag});

    Plotter::DataCollection dcData_2015B_nj_bl("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015B_bl});
    Plotter::DataCollection dcData_2015CD_nj_bl("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_bl, dsData_2015D_bl});
    Plotter::DataCollection dcDY_nj_bl("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_bl, dstt2l_bl});


    Plotter::DataCollection dcData_2015B_ht_nosel("data",  "HTZinv", {dsData_2015B_nosel});
    Plotter::DataCollection dcData_2015CD_ht_nosel("data",  "HTZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_ht_nosel("stack",  "HTZinv", {dsDY_nosel, dstt2l_nosel});

    Plotter::DataCollection dcData_2015B_ht_blnotag("data",  "HTZinv", {dsData_2015B_blnotag});
    Plotter::DataCollection dcData_2015CD_ht_blnotag("data",  "HTZinv", {dsData_2015C_blnotag, dsData_2015D_blnotag});
    Plotter::DataCollection dcDY_ht_blnotag("stack",  "HTZinv", {dsDY_blnotag, dstt2l_blnotag});

    Plotter::DataCollection dcData_2015B_ht_bl("data",  "HTZinv", {dsData_2015B_bl});
    Plotter::DataCollection dcData_2015CD_ht_bl("data",  "HTZinv", {dsData_2015C_bl, dsData_2015D_bl});
    Plotter::DataCollection dcDY_ht_bl("stack",  "HTZinv", {dsDY_bl, dstt2l_bl});


    vh.push_back(PHS("DataMC_2015B_met_nosel", {dcData_2015B_met_nosel, dcDY_met_nosel},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",    ""));
    vh.push_back(PHS("DataMC_2015B_ht_nosel",  {dcData_2015B_ht_nosel, dcDY_ht_nosel},    {1, 2}, "", 150, 0, 1500,   true, false,  "ht",     ""));
    vh.push_back(PHS("DataMC_2015B_nt_nosel",  {dcData_2015B_nt_nosel, dcDY_nt_nosel},    {1, 2}, "", 5, 0, 5,        true, false,  "ntop",   ""));
    vh.push_back(PHS("DataMC_2015B_nb_nosel",  {dcData_2015B_nb_nosel, dcDY_nb_nosel},    {1, 2}, "", 10, 0, 10,      true, false,  "nb",     ""));
    vh.push_back(PHS("DataMC_2015B_nj_nosel",  {dcData_2015B_nj_nosel, dcDY_nj_nosel},    {1, 2}, "", 20, 0, 20,      true, false,  "nj",     ""));
    vh.push_back(PHS("DataMC_2015B_mht_nosel",  {PDC("data", "cleanMHt", {dsData_2015B_nosel}), PDC("stack", "cleanMHt", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015B_jpt_nosel",  {PDC("data", "cleanJetVec(pt)", {dsData_2015B_nosel}), PDC("stack", "cleanJetVec(pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j1pt_nosel",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015B_nosel}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j2pt_nosel",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015B_nosel}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j3pt_nosel",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015B_nosel}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mt2_nosel",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015B_nosel}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015B_mupt_nosel",  {PDC("data", "cutMuVec(pt)", {dsData_2015B_nosel}), PDC("stack", "cutMuVec(pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu1pt_nosel",  {PDC("data", "cutMuPt1", {dsData_2015B_nosel}), PDC("stack", "cutMuPt1", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu2pt_nosel",  {PDC("data", "cutMuPt2", {dsData_2015B_nosel}), PDC("stack", "cutMuPt2", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    
    vh.push_back(PHS("DataMC_2015B_met_baselineNoTag",  {dcData_2015B_met_blnotag, dcDY_met_blnotag}, {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015B_ht_baselineNoTag",  {dcData_2015B_ht_blnotag,   dcDY_ht_blnotag},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015B_nt_baselineNoTag",  {dcData_2015B_nt_blnotag,   dcDY_nt_blnotag},  {1, 2}, "", 5, 0, 5,        true, false,  "ntop",         ""));
    vh.push_back(PHS("DataMC_2015B_nb_baselineNoTag",  {dcData_2015B_nb_blnotag,   dcDY_nb_blnotag},  {1, 2}, "", 10, 0, 10,      true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015B_nj_baselineNoTag",  {dcData_2015B_nj_blnotag,   dcDY_nj_blnotag},  {1, 2}, "", 20, 0, 20,      true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015B_mht_baselineNoTag",  {PDC("data", "cleanMHt", {dsData_2015B_blnotag}), PDC("stack", "cleanMHt", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015B_jpt_baselineNoTag",  {PDC("data", "cleanJetVec(pt)", {dsData_2015B_blnotag}), PDC("stack", "cleanJetVec(pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j1pt_baselineNoTag",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015B_blnotag}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j2pt_baselineNoTag",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015B_blnotag}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j3pt_baselineNoTag",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015B_blnotag}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mt2_baselineNoTag",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015B_blnotag}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015B_mupt_baselineNoTag",  {PDC("data", "cutMuVec(pt)", {dsData_2015B_blnotag}), PDC("stack", "cutMuVec(pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu1pt_baselineNoTag",  {PDC("data", "cutMuPt1", {dsData_2015B_blnotag}), PDC("stack", "cutMuPt1", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu2pt_baselineNoTag",  {PDC("data", "cutMuPt2", {dsData_2015B_blnotag}), PDC("stack", "cutMuPt2", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    
    vh.push_back(PHS("DataMC_2015B_met_baseline",  {dcData_2015B_met_bl, dcDY_met_bl}, {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015B_ht_baseline",  {dcData_2015B_ht_bl,   dcDY_ht_bl},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015B_nt_baseline",  {dcData_2015B_nt_bl,   dcDY_nt_bl},  {1, 2}, "", 5, 0, 5,        true, false,  "ntop",         ""));
    vh.push_back(PHS("DataMC_2015B_nb_baseline",  {dcData_2015B_nb_bl,   dcDY_nb_bl},  {1, 2}, "", 10, 0, 10,      true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015B_nj_baseline",  {dcData_2015B_nj_bl,   dcDY_nj_bl},  {1, 2}, "", 20, 0, 20,      true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015B_mht_baseline",  {PDC("data", "cleanMHt", {dsData_2015B_bl}), PDC("stack", "cleanMHt", {dsDY_bl, dstt2l_bl})},  
    		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015B_jpt_baseline",  {PDC("data", "cleanJetVec(pt)", {dsData_2015B_bl}), PDC("stack", "cleanJetVec(pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j1pt_baseline",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015B_bl}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j2pt_baseline",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015B_bl}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_j3pt_baseline",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015B_bl}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mt2_baseline",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015B_bl}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015B_mupt_baseline",  {PDC("data", "cutMuVec(pt)", {dsData_2015B_bl}), PDC("stack", "cutMuVec(pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu1pt_baseline",  {PDC("data", "cutMuPt1", {dsData_2015B_bl}), PDC("stack", "cutMuPt1", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015B_mu2pt_baseline",  {PDC("data", "cutMuPt2", {dsData_2015B_bl}), PDC("stack", "cutMuPt2", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));


    vh.push_back(PHS("DataMC_2015CD_met_nosel", {dcData_2015CD_met_nosel, dcDY_met_nosel},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015CD_ht_nosel",  {dcData_2015CD_ht_nosel, dcDY_ht_nosel},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_nosel",  {dcData_2015CD_nt_nosel, dcDY_nt_nosel},  {1, 2}, "", 5, 0, 5,   true, false,  "ntop",         ""));
    vh.push_back(PHS("DataMC_2015CD_nb_nosel",  {dcData_2015CD_nb_nosel, dcDY_nb_nosel},  {1, 2}, "", 10, 0, 10,   true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_nosel",  {dcData_2015CD_nj_nosel, dcDY_nj_nosel},  {1, 2}, "", 20, 0, 20,   true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_nosel",  {PDC("data", "cleanMHt", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cleanMHt", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_nosel",  {PDC("data", "cleanJetVec(pt)", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cleanJetVec(pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_nosel",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_nosel",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_nosel",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_nosel",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_nosel",  {PDC("data", "cutMuVec(pt)", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cutMuVec(pt)", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_nosel",  {PDC("data", "cutMuPt1", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cutMuPt1", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_nosel",  {PDC("data", "cutMuPt2", {dsData_2015C_nosel, dsData_2015D_nosel}), PDC("stack", "cutMuPt2", {dsDY_nosel, dstt2l_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    
    vh.push_back(PHS("DataMC_2015CD_met_baselineNoTag", {dcData_2015CD_met_blnotag, dcDY_met_blnotag},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baselineNoTag",  {dcData_2015CD_ht_blnotag, dcDY_ht_blnotag},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baselineNoTag",  {dcData_2015CD_nt_blnotag, dcDY_nt_blnotag},  {1, 2}, "", 5, 0, 5,   true, false,  "ntop",         ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baselineNoTag",  {dcData_2015CD_nb_blnotag, dcDY_nb_blnotag},  {1, 2}, "", 10, 0, 10,   true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baselineNoTag",  {dcData_2015CD_nj_blnotag, dcDY_nj_blnotag},  {1, 2}, "", 20, 0, 20,   true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baselineNoTag",  {PDC("data", "cleanMHt", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cleanMHt", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baselineNoTag",  {PDC("data", "cleanJetVec(pt)", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cleanJetVec(pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baselineNoTag",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baselineNoTag",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baselineNoTag",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baselineNoTag",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baselineNoTag",  {PDC("data", "cutMuVec(pt)", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cutMuVec(pt)", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baselineNoTag",  {PDC("data", "cutMuPt1", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cutMuPt1", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baselineNoTag",  {PDC("data", "cutMuPt2", {dsData_2015C_blnotag, dsData_2015D_blnotag}), PDC("stack", "cutMuPt2", {dsDY_blnotag, dstt2l_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));

    vh.push_back(PHS("DataMC_2015CD_met_baseline", {dcData_2015CD_met_bl, dcDY_met_bl},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baseline",  {dcData_2015CD_ht_bl, dcDY_ht_bl},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baseline",  {dcData_2015CD_nt_bl, dcDY_nt_bl},  {1, 2}, "", 5, 0, 5,   true, false,  "ntop",         ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baseline",  {dcData_2015CD_nb_bl, dcDY_nb_bl},  {1, 2}, "", 10, 0, 10,   true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baseline",  {dcData_2015CD_nj_bl, dcDY_nj_bl},  {1, 2}, "", 20, 0, 20,   true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baseline",  {PDC("data", "cleanMHt", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cleanMHt", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baseline",  {PDC("data", "cleanJetVec(pt)", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cleanJetVec(pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baseline",  {PDC("data", "cleanJetVec[0](pt)", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baseline",  {PDC("data", "cleanJetVec[1](pt)", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baseline",  {PDC("data", "cleanJetVec[2](pt)", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baseline",  {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baseline",  {PDC("data", "cutMuVec(pt)", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cutMuVec(pt)", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baseline",  {PDC("data", "cutMuPt1", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cutMuPt1", {dsDY_bl, dstt2l_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baseline",  {PDC("data", "cutMuPt2", {dsData_2015C_bl, dsData_2015D_bl}), PDC("stack", "cutMuPt2", {dsDY_bl, dstt2l_bl})},  
    		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setPlotDir(plotDir);
    if(doSave)  plotter.saveHists();
    if(doPlots) plotter.plot();
}
