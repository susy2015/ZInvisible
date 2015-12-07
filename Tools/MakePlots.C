#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
    using namespace std;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"plot",             no_argument, 0, 'p'},
        {"savehist",         no_argument, 0, 's'},
        {"savetuple",        no_argument, 0, 't'},
        {"fromFile",         no_argument, 0, 'f'},
        {"condor",           no_argument, 0, 'c'},
        {"histFile",   required_argument, 0, 'H'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
	{"luminosity", required_argument, 0, 'L'}
    };

    bool doPlots = true, doSave = true, doTuple = true, fromTuple = true, runOnCondor = false;
    string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi = AnaSamples::luminosity;

    while((opt = getopt_long(argc, argv, "pstfcH:D:N:M:E:P:L:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'p':
            if(doPlots) doSave  = doTuple = false;
            else        doPlots = true;
            break;

        case 's':
            if(doSave) doPlots = doTuple = false;
            else       doSave  = true;
            break;

        case 't':
            if(doTuple) doPlots = doSave = false;
            else        doTuple  = true;
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
        sampleloc = "condor";
    }

    AnaSamples::SampleSet        ss(sampleloc, lumi);
    AnaSamples::SampleCollection sc(ss);

    const double zAcc = 1.0;
//    const double zAcc = 0.5954;
//    const double zAcc = 0.855;
    const double znunu_mumu_ratio = 5.942;
    const double znunu_ee_ratio   = 5.942;

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_400to600"]};
        fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
    }
    else
    {
        if(ss[dataSets] != ss.null())
        {
            fileMap[dataSets] = {ss[dataSets]};
            for(const auto& colls : ss[dataSets].getCollections())
            {
                fileMap[colls] = {ss[dataSets]};
            }
        }
        else if(sc[dataSets] != sc.null())
        {
            fileMap[dataSets] = {sc[dataSets]};
            int i = 0;
            for(const auto& fs : sc[dataSets])
            {
                fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
            }
        }
    }


    vector<Plotter::HistSummary> vh;

    

    Plotter::DatasetSummary dsDY_ll_inc(    "DY#rightarrow#mu#mu Inc",               fileMap["IncDY"], "pdgIdZDec=13;passMuZinvSel", "");

    Plotter::DatasetSummary dsDY_ll_noSel(      "DY#rightarrow#mu#mu no sel",            fileMap["DYJetsToLL"],   "pdgIdZDec=13", "");
    Plotter::DatasetSummary dsDY_ll_gen_noSel(  "DY#rightarrow#mu#mu Gen no Sel",        fileMap["DYJetsToLL"],   "pdgIdZDec=13", "");

    Plotter::DatasetSummary dsTT_mm("t#bar{t}#rightarrow#mu#mu", fileMap["TTbarDiLep"], "", "");
    Plotter::DatasetSummary dsTT_ee("t#bar{t}#rightarrowee",     fileMap["TTbarDiLep"], "", "");
    Plotter::DatasetSummary dsTT_em("t#bar{t}#rightarrowe#mu",   fileMap["TTbarDiLep"], "", "");

    Plotter::DatasetSummary dsTT_mm_lepSel("t#bar{t}#rightarrow#mu#mu", fileMap["TTbarDiLep"], "passMuZinvSel", "");
    Plotter::DatasetSummary dsTT_ee_lepSel("t#bar{t}#rightarrowee",     fileMap["TTbarDiLep"], "passElecZinvSel", "");
    Plotter::DatasetSummary dsTT_em_lepSel("t#bar{t}#rightarrowe#mu",   fileMap["TTbarDiLep"], "passElMuZinvSel", "");

    Plotter::DatasetSummary dsDY_ee(               "DY#rightarrowee",                 fileMap["DYJetsToLL"],   "pdgIdZDec=11;passElecZinvSel", "");
    Plotter::DatasetSummary dsDY_ee_scaled(        "DY#rightarrowee",                 fileMap["DYJetsToLL"],   "pdgIdZDec=11;passElecZinvSel", "",                        znunu_ee_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ee_zAcc_scaled(   "DY#rightarrowee no e, Z eff+acc", fileMap["DYJetsToLL"],   "pdgIdZDec=11;passElecZinvSel", "zEffWgtElec;zAccWgtElec", znunu_ee_ratio / zAcc);
    Plotter::DatasetSummary dsDY_ee_zWeight_scaled("DY#rightarrowee no e, Z eff",     fileMap["DYJetsToLL"],   "pdgIdZDec=11;passElecZinvSel", "zEffWgtElec",             znunu_ee_ratio / zAcc);

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

    Plotter::DatasetSummary dsDY_nunu(            "Z#rightarrow#nu#nu",                fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_SB0b(       "Z#rightarrow#nu#nu, N(b) = 0",      fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_nunu_SBMb(       "Z#rightarrow#nu#nu, Direct MC",     fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    Plotter::DatasetSummary dsDY_test(            "test",                              fileMap["ZJetsToNuNu"], "passLeptVeto", "");

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
    vh.push_back(PHS("cleanJetPt",                 {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                      100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJetPt_baseline",        {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",      100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJetPt_baselineNoTag",   {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ",               {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                      100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ_baseline",      {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",      100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("mindPhigenZJ_baselineNoTag", {PDC("single", "mindPhiMetJ",      {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt",               {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                           100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt_baseline",      {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",           100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet1pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",      100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt",               {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                           100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt_baseline",      {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",           100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet2pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",      100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt",               {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",                           100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt_baseline",      {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",           100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
    vh.push_back(PHS("cleanJet3pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",      100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));

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

    vh.push_back(PHS("ttbar_mht",         {PDC("single", "mht",                    {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "M(H_{t}) [GeV]",          "Norm Events"));
    vh.push_back(PHS("ttbar_mt2",         {PDC("single", "best_had_brJet_MT2Zinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MT2 [GeV]",               "Norm Events"));
    vh.push_back(PHS("ttbar_met",         {PDC("single", "met",                    {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MET [GeV]",               "Norm Events"));
    vh.push_back(PHS("ttbar_cleanMet",    {PDC("single", "cleanMetPt",             {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "MET [GeV]",               "Norm Events"));
    vh.push_back(PHS("ttbar_bestZPt",     {PDC("single", "bestRecoZPt",            {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]",           "Norm Events"));
    vh.push_back(PHS("ttbar_ht",          {PDC("single", "HT",                     {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",             "Norm Events"));
    vh.push_back(PHS("ttbar_cleanht",     {PDC("single", "HTZinv",                 {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",             "Norm Events"));
    vh.push_back(PHS("ttbar_nMuons",      {PDC("single", "cutMuVec(size)",         {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(#mu)",                  "Events"));
    vh.push_back(PHS("ttbar_jetPt",       {PDC("single", "jetsLVec(pt)",           {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Jet p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_cleannJets",  {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 20,  0, 20,    true,  false, "N(jets)",                 "Events"));
    vh.push_back(PHS("ttbar_muPt",        {PDC("single", "cutMuVec(pt)",           {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_mu1Pt",       {PDC("single", "cutMuVec[0](pt)",        {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_mu2Pt",       {PDC("single", "cutMuVec[1](pt)",        {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_elecPt",      {PDC("single", "cutElecVec(pt)",         {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_elec1Pt",     {PDC("single", "cutElecVec[0](pt)",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_elec2Pt",     {PDC("single", "cutElecVec[1](pt)",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",         "Norm Events"));
    vh.push_back(PHS("ttbar_nb",          {PDC("single", "cntCSVSZinv",            {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(b)",                    "Events"));
    vh.push_back(PHS("ttbar_nt",          {PDC("single", "nTopCandSortedCnt",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10,    true,  false, "N(t)",                    "Events"));

    vh.push_back(PHS("ttbar_baseline_mht",         {PDC("single", "mht",                    {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "M(H_{t}) [GeV]",  "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_mt2",         {PDC("single", "best_had_brJet_MT2Zinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "MT2 [GeV]",       "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_met",         {PDC("single", "met",                    {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "MET [GeV]",       "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_cleanMet",    {PDC("single", "cleanMetPt",             {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "MET [GeV]",       "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_bestZPt",     {PDC("single", "bestRecoZPt",            {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]",   "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_ht",          {PDC("single", "HT",                     {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",     "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_cleanht",     {PDC("single", "HTZinv",                 {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 2000,  true,  true,  "H_{T} [GeV]",     "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_nMuons",      {PDC("single", "cutMuVec(size)",         {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10,    true,  false, "N(#mu)",          "Events"));
    vh.push_back(PHS("ttbar_baseline_jetPt",       {PDC("single", "jetsLVec(pt)",           {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "Jet p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_cleannJets",  {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 20,  0, 20,    true,  false, "N(jets)",         "Events"));
    vh.push_back(PHS("ttbar_baseline_muPt",        {PDC("single", "cutMuVec(pt)",           {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_mu1Pt",       {PDC("single", "cutMuVec[0](pt)",        {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_mu2Pt",       {PDC("single", "cutMuVec[1](pt)",        {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_elecPt",      {PDC("single", "cutElecVec(pt)",         {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_elec1Pt",     {PDC("single", "cutElecVec[0](pt)",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_elec2Pt",     {PDC("single", "cutElecVec[1](pt)",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("ttbar_baseline_nb",          {PDC("single", "cntCSVSZinv",            {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10,    true,  false, "N(b)",            "Events"));
    vh.push_back(PHS("ttbar_baseline_nt",          {PDC("single", "nTopCandSortedCnt",      {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10,    true,  false, "N(t)",            "Events"));

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

    Plotter::DataCollection scaled_cleanht(     "single", {{"HTZinv", dsDY_nunu},                 {"HTZinv", dsDY_ll_zAcc_scaled}                 });
    Plotter::DataCollection scaled_cleanmt2(    "single", {{"best_had_brJet_MT2Zinv", dsDY_nunu}, {"best_had_brJet_MT2Zinv", dsDY_ll_zAcc_scaled} });
    Plotter::DataCollection scaled_cleanmht(    "single", {{"cleanMHt", dsDY_nunu},               {"cleanMHt", dsDY_ll_zAcc_scaled}               });
    Plotter::DataCollection scaled_cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu},            {"cleanMHtPhi", dsDY_ll_zAcc_scaled}            });
    Plotter::DataCollection scaled_cleanMet(    "single", {{"cleanMetPt", dsDY_nunu},             {"cleanMetPt", dsDY_ll_zAcc_scaled}             });
    Plotter::DataCollection scaled_nCleanJet(   "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nTop(        "single", {{"nTopCandSortedCntZinv", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nBottom(     "single", {{"cntCSVSZinv", dsDY_nunu},            {"cntCSVSZinv", dsDY_ll_zAcc_scaled}            });

    Plotter::DataCollection scaled_elec_ht(   "single", "ht",  {{dsDY_nunu},            {dsDY_ee_zWeight_scaled},               {dsDY_ee_zAcc_scaled}               });
    Plotter::DataCollection scaled_elec_mht(  "single", "mht", {{dsDY_nunu},            {dsDY_ee_zWeight_scaled},               {dsDY_ee_zAcc_scaled}               });
    Plotter::DataCollection scaled_elec_met(  "single", "met", {{dsDY_nunu},            {dsDY_ee_zWeight_scaled},               {dsDY_ee_zAcc_scaled}               });

    Plotter::DataCollection scaled_elec_cleanht(     "single", {{"HTZinv",                 dsDY_nunu},                 {"HTZinv", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_cleanmt2(    "single", {{"best_had_brJet_MT2Zinv", dsDY_nunu}, {"best_had_brJet_MT2Zinv", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_cleanmht(    "single", {{"cleanMHt",               dsDY_nunu},               {"cleanMHt", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_cleanmhtphi( "single", {{"cleanMHtPhi",            dsDY_nunu},            {"cleanMHtPhi", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_cleanMet(    "single", {{"cleanMetPt",             dsDY_nunu},             {"cleanMetPt", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_nCleanJet(   "single", {{"cntNJetsPt30Eta24Zinv",  dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_nTop(        "single", {{"nTopCandSortedCntZinv",  dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ee_zAcc_scaled} });
    Plotter::DataCollection scaled_elec_nBottom(     "single", {{"cntCSVSZinv",            dsDY_nunu},            {"cntCSVSZinv", dsDY_ee_zAcc_scaled} });

    Plotter::DataCollection scaled_nSearchBin(     "single", {{"nSearchBin", dsDY_nunu},             {"nSearchBin", dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nSearchBin_stat("single", {{"nSearchBin", dsDY_ll},               {"nb0BinsNW",  dsDY_ll}              });
    Plotter::DataCollection scaled_nSearchBinNb0(  "single", {{"nSearchBin", dsDY_nunu},             {"nb0Bins",    dsDY_ll_zAcc_scaled}  });
    Plotter::DataCollection scaled_nSearchBinNb0ZZ("single", {{"nSearchBin", dsDY_nunu_SBMb},        {"nb0Bins",    dsDY_nunu_SB0b}       });

    Plotter::DataCollection scaled_Elec_nSearchBin(     "single", {{"nSearchBin", dsDY_nunu},             {"nSearchBin", dsDY_ee_zAcc_scaled}  });
    Plotter::DataCollection scaled_Elec_nSearchBin_stat("single", {{"nSearchBin", dsDY_ee},               {"nb0BinsNW",  dsDY_ee}              });
    Plotter::DataCollection scaled_Elec_nSearchBinNb0(  "single", {{"nSearchBin", dsDY_nunu},             {"nb0Bins",    dsDY_ee_zAcc_scaled}  });
    Plotter::DataCollection scaled_Elec_nSearchBinNb0ZZ("single", {{"nSearchBin", dsDY_nunu_SBMb},        {"nb0Bins",    dsDY_nunu_SB0b}       });

    Plotter::DataCollection scaled_nSearchBin_Sig_mm("stack", {{"nSearchBin", dsDY_ll}, {"nSearchBin",  dsTT_mm_lepSel} });
    Plotter::DataCollection scaled_nSearchBin_Con_mm("stack", {{"nb0BinsNW", dsDY_ll},  {"nb0BinsNW",   dsTT_mm_lepSel} });
    Plotter::DataCollection scaled_nSearchBin_Sig_ee("stack", {{"nSearchBin", dsDY_ee}, {"nSearchBin",  dsTT_ee_lepSel} });
    Plotter::DataCollection scaled_nSearchBin_Con_ee("stack", {{"nb0BinsNW", dsDY_ee},  {"nb0BinsNW",   dsTT_ee_lepSel} });


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
    vh.push_back(PHS("nSearchBin_stat",          {scaled_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_stat_log",      {scaled_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",               45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
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

    vh.push_back(PHS("elec_nSearchBin",               {scaled_Elec_nSearchBin},      {2, 1}, "passBaselineZinv",                                      45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBin_log",           {scaled_Elec_nSearchBin},      {2, 1}, "passBaselineZinv",                                      45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBin_stat",          {scaled_Elec_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBin_stat_log",      {scaled_Elec_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBinnb0",            {scaled_Elec_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBinnb0_log",        {scaled_Elec_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBinnb0ZZ",          {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBinnb0ZZ_log",      {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("elec_nSearchBin_pull",          {scaled_Elec_nSearchBin},      {2, 1}, "passBaselineZinv",                                      45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("elec_nSearchBin_pull_log",      {scaled_Elec_nSearchBin},      {2, 1}, "passBaselineZinv",                                      45,  0,     45,   true,  false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("elec_nSearchBinnb0_pull",       {scaled_Elec_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("elec_nSearchBinnb0_pull_log",   {scaled_Elec_nSearchBinNb0},   {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("elec_nSearchBinnb0ZZ_pull",     {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("elec_nSearchBinnb0ZZ_pull_log", {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", false));

    vh.push_back(PHS("nSearchBin_Sig_mm",          {scaled_nSearchBin_Sig_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Sig_mm_log",      {scaled_nSearchBin_Sig_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Con_mm",          {scaled_nSearchBin_Con_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Con_mm_log",      {scaled_nSearchBin_Con_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Sig_ee",          {scaled_nSearchBin_Sig_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Sig_ee_log",      {scaled_nSearchBin_Sig_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Con_ee",          {scaled_nSearchBin_Con_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("nSearchBin_Con_ee_log",      {scaled_nSearchBin_Con_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0",        45,  0,     45,   true,  false,  "Search Bin",     "Events", true));

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

    vh.push_back(PHS("elec_mT2Zinv",                {scaled_elec_cleanmt2},      {2, 1}, "",                                                             100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_baseline_mT2Zinv",       {scaled_elec_cleanmt2},      {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_baselineNoTag_mT2Zinv",  {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_0_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_1_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_2_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_3_mT2Zinv",           {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_0_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_0_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_0_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_0_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_1_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_1_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_1_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_1_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_2_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_2_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_2_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_2_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_3_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_3_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_3_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_3_mT2Zinv",      {scaled_elec_cleanmt2},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  false,  "M_{T2} [GeV]",   "Events"));
    vh.push_back(PHS("elec_cleanht",                {scaled_elec_cleanht},       {2, 1}, "",                                                             100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_baseline_cleanht",       {scaled_elec_cleanht},       {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_baselineNoTag_cleanht",  {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_baseline_ht",            {scaled_elec_ht},            {2, 1}, "passBaselineZinv",                                             50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_0_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_1_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_2_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_3_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nt_0_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nt_1_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nt_2_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nt_3_cleanht",           {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_0_nt_0_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_1_nt_0_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_2_nt_0_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_3_nt_0_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_0_nt_1_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_1_nt_1_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_2_nt_1_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_3_nt_1_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_0_nt_2_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_1_nt_2_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_2_nt_2_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_3_nt_2_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_0_nt_3_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_1_nt_3_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_2_nt_3_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_nb_3_nt_3_cleanht",      {scaled_elec_cleanht},       {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     2000, true,  false,  "H_{T} [GeV]",    "Events"));
    vh.push_back(PHS("elec_cleanmht",               {scaled_elec_cleanmht},      {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_baseline_cleanmht",      {scaled_elec_cleanmht},      {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_baselineNoTag_cleanmht", {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_baseline_mht",           {scaled_elec_mht},           {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_0_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_1_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_2_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nt_3_cleanmht",          {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_0_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_0_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_0_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_0_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_1_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_1_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_1_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_1_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_2_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_2_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_2_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_2_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_0_nt_3_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_1_nt_3_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_2_nt_3_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_nb_3_nt_3_cleanmht",     {scaled_elec_cleanmht},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MH_{T} [GeV]",   "Events"));
    vh.push_back(PHS("elec_cleanmet",               {scaled_elec_cleanMet},      {2, 1}, "",                                                             100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_cleanmet_gr",            {scaled_elec_cleanMet},      {2, 1}, "",                                                             200, 0,     3000, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_baseline_cleanmet",      {scaled_elec_cleanMet},      {2, 1}, "passBaselineZinv",                                             50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_baselineNoTag_cleanmet", {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv",                                        50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_0_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_1_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_2_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_3_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nt_0_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0",                  100, 0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nt_1_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=1",                50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nt_2_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv=2",                25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nt_3_cleanmet",          {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;nTopCandSortedCntZinv>2",                15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_0_nt_0_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0",    50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_1_nt_0_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0",    25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_2_nt_0_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0",    15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_3_nt_0_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0",    10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_0_nt_1_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_1_nt_1_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_2_nt_1_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_3_nt_1_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_0_nt_2_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_1_nt_2_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_2_nt_2_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_3_nt_2_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_0_nt_3_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2",  50,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_1_nt_3_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2",  25,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_2_nt_3_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2",  15,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nb_3_nt_3_cleanmet",     {scaled_elec_cleanMet},      {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2",  10,  0,     1500, true,  false,  "MET [GeV]",      "Events"));
    vh.push_back(PHS("elec_nTop",                   {scaled_elec_nTop},          {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_baseline_nTop",          {scaled_elec_nTop},          {2, 1}, "passBaselineZinv",                                             10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_baselineNoTag_nTop",     {scaled_elec_nTop},          {2, 1}, "passBaselineNoTagZinv",                                        10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_nb_0_nTop",              {scaled_elec_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_nb_1_nTop",              {scaled_elec_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=1",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_nb_2_nTop",              {scaled_elec_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_nb_3_nTop",              {scaled_elec_nTop},          {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv>2",                          10,  0,       10, true,  false,  "N(t)",           "Events"));
    vh.push_back(PHS("elec_cleannJet",              {scaled_elec_nCleanJet},     {2, 1}, "",                                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("elec_baseline_cleannJet",     {scaled_elec_nCleanJet},     {2, 1}, "passBaselineZinv",                                             20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("elec_baselineNoTag_cleannJet",{scaled_elec_nCleanJet},     {2, 1}, "passBaselineNoTagZinv",                                        20,  0,       20, true,  false,  "N(jet)",         "Events"));
    vh.push_back(PHS("elec_nBottom",                {scaled_elec_nBottom},       {2, 1}, "",                                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("elec_baseline_nBottom",       {scaled_elec_nBottom},       {2, 1}, "passBaselineZinv",                                             10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("elec_baselineNoTAg_nBottom",  {scaled_elec_nBottom},       {2, 1}, "passBaselineNoTagZinv",                                        10,  0,       10, true,  false,  "N(b)",           "Events"));
    vh.push_back(PHS("elec_cleanmhtphi",            {scaled_elec_cleanmhtphi},   {2, 1}, "",                                                             100, -3.14, 3.14, false, false,  "#phi(MH_{T})",   "Events"));
    vh.push_back(PHS("elec_baseline_cleanmhtphi",   {scaled_elec_cleanmhtphi},   {2, 1}, "passBaselineZinv",                                             100, -3.14, 3.14, false, false,  "#phi(MH_{T})",   "Events"));

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

    

    // -----------------
    // - Data/MC plots -
    // -----------------

    // Datasetsummaries we are using
    // no weight
    Plotter::DatasetSummary dsData_SingleMuon("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "");
    Plotter::DatasetSummary dsData_DoubleEG(  "Data",       fileMap["Data_DoubleEG"],   "passElecTrigger", "");
    Plotter::DatasetSummary dsDY(             "DY",         fileMap["DYJetsToLL"],      "",                "nJetWgtDYZ");
    Plotter::DatasetSummary dsDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",       "nJetWgtDYZ");
    Plotter::DatasetSummary dstt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",                "nJetWgtDYZ");
    Plotter::DatasetSummary dstW(             "single top", fileMap["tW"],              "",                "nJetWgtDYZ");
    Plotter::DatasetSummary dsttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",                "nJetWgtDYZ;genWeight");
    Plotter::DatasetSummary dsVV(             "Diboson",    fileMap["Diboson"],         "",                "nJetWgtDYZ");
    Plotter::DatasetSummary dsRare(           "Rare",       fileMap["Rare"],            "",                "nJetWgtDYZ;genWeight");
    std::vector<Plotter::DatasetSummary> stack_MC = {dsDY, dsDYInc, dstt2l, dstW, dsttZ, dsVV, dsRare};
    // 1 fake b weight
    Plotter::DatasetSummary dsData_SingleMuon_w1b("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "weight1fakebComb");
    Plotter::DatasetSummary dsData_DoubleEG_w1b(  "Data",       fileMap["Data_DoubleEG"],   "passElecTrigger", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_w1b(             "DY",         fileMap["DYJetsToLL"],      "", "weight1fakebComb");
    Plotter::DatasetSummary dsDYInc_w1b(          "DY HT<100",  fileMap["IncDY"],           "genHT<100", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_w1b(           "t#bar{t}",   fileMap["TTbarNoHad"],      "", "weight1fakebComb");
    Plotter::DatasetSummary dstW_w1b(             "single top", fileMap["tW"],              "", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_w1b(            "t#bar{t}Z",  fileMap["TTZ"],             "", "genWeight;weight1fakebComb");
    Plotter::DatasetSummary dsVV_w1b(             "Diboson",    fileMap["Diboson"],         "", "weight1fakebComb");
    Plotter::DatasetSummary dsRare_w1b(           "Rare",       fileMap["Rare"],            "", "genWeight;weight1fakebComb");
    std::vector<Plotter::DatasetSummary> stack_MC_w1b = {dsDY_w1b, dsDYInc_w1b, dstt2l_w1b, dstW_w1b, dsttZ_w1b, dsVV_w1b, dsRare_w1b};
    // 2 fake b weight
    Plotter::DatasetSummary dsData_SingleMuon_w2b("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "weight2fakebComb");
    Plotter::DatasetSummary dsData_DoubleEG_w2b(  "Data",       fileMap["Data_DoubleEG"],   "passElecTrigger", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_w2b(             "DY",         fileMap["DYJetsToLL"],      "", "weight2fakebComb");
    Plotter::DatasetSummary dsDYInc_w2b(          "DY HT<100",  fileMap["IncDY"],           "genHT<100", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_w2b(           "t#bar{t}",   fileMap["TTbarNoHad"],      "", "weight2fakebComb");
    Plotter::DatasetSummary dstW_w2b(             "single top", fileMap["tW"],              "", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_w2b(            "t#bar{t}Z",  fileMap["TTZ"],             "", "genWeight;weight2fakebComb");
    Plotter::DatasetSummary dsVV_w2b(             "Dibson",     fileMap["Diboson"],         "", "weight2fakebComb");
    Plotter::DatasetSummary dsRare_w2b(           "Rare",       fileMap["Rare"],            "", "genWeight;weight2fakebComb");
    std::vector<Plotter::DatasetSummary> stack_MC_w2b = {dsDY_w2b, dsDYInc_w2b, dstt2l_w2b, dstW_w2b, dsttZ_w2b, dsVV_w2b, dsRare_w2b};
    // 3 fake b weight
    Plotter::DatasetSummary dsData_SingleMuon_w3b("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "weight3fakebComb");
    Plotter::DatasetSummary dsData_DoubleEG_w3b(  "Data",       fileMap["Data_DoubleEG"],   "passElecTrigger", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_w3b(             "DY",         fileMap["DYJetsToLL"],      "", "weight3fakebComb");
    Plotter::DatasetSummary dsDYInc_w3b(          "DY HT<100",  fileMap["IncDY"],           "genHT<100", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_w3b(           "t#bar{t}",   fileMap["TTbarNoHad"],      "", "weight3fakebComb");
    Plotter::DatasetSummary dstW_w3b(             "single top", fileMap["tW"],              "", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_w3b(            "t#bar{t}Z",  fileMap["TTZ"],             "", "genWeight;weight3fakebComb");
    Plotter::DatasetSummary dsVV_w3b(             "Diboson",    fileMap["Diboson"],         "", "weight3fakebComb");
    Plotter::DatasetSummary dsRare_w3b(           "Rare",       fileMap["Rare"],            "", "genWeight;weight3fakebComb");
    std::vector<Plotter::DatasetSummary> stack_MC_w3b = {dsDY_w3b, dsDYInc_w3b, dstt2l_w3b, dstW_w3b, dsttZ_w3b, dsVV_w3b, dsRare_w3b};

    // Collections for all variables, no cuts applied yet
    // met
    Plotter::DataCollection dcData_SingleMuon_met("data",   "cleanMetPt", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_met(  "data",   "cleanMetPt", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_met(             "stack",  "cleanMetPt", stack_MC);
    // ntops
    Plotter::DataCollection dcData_SingleMuon_nt("data",   "nTopCandSortedCntZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_nt(  "data",   "nTopCandSortedCntZinv", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_nt(             "stack",  "nTopCandSortedCntZinv", stack_MC);
    // ntops 1b fake
    Plotter::DataCollection dcData_SingleMuon_nt1b("data",   "nTopCandSortedCntZinv1b", {dsData_SingleMuon_w1b});
    Plotter::DataCollection dcData_DoubleEG_nt1b(  "data",   "nTopCandSortedCntZinv1b", {dsData_DoubleEG_w1b});
    Plotter::DataCollection dcMC_nt1b(             "stack",  "nTopCandSortedCntZinv1b", stack_MC_w1b);
    // ntops 2b fake
    Plotter::DataCollection dcData_SingleMuon_nt2b("data",   "nTopCandSortedCntZinv2b", {dsData_SingleMuon_w2b});
    Plotter::DataCollection dcData_DoubleEG_nt2b(  "data",   "nTopCandSortedCntZinv2b", {dsData_DoubleEG_w2b});
    Plotter::DataCollection dcMC_nt2b(             "stack",  "nTopCandSortedCntZinv2b", stack_MC_w2b);
    // ntops 3b fake
    Plotter::DataCollection dcData_SingleMuon_nt3b("data",   "nTopCandSortedCntZinv3b", {dsData_SingleMuon_w3b});
    Plotter::DataCollection dcData_DoubleEG_nt3b(  "data",   "nTopCandSortedCntZinv3b", {dsData_DoubleEG_w3b});
    Plotter::DataCollection dcMC_nt3b(             "stack",  "nTopCandSortedCntZinv3b", stack_MC_w3b);
    // MT2
    Plotter::DataCollection dcData_SingleMuon_mt2("data",   "best_had_brJet_MT2Zinv", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mt2(  "data",   "best_had_brJet_MT2Zinv", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mt2(             "stack",  "best_had_brJet_MT2Zinv", stack_MC);
    // MT2 1b fake
    Plotter::DataCollection dcData_SingleMuon_mt21b("data",   "best_had_brJet_MT2Zinv1b", {dsData_SingleMuon_w1b});
    Plotter::DataCollection dcData_DoubleEG_mt21b(  "data",   "best_had_brJet_MT2Zinv1b", {dsData_DoubleEG_w1b});
    Plotter::DataCollection dcMC_mt21b(             "stack",  "best_had_brJet_MT2Zinv1b", stack_MC_w1b);
    // MT2 2b fake
    Plotter::DataCollection dcData_SingleMuon_mt22b("data",   "best_had_brJet_MT2Zinv2b", {dsData_SingleMuon_w2b});
    Plotter::DataCollection dcData_DoubleEG_mt22b(  "data",   "best_had_brJet_MT2Zinv2b", {dsData_DoubleEG_w2b});
    Plotter::DataCollection dcMC_mt22b(             "stack",  "best_had_brJet_MT2Zinv2b", stack_MC_w2b);
    // MT2 3b fake
    Plotter::DataCollection dcData_SingleMuon_mt23b("data",   "best_had_brJet_MT2Zinv3b", {dsData_SingleMuon_w3b});
    Plotter::DataCollection dcData_DoubleEG_mt23b(  "data",   "best_had_brJet_MT2Zinv3b", {dsData_DoubleEG_w3b});
    Plotter::DataCollection dcMC_mt23b(             "stack",  "best_had_brJet_MT2Zinv3b", stack_MC_w3b);
    // nb
    Plotter::DataCollection dcData_SingleMuon_nb("data",   "cntCSVSZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_nb(  "data",   "cntCSVSZinv", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_nb(             "stack",  "cntCSVSZinv", stack_MC);
    // nj
    Plotter::DataCollection dcData_SingleMuon_nj("data",   "cntNJetsPt30Eta24Zinv", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_nj(  "data",   "cntNJetsPt30Eta24Zinv", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_nj(             "stack",  "cntNJetsPt30Eta24Zinv", stack_MC);
    // ht
    Plotter::DataCollection dcData_SingleMuon_ht("data",   "HTZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_ht(  "data",   "HTZinv", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_ht(             "stack",  "HTZinv", stack_MC);
    // mht
    Plotter::DataCollection dcData_SingleMuon_mht("data",   "cleanMHt", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mht(  "data",   "cleanMHt", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mht(             "stack",  "cleanMHt", stack_MC);
    // jpt
    Plotter::DataCollection dcData_SingleMuon_jpt("data",   "jetsLVecLepCleaned(pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_jpt(  "data",   "jetsLVecLepCleaned(pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_jpt(             "stack",  "jetsLVecLepCleaned(pt)", stack_MC);
    // j1pt
    Plotter::DataCollection dcData_SingleMuon_j1pt("data",   "jetsLVecLepCleaned[0](pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_j1pt(  "data",   "jetsLVecLepCleaned[0](pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_j1pt(             "stack",  "jetsLVecLepCleaned[0](pt)", stack_MC);
    // j2pt
    Plotter::DataCollection dcData_SingleMuon_j2pt("data",   "jetsLVecLepCleaned[1](pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_j2pt(  "data",   "jetsLVecLepCleaned[1](pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_j2pt(             "stack",  "jetsLVecLepCleaned[1](pt)", stack_MC);
    // j3pt
    Plotter::DataCollection dcData_SingleMuon_j3pt("data",   "jetsLVecLepCleaned[2](pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_j3pt(  "data",   "jetsLVecLepCleaned[2](pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_j3pt(             "stack",  "jetsLVecLepCleaned[2](pt)", stack_MC);
    // mupt
    Plotter::DataCollection dcData_SingleMuon_mupt("data",   "cutMuVec(pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mupt(  "data",   "cutMuVec(pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mupt(             "stack",  "cutMuVec(pt)", stack_MC);
    // mu1pt
    Plotter::DataCollection dcData_SingleMuon_mu1pt("data",   "cutMuPt1", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mu1pt(  "data",   "cutMuPt1", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mu1pt(             "stack",  "cutMuPt1", stack_MC);
    // mu2pt
    Plotter::DataCollection dcData_SingleMuon_mu2pt("data",   "cutMuPt2", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mu2pt(  "data",   "cutMuPt2", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mu2pt(             "stack",  "cutMuPt2", stack_MC);
    // elpt
    Plotter::DataCollection dcData_SingleMuon_elpt("data",   "cutElecVec(pt)", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_elpt(  "data",   "cutElecVec(pt)", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_elpt(             "stack",  "cutElecVec(pt)", stack_MC);
    // el1pt
    Plotter::DataCollection dcData_SingleMuon_el1pt("data",   "cutElecPt1", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_el1pt(  "data",   "cutElecPt1", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_el1pt(             "stack",  "cutElecPt1", stack_MC);
    // el2pt
    Plotter::DataCollection dcData_SingleMuon_el2pt("data",   "cutElecPt2", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_el2pt(  "data",   "cutElecPt2", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_el2pt(             "stack",  "cutElecPt2", stack_MC);
    // mll
    Plotter::DataCollection dcData_SingleMuon_mll("data",   "bestRecoZM", {dsData_SingleMuon});
    Plotter::DataCollection dcData_DoubleEG_mll(  "data",   "bestRecoZM", {dsData_DoubleEG});
    Plotter::DataCollection dcMC_mll(             "stack",  "bestRecoZM", stack_MC);


    // Define the cutlevels
    std::string s_nosel                   = "passNoiseEventFilterZinv";
    // dimu				  
    std::string s_2mu                     = "passNoiseEventFilterZinv;passDiMuSel";
    std::string s_muZinv                  = "passNoiseEventFilterZinv;passMuZinvSel";
    std::string s_muZinv_ht200            = "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200";
    std::string s_muZinv_loose0           = "passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_muZinv_loose50          = "passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_muZinv_loose100         = "passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_muZinv_loose200         = "passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    std::string s_muZinv_blnotagmt2       = "passMuZinvSel;passBaselineNoTagMT2Zinv";
    std::string s_muZinv_bl               = "passMuZinvSel;passBaselineZinv";
    std::string s_muZinv_0b               = "passMuZinvSel;cntCSVSZinv=0";
    std::string s_muZinv_0b_ht200         = "passMuZinvSel;cntCSVSZinv=0;HTZinv>200";
    std::string s_muZinv_0b_loose0        = "passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_muZinv_0b_loose50       = "passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_muZinv_0b_loose100      = "passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_muZinv_0b_loose200      = "passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    std::string s_muZinv_0b_blnotagmt2    = "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv";
    std::string s_muZinv_0b_blnotag       = "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv";
    std::string s_muZinv_g1b_loose0       = "passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_muZinv_g1b_loose50      = "passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_muZinv_g1b_loose100     = "passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_muZinv_g1b_loose200     = "passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    // diele
    std::string s_2el                     = "passNoiseEventFilterZinv;passDiElecSel";
    std::string s_elZinv                  = "passNoiseEventFilterZinv;passElecZinvSel";
    std::string s_elZinv_ht200            = "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200";
    std::string s_elZinv_blnotagmt2       = "passElecZinvSel;passBaselineNoTagMT2Zinv";
    std::string s_elZinv_bl               = "passElecZinvSel;passBaselineZinv";
    std::string s_elZinv_0b               = "passElecZinvSel;cntCSVSZinv=0";
    std::string s_elZinv_0b_ht200         = "passElecZinvSel;cntCSVSZinv=0;HTZinv>200";
    //std::string s_elZinv_0b_loose0        = "passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv";
    //std::string s_elZinv_0b_loose50       = "passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    //std::string s_elZinv_0b_loose100      = "passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    //std::string s_elZinv_0b_loose200      = "passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    std::string s_elZinv_0b_blnotagmt2    = "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv";
    std::string s_elZinv_0b_blnotag       = "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv";
    // elmu
    std::string s_elmu                    = "passNoiseEventFilterZinv;passElMuSel";
    std::string s_elmuZinv                = "passNoiseEventFilterZinv;passElMuZinvSel";
    std::string s_elmuZinv_ht200          = "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200";
    std::string s_elmuZinv_loose0         = "passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_elmuZinv_loose50        = "passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_elmuZinv_loose100       = "passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_elmuZinv_loose200       = "passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    std::string s_elmuZinv_blnotagmt2     = "passElMuZinvSel;passBaselineNoTagMT2Zinv";
    std::string s_elmuZinv_bl             = "passElMuZinvSel;passBaselineZinv";
    std::string s_elmuZinv_0b             = "passElMuZinvSel;cntCSVSZinv=0";
    std::string s_elmuZinv_0b_ht200       = "passElMuZinvSel;cntCSVSZinv=0;HTZinv>200";
    std::string s_elmuZinv_0b_loose0      = "passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_elmuZinv_0b_loose50     = "passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_elmuZinv_0b_loose100    = "passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_elmuZinv_0b_loose200    = "passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    std::string s_elmuZinv_0b_blnotagmt2  = "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv";
    std::string s_elmuZinv_0b_blnotag     = "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv";
    std::string s_elmuZinv_g1b_loose0     = "passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_elmuZinv_g1b_loose50    = "passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_elmuZinv_g1b_loose100   = "passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_elmuZinv_g1b_loose200   = "passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";
    
    //interlude for MC checks

    std::string s_znunu_baseline      = "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0";  //THIS IS NOT PERFECT, clearly nTopCandSortedCntZinv and MT2 cuts need to be specialized by sample
    std::string s_znunu_loose0        = "passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv";
    std::string s_znunu_loose50       = "passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50";
    std::string s_znunu_loose100      = "passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100";
    std::string s_znunu_loose200      = "passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv";

    Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt1b( "Z#rightarrow#nu#nu, N(b) = 0, 1 fake",     fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "nJet1bfakeWgt");
    Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt2b( "Z#rightarrow#nu#nu, N(b) = 0, 2 fake",     fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "nJet2bfakeWgt");
    Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt3b( "Z#rightarrow#nu#nu, N(b) = 0, 3 fake",     fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "nJet3bfakeWgt");
    Plotter::DatasetSummary dsDY_nunu_SB0b_2(     "Z#rightarrow#nu#nu, N(b) = 0",             fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_nunu_SB1b(       "Z#rightarrow#nu#nu, Direct MC, N(b) = 1",  fileMap["ZJetsToNuNu"], "cntCSVSZinv=1", "");
    Plotter::DatasetSummary dsDY_nunu_SB2b(       "Z#rightarrow#nu#nu, Direct MC, N(b) = 2",  fileMap["ZJetsToNuNu"], "cntCSVSZinv=2", "");
    Plotter::DatasetSummary dsDY_nunu_SB3b(       "Z#rightarrow#nu#nu, Direct MC, N(b) >= 3", fileMap["ZJetsToNuNu"], "cntCSVSZinv>2", "");

    Plotter::DatasetSummary dsDY_nunu_njet(       "Z#rightarrow#nu#nu Data/MC weight",        fileMap["ZJetsToNuNu"], "passLeptVeto", "nJetWgtDYZ");

    //norm plots 
    Plotter::DataCollection dcMC_nunu_nj_1b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB1b});
    Plotter::DataCollection dcMC_nunu_nj_2b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB2b});
    Plotter::DataCollection dcMC_nunu_nj_3b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB3b});

    //1 fake b
    Plotter::DataCollection dcMC_nunu_Wgt1b_met("single", "cleanMetPt",                {dsDY_nunu_SB0b_Wgt1b,                             dsDY_nunu_SB1b});
    Plotter::DataCollection dcMC_nunu_Wgt1b_ht( "single", "HTZinv",                    {dsDY_nunu_SB0b_Wgt1b,                             dsDY_nunu_SB1b});
    Plotter::DataCollection dcMC_nunu_Wgt1b_nb( "single", "cntCSVSZinv",               {dsDY_nunu_SB0b_Wgt1b,                             dsDY_nunu_SB1b});
    Plotter::DataCollection dcMC_nunu_Wgt1b_nj( "single", "cntNJetsPt30Eta24Zinv",     {dsDY_nunu_SB0b_Wgt1b,                             dsDY_nunu_SB1b});
    Plotter::DataCollection dcMC_nunu_Wgt1b_nt( "single", {{ "nTopCandSortedCntZinv1b", dsDY_nunu_SB0b_Wgt1b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB1b}});
    Plotter::DataCollection dcMC_nunu_Wgt1b_mt2("single", {{"best_had_brJet_MT2Zinv1b", dsDY_nunu_SB0b_Wgt1b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB1b}});
    //2 fake b
    Plotter::DataCollection dcMC_nunu_Wgt2b_met("single", "cleanMetPt",                {dsDY_nunu_SB0b_Wgt2b,                             dsDY_nunu_SB2b});
    Plotter::DataCollection dcMC_nunu_Wgt2b_ht( "single", "HTZinv",                    {dsDY_nunu_SB0b_Wgt2b,                             dsDY_nunu_SB2b});
    Plotter::DataCollection dcMC_nunu_Wgt2b_nb( "single", "cntCSVSZinv",               {dsDY_nunu_SB0b_Wgt2b,                             dsDY_nunu_SB2b});
    Plotter::DataCollection dcMC_nunu_Wgt2b_nj( "single", "cntNJetsPt30Eta24Zinv",     {dsDY_nunu_SB0b_Wgt2b,                             dsDY_nunu_SB2b});
    Plotter::DataCollection dcMC_nunu_Wgt2b_nt( "single", {{ "nTopCandSortedCntZinv2b", dsDY_nunu_SB0b_Wgt2b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB2b}});
    Plotter::DataCollection dcMC_nunu_Wgt2b_mt2("single", {{"best_had_brJet_MT2Zinv2b", dsDY_nunu_SB0b_Wgt2b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB2b}});
    //3 fake b
    Plotter::DataCollection dcMC_nunu_Wgt3b_met("single", "cleanMetPt",                {dsDY_nunu_SB0b_Wgt3b,                             dsDY_nunu_SB3b});
    Plotter::DataCollection dcMC_nunu_Wgt3b_ht( "single", "HTZinv",                    {dsDY_nunu_SB0b_Wgt3b,                             dsDY_nunu_SB3b});
    Plotter::DataCollection dcMC_nunu_Wgt3b_nb( "single", "cntCSVSZinv",               {dsDY_nunu_SB0b_Wgt3b,                             dsDY_nunu_SB3b});
    Plotter::DataCollection dcMC_nunu_Wgt3b_nj( "single", "cntNJetsPt30Eta24Zinv",     {dsDY_nunu_SB0b_Wgt3b,                             dsDY_nunu_SB3b});
    Plotter::DataCollection dcMC_nunu_Wgt3b_nt( "single", {{ "nTopCandSortedCntZinv3b", dsDY_nunu_SB0b_Wgt3b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB3b}});
    Plotter::DataCollection dcMC_nunu_Wgt3b_mt2("single", {{"best_had_brJet_MT2Zinv3b", dsDY_nunu_SB0b_Wgt3b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB3b}});

    Plotter::DataCollection njetw_nSearchBin(  "single", {{"nSearchBin",                      dsDY_nunu_njet}, {"nSearchBin",             dsDY_nunu}  });

    vh.push_back(PHS("NJetWgt_nSearchBin",            {njetw_nSearchBin}, {2, 1}, "passBaselineZinv",   45,  0,     45,   false, false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("NJetWgt_nSearchBin_log",        {njetw_nSearchBin}, {2, 1}, "passBaselineZinv",   45,  0,     45,   true,  false,  "Search Bin",     "Events", true));
    vh.push_back(PHS("NJetWgt_nSearchBin_pull",       {njetw_nSearchBin}, {2, 1}, "passBaselineZinv",   45,  0,     45,   false, false,  "Search Bin",     "Events", false));
    vh.push_back(PHS("NJetWgt_nSearchBin_pull_log",   {njetw_nSearchBin}, {2, 1}, "passBaselineZinv",   45,  0,     45,   true,  false,  "Search Bin",     "Events", false));

    //Baseline cuts are flawed here, fix!
    vh.push_back(PHS( "ClosureNb_nj_nw_1fakeb_baseline",   {dcMC_nunu_nj_1b},     {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_2fakeb_baseline",   {dcMC_nunu_nj_2b},     {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_3fakeb_baseline",   {dcMC_nunu_nj_3b},     {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_1fakeb_baseline",   {dcMC_nunu_Wgt1b_met}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_1fakeb_baseline",   {dcMC_nunu_Wgt1b_ht},  {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_1fakeb_baseline",   {dcMC_nunu_Wgt1b_nt},  {1, 2}, s_znunu_baseline,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_1fakeb_baseline",   {dcMC_nunu_Wgt1b_nb},  {1, 2}, s_znunu_baseline,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_1fakeb_baseline",   {dcMC_nunu_Wgt1b_mt2}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_1fakeb_baseline",   {dcMC_nunu_Wgt1b_nj},  {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_2fakeb_baseline",   {dcMC_nunu_Wgt2b_met}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_2fakeb_baseline",   {dcMC_nunu_Wgt2b_ht},  {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_2fakeb_baseline",   {dcMC_nunu_Wgt2b_nt},  {1, 2}, s_znunu_baseline,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_2fakeb_baseline",   {dcMC_nunu_Wgt2b_nb},  {1, 2}, s_znunu_baseline,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_2fakeb_baseline",   {dcMC_nunu_Wgt2b_mt2}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_2fakeb_baseline",   {dcMC_nunu_Wgt2b_nj},  {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_3fakeb_baseline",   {dcMC_nunu_Wgt3b_met}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_3fakeb_baseline",   {dcMC_nunu_Wgt3b_ht},  {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_3fakeb_baseline",   {dcMC_nunu_Wgt3b_nt},  {1, 2}, s_znunu_baseline,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_3fakeb_baseline",   {dcMC_nunu_Wgt3b_nb},  {1, 2}, s_znunu_baseline,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_3fakeb_baseline",   {dcMC_nunu_Wgt3b_mt2}, {1, 2}, s_znunu_baseline,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_3fakeb_baseline",   {dcMC_nunu_Wgt3b_nj},  {1, 2}, s_znunu_baseline,  20, 0,   20,   true, false,  "Nj",   ""));

    vh.push_back(PHS( "ClosureNb_nj_nw_1fakeb_loose0",   {dcMC_nunu_nj_1b},     {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_2fakeb_loose0",   {dcMC_nunu_nj_2b},     {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_3fakeb_loose0",   {dcMC_nunu_nj_3b},     {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_1fakeb_loose0",   {dcMC_nunu_Wgt1b_met}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_1fakeb_loose0",   {dcMC_nunu_Wgt1b_ht},  {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_1fakeb_loose0",   {dcMC_nunu_Wgt1b_nt},  {1, 2}, s_znunu_loose0,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_1fakeb_loose0",   {dcMC_nunu_Wgt1b_nb},  {1, 2}, s_znunu_loose0,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_1fakeb_loose0",   {dcMC_nunu_Wgt1b_mt2}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_1fakeb_loose0",   {dcMC_nunu_Wgt1b_nj},  {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_2fakeb_loose0",   {dcMC_nunu_Wgt2b_met}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_2fakeb_loose0",   {dcMC_nunu_Wgt2b_ht},  {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_2fakeb_loose0",   {dcMC_nunu_Wgt2b_nt},  {1, 2}, s_znunu_loose0,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_2fakeb_loose0",   {dcMC_nunu_Wgt2b_nb},  {1, 2}, s_znunu_loose0,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_2fakeb_loose0",   {dcMC_nunu_Wgt2b_mt2}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_2fakeb_loose0",   {dcMC_nunu_Wgt2b_nj},  {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_3fakeb_loose0",   {dcMC_nunu_Wgt3b_met}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_3fakeb_loose0",   {dcMC_nunu_Wgt3b_ht},  {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_3fakeb_loose0",   {dcMC_nunu_Wgt3b_nt},  {1, 2}, s_znunu_loose0,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_3fakeb_loose0",   {dcMC_nunu_Wgt3b_nb},  {1, 2}, s_znunu_loose0,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_3fakeb_loose0",   {dcMC_nunu_Wgt3b_mt2}, {1, 2}, s_znunu_loose0,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_3fakeb_loose0",   {dcMC_nunu_Wgt3b_nj},  {1, 2}, s_znunu_loose0,  20, 0,   20,   true, false,  "Nj",   ""));

    vh.push_back(PHS( "ClosureNb_nj_nw_1fakeb_loose50",   {dcMC_nunu_nj_1b},     {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_2fakeb_loose50",   {dcMC_nunu_nj_2b},     {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_3fakeb_loose50",   {dcMC_nunu_nj_3b},     {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_1fakeb_loose50",   {dcMC_nunu_Wgt1b_met}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_1fakeb_loose50",   {dcMC_nunu_Wgt1b_ht},  {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_1fakeb_loose50",   {dcMC_nunu_Wgt1b_nt},  {1, 2}, s_znunu_loose50,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_1fakeb_loose50",   {dcMC_nunu_Wgt1b_nb},  {1, 2}, s_znunu_loose50,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_1fakeb_loose50",   {dcMC_nunu_Wgt1b_mt2}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_1fakeb_loose50",   {dcMC_nunu_Wgt1b_nj},  {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_2fakeb_loose50",   {dcMC_nunu_Wgt2b_met}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_2fakeb_loose50",   {dcMC_nunu_Wgt2b_ht},  {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_2fakeb_loose50",   {dcMC_nunu_Wgt2b_nt},  {1, 2}, s_znunu_loose50,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_2fakeb_loose50",   {dcMC_nunu_Wgt2b_nb},  {1, 2}, s_znunu_loose50,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_2fakeb_loose50",   {dcMC_nunu_Wgt2b_mt2}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_2fakeb_loose50",   {dcMC_nunu_Wgt2b_nj},  {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_3fakeb_loose50",   {dcMC_nunu_Wgt3b_met}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_3fakeb_loose50",   {dcMC_nunu_Wgt3b_ht},  {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_3fakeb_loose50",   {dcMC_nunu_Wgt3b_nt},  {1, 2}, s_znunu_loose50,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_3fakeb_loose50",   {dcMC_nunu_Wgt3b_nb},  {1, 2}, s_znunu_loose50,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_3fakeb_loose50",   {dcMC_nunu_Wgt3b_mt2}, {1, 2}, s_znunu_loose50,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_3fakeb_loose50",   {dcMC_nunu_Wgt3b_nj},  {1, 2}, s_znunu_loose50,  20, 0,   20,   true, false,  "Nj",   ""));

    vh.push_back(PHS( "ClosureNb_nj_nw_1fakeb_loose100",   {dcMC_nunu_nj_1b},     {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_2fakeb_loose100",   {dcMC_nunu_nj_2b},     {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_3fakeb_loose100",   {dcMC_nunu_nj_3b},     {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_1fakeb_loose100",   {dcMC_nunu_Wgt1b_met}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_1fakeb_loose100",   {dcMC_nunu_Wgt1b_ht},  {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_1fakeb_loose100",   {dcMC_nunu_Wgt1b_nt},  {1, 2}, s_znunu_loose100,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_1fakeb_loose100",   {dcMC_nunu_Wgt1b_nb},  {1, 2}, s_znunu_loose100,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_1fakeb_loose100",   {dcMC_nunu_Wgt1b_mt2}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_1fakeb_loose100",   {dcMC_nunu_Wgt1b_nj},  {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_2fakeb_loose100",   {dcMC_nunu_Wgt2b_met}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_2fakeb_loose100",   {dcMC_nunu_Wgt2b_ht},  {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_2fakeb_loose100",   {dcMC_nunu_Wgt2b_nt},  {1, 2}, s_znunu_loose100,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_2fakeb_loose100",   {dcMC_nunu_Wgt2b_nb},  {1, 2}, s_znunu_loose100,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_2fakeb_loose100",   {dcMC_nunu_Wgt2b_mt2}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_2fakeb_loose100",   {dcMC_nunu_Wgt2b_nj},  {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_3fakeb_loose100",   {dcMC_nunu_Wgt3b_met}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_3fakeb_loose100",   {dcMC_nunu_Wgt3b_ht},  {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_3fakeb_loose100",   {dcMC_nunu_Wgt3b_nt},  {1, 2}, s_znunu_loose100,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_3fakeb_loose100",   {dcMC_nunu_Wgt3b_nb},  {1, 2}, s_znunu_loose100,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_3fakeb_loose100",   {dcMC_nunu_Wgt3b_mt2}, {1, 2}, s_znunu_loose100,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_3fakeb_loose100",   {dcMC_nunu_Wgt3b_nj},  {1, 2}, s_znunu_loose100,  20, 0,   20,   true, false,  "Nj",   ""));

    vh.push_back(PHS( "ClosureNb_nj_nw_1fakeb_loose200",   {dcMC_nunu_nj_1b},     {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_2fakeb_loose200",   {dcMC_nunu_nj_2b},     {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS( "ClosureNb_nj_nw_3fakeb_loose200",   {dcMC_nunu_nj_3b},     {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_1fakeb_loose200",   {dcMC_nunu_Wgt1b_met}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_1fakeb_loose200",   {dcMC_nunu_Wgt1b_ht},  {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_1fakeb_loose200",   {dcMC_nunu_Wgt1b_nt},  {1, 2}, s_znunu_loose200,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_1fakeb_loose200",   {dcMC_nunu_Wgt1b_nb},  {1, 2}, s_znunu_loose200,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_1fakeb_loose200",   {dcMC_nunu_Wgt1b_mt2}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_1fakeb_loose200",   {dcMC_nunu_Wgt1b_nj},  {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_2fakeb_loose200",   {dcMC_nunu_Wgt2b_met}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_2fakeb_loose200",   {dcMC_nunu_Wgt2b_ht},  {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_2fakeb_loose200",   {dcMC_nunu_Wgt2b_nt},  {1, 2}, s_znunu_loose200,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_2fakeb_loose200",   {dcMC_nunu_Wgt2b_nb},  {1, 2}, s_znunu_loose200,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_2fakeb_loose200",   {dcMC_nunu_Wgt2b_mt2}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_2fakeb_loose200",   {dcMC_nunu_Wgt2b_nj},  {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    vh.push_back(PHS(   "ClosureNb_met_3fakeb_loose200",   {dcMC_nunu_Wgt3b_met}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS(    "ClosureNb_ht_3fakeb_loose200",   {dcMC_nunu_Wgt3b_ht},  {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "ht",   ""));
    vh.push_back(PHS(    "ClosureNb_nt_3fakeb_loose200",   {dcMC_nunu_Wgt3b_nt},  {1, 2}, s_znunu_loose200,   5, 0,    5,   true, false,  "Ntop", ""));
    vh.push_back(PHS(    "ClosureNb_nb_3fakeb_loose200",   {dcMC_nunu_Wgt3b_nb},  {1, 2}, s_znunu_loose200,  10, 0,   10,   true, false,  "Nb",   ""));
    vh.push_back(PHS(   "ClosureNb_mt2_3fakeb_loose200",   {dcMC_nunu_Wgt3b_mt2}, {1, 2}, s_znunu_loose200,  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS(    "ClosureNb_nj_3fakeb_loose200",   {dcMC_nunu_Wgt3b_nj},  {1, 2}, s_znunu_loose200,  20, 0,   20,   true, false,  "Nj",   ""));
    //MC interlude over

    // Define all the histograms
    // met
    vh.push_back(PHS("DataMC_SingleMuon_met_nosel",                  {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_2mu",                    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv",                 {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_ht200",           {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_loose0",          {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_loose50",         {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_loose100",        {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_loose200",        {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_blnotagmt2",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_bl",              {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b",              {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_ht200",        {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_loose0",       {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_loose0,       50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_loose50",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_loose50,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_loose100",     {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_loose100,     50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_loose200",     {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_loose200,     50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_0b_blnotag",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_g1b_loose0",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_g1b_loose0,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_g1b_loose50",     {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_g1b_loose50,     50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_g1b_loose100",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_g1b_loose100,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_muZinv_g1b_loose200",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_muZinv_g1b_loose200,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_2el",                      {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv",                   {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_ht200",             {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_blnotagmt2",        {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_bl",                {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_0b",                {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_0b_ht200",          {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_met_elZinv_0b_blnotag",        {dcData_DoubleEG_met, dcMC_met},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmu",                   {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv",               {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_ht200",         {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_loose0",        {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_loose0,        50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_loose50",       {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_loose50,       50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_loose100",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_loose100,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_loose200",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_loose200,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_blnotagmt2",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_bl",            {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b",            {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_ht200",      {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_loose0",     {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_loose0,     50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_loose50",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_loose50,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_loose100",   {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_loose100,   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_loose200",   {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_loose200,   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_0b_blnotag",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_g1b_loose0",    {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_g1b_loose0,    50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_g1b_loose50",   {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_g1b_loose50,   50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_g1b_loose100",  {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_g1b_loose100,  50, 0, 1500,   true, false,  "met",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_met_elmuZinv_g1b_loose200",  {dcData_SingleMuon_met, dcMC_met}, {1, 2}, s_elmuZinv_g1b_loose200,  50, 0, 1500,   true, false,  "met",  ""));

    // ht
    vh.push_back(PHS("DataMC_SingleMuon_ht_nosel",                  {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_2mu",                    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv",                 {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_ht200",           {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_loose0",          {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_loose50",         {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_loose100",        {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_loose200",        {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_blnotagmt2",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_bl",              {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b",              {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_ht200",        {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_loose0",       {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_loose0,       50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_loose50",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_loose50,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_loose100",     {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_loose100,     50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_loose200",     {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_loose200,     50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_0b_blnotag",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_g1b_loose0",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_g1b_loose0,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_g1b_loose50",     {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_g1b_loose50,     50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_g1b_loose100",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_g1b_loose100,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_muZinv_g1b_loose200",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_muZinv_g1b_loose200,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_2el",                      {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv",                   {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_ht200",             {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_blnotagmt2",        {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_bl",                {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_0b",                {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_0b_ht200",          {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_elZinv_0b_blnotag",        {dcData_DoubleEG_ht, dcMC_ht},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmu",                   {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv",               {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_ht200",         {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_loose0",        {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_loose0,        50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_loose50",       {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_loose50,       50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_loose100",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_loose100,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_loose200",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_loose200,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_blnotagmt2",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_bl",            {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b",            {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_ht200",      {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_loose0",     {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_loose0,     50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_loose50",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_loose50,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_loose100",   {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_loose100,   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_loose200",   {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_loose200,   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_0b_blnotag",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_g1b_loose0",    {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_g1b_loose0,    50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_g1b_loose50",   {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_g1b_loose50,   50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_g1b_loose100",  {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_g1b_loose100,  50, 0, 1500,   true, false,  "ht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_ht_elmuZinv_g1b_loose200",  {dcData_SingleMuon_ht, dcMC_ht}, {1, 2}, s_elmuZinv_g1b_loose200,  50, 0, 1500,   true, false,  "ht",  ""));

    // mht
    vh.push_back(PHS("DataMC_SingleMuon_mht_nosel",                  {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_2mu",                    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv",                 {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_ht200",           {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_loose0",          {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_loose50",         {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_loose100",        {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_loose200",        {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_loose200,       50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_blnotagmt2",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_bl",              {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b",              {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_ht200",        {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_loose0",       {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_loose0,       50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_loose50",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_loose50,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_loose100",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_loose100,     50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_loose200",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_loose200,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_0b_blnotag",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_g1b_loose0",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_g1b_loose0,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_g1b_loose50",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_g1b_loose50,     50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_g1b_loose100",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_g1b_loose100,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_muZinv_g1b_loose200",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_muZinv_g1b_loose200,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_2el",                      {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv",                   {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_ht200",             {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_blnotagmt2",        {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_bl",                {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_0b",                {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_0b_ht200",          {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_elZinv_0b_blnotag",        {dcData_DoubleEG_mht, dcMC_mht},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmu",                   {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv",               {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_ht200",         {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_loose0",          {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_loose50",         {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_loose100",        {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_loose200",        {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_loose200,       50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_bl",            {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b",            {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_ht200",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_loose0",       {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_loose0,       50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_loose50",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_loose50,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_loose100",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_loose100,     50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_loose200",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_loose200,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_g1b_loose0",      {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_g1b_loose0,      50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_g1b_loose50",     {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_g1b_loose50,     50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_g1b_loose100",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_g1b_loose100,    50, 0, 1500,   true, false,  "mht",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mht_elmuZinv_g1b_loose200",    {dcData_SingleMuon_mht, dcMC_mht}, {1, 2}, s_elmuZinv_g1b_loose200,    50, 0, 1500,   true, false,  "mht",  ""));

    // nt
    vh.push_back(PHS("DataMC_SingleMuon_nt_nosel",                  {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_nosel,                  5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_2mu",                    {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_2mu,                    5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv",                 {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv,                 5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_ht200",           {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_ht200,           5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_blnotagmt2",      {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_bl",              {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_bl,              5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b",              {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b,              5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_ht200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_0b_blnotag",      {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_g1b_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_g1b_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_g1b_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_muZinv_g1b_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_muZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_2el",                      {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_2el,                    5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv",                   {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv,                 5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_ht200",             {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_ht200,           5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_blnotagmt2",        {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_bl",                {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_bl,              5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_0b",                {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_0b,              5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_0b_ht200",          {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_elZinv_0b_blnotag",        {dcData_DoubleEG_nt, dcMC_nt},   {1, 2}, s_elZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmu",                   {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmu,                   5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv",               {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv,               5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_ht200",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_ht200,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_blnotagmt2,    5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_bl",            {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_bl,            5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b",            {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b,            5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_ht200",      {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_ht200,      5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_0b_blnotag,    5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nt, dcMC_nt}, {1, 2}, s_elmuZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop",  ""));

    // nt1b
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_nosel",                  {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_nosel,                  5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_2mu",                    {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_2mu,                    5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv",                 {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv,                 5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_ht200",           {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_ht200,           5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_blnotagmt2",      {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_bl",              {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_bl,              5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b",              {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b,              5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_ht200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_0b_blnotag",      {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_g1b_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_g1b_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_g1b_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_muZinv_g1b_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_muZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_2el",                      {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_2el,                    5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv",                   {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv,                 5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_ht200",             {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_ht200,           5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_blnotagmt2",        {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_bl",                {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_bl,              5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_0b",                {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_0b,              5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_0b_ht200",          {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_elZinv_0b_blnotag",        {dcData_DoubleEG_nt1b, dcMC_nt1b},   {1, 2}, s_elZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmu",                   {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmu,                   5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv",               {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv,               5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_ht200",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_ht200,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_blnotagmt2,    5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_bl",            {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_bl,            5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b",            {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b,            5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_ht200",      {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_ht200,      5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_0b_blnotag,    5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt1b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nt1b, dcMC_nt1b}, {1, 2}, s_elmuZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (1b fake)",  ""));

    // nt2b
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_nosel",                  {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_nosel,                  5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_2mu",                    {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_2mu,                    5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv",                 {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv,                 5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_ht200",           {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_ht200,           5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_blnotagmt2",      {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_bl",              {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_bl,              5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b",              {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b,              5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_ht200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_0b_blnotag",      {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_g1b_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_g1b_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_g1b_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_muZinv_g1b_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_muZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_2el",                      {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_2el,                    5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv",                   {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv,                 5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_ht200",             {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_ht200,           5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_blnotagmt2",        {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_bl",                {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_bl,              5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_0b",                {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_0b,              5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_0b_ht200",          {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_elZinv_0b_blnotag",        {dcData_DoubleEG_nt2b, dcMC_nt2b},   {1, 2}, s_elZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmu",                   {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmu,                   5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv",               {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv,               5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_ht200",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_ht200,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_blnotagmt2,    5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_bl",            {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_bl,            5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b",            {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b,            5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_ht200",      {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_ht200,      5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_0b_blnotag,    5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt2b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nt2b, dcMC_nt2b}, {1, 2}, s_elmuZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (2b fake)",  ""));

    // nt3b
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_nosel",                  {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_nosel,                  5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_2mu",                    {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_2mu,                    5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv",                 {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv,                 5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_ht200",           {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_ht200,           5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_blnotagmt2",      {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_bl",              {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_bl,              5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b",              {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b,              5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_ht200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_0b_blnotag",      {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_g1b_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_g1b_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_g1b_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_muZinv_g1b_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_muZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_2el",                      {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_2el,                    5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv",                   {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv,                 5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_ht200",             {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_ht200,           5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_blnotagmt2",        {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_blnotagmt2,      5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_bl",                {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_bl,              5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_0b",                {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_0b,              5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_0b_ht200",          {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_0b_ht200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_0b_blnotagmt2,   5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_elZinv_0b_blnotag",        {dcData_DoubleEG_nt3b, dcMC_nt3b},   {1, 2}, s_elZinv_0b_blnotag,      5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmu",                   {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmu,                   5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv",               {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv,               5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_ht200",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_ht200,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_blnotagmt2,    5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_bl",            {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_bl,            5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b",            {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b,            5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_ht200",      {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_ht200,      5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_0b_blnotag,    5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_g1b_loose0,          5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_g1b_loose50,         5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_g1b_loose100,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nt3b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nt3b, dcMC_nt3b}, {1, 2}, s_elmuZinv_g1b_loose200,        5, 0, 5,   true, false,  "Ntop (3b fake)",  ""));

    // mt2
    vh.push_back(PHS("DataMC_SingleMuon_mt2_nosel",                  {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_2mu",                    {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv",                 {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_ht200",           {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_blnotagmt2",      {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_bl",              {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b",              {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_ht200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_0b_blnotag",      {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_g1b_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_g1b_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_g1b_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_muZinv_g1b_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_2el",                      {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv",                   {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_ht200",             {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_blnotagmt2",        {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_bl",                {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_0b",                {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_0b_ht200",          {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_elZinv_0b_blnotag",        {dcData_DoubleEG_mt2, dcMC_mt2},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmu",                   {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv",               {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_ht200",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_bl",            {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b",            {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_ht200",      {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "mt2",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt2_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mt2, dcMC_mt2}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2",  ""));

    // mt21b
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_nosel",                  {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_2mu",                    {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv",                 {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_ht200",           {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_blnotagmt2",      {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_bl",              {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b",              {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_ht200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_0b_blnotag",      {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_g1b_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_g1b_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_g1b_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_muZinv_g1b_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_2el",                      {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv",                   {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_ht200",             {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_blnotagmt2",        {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_bl",                {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_0b",                {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_0b_ht200",          {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_elZinv_0b_blnotag",        {dcData_DoubleEG_mt21b, dcMC_mt21b},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmu",                   {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv",               {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_ht200",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_bl",            {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b",            {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_ht200",      {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt21b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mt21b, dcMC_mt21b}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2 (1b fake)",  ""));

    // mt22b
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_nosel",                  {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_2mu",                    {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv",                 {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_ht200",           {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_blnotagmt2",      {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_bl",              {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b",              {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_ht200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_0b_blnotag",      {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_g1b_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_g1b_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_g1b_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_g1b_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_g1b_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_muZinv_g1b_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_muZinv_g1b_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_2el",                      {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv",                   {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_ht200",             {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_blnotagmt2",        {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_bl",                {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_0b",                {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_0b_ht200",          {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_elZinv_0b_blnotag",        {dcData_DoubleEG_mt22b, dcMC_mt22b},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmu",                   {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv",               {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_ht200",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_bl",            {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b",            {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_ht200",      {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_g1b_loose50,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_g1b_loose100,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt22b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mt22b, dcMC_mt22b}, {1, 2}, s_elmuZinv_g1b_loose200,          50, 0, 1500,   true, false,  "mt2 (2b fake)",  ""));

    // mt23b
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_nosel",                  {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_2mu",                    {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv",                 {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_ht200",           {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_g1b_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_g1b_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_g1b_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_g1b_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_g1b_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_blnotagmt2",      {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_bl",              {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b",              {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_ht200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_muZinv_0b_blnotag",      {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_2el",                      {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv",                   {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_ht200",             {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_blnotagmt2",        {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_bl",                {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_0b",                {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_0b_ht200",          {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_elZinv_0b_blnotag",        {dcData_DoubleEG_mt23b, dcMC_mt23b},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmu",                   {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv",               {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_ht200",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_g1b_loose100,         50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_bl",            {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b",            {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_ht200",      {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mt23b_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mt23b, dcMC_mt23b}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "mt2 (3b fake)",  ""));

    // nb
    vh.push_back(PHS("DataMC_SingleMuon_nb_nosel",                  {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_nosel,                  10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_2mu",                    {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_2mu,                    10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv",                 {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv,                 10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_ht200",           {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_ht200,           10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_loose200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_blnotagmt2",      {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_blnotagmt2,      10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_bl",              {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_bl,              10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b",              {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b,              10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_ht200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_ht200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_loose200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_blnotagmt2,   10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_0b_blnotag",      {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_0b_blnotag,      10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_g1b_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_g1b_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_g1b_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_g1b_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_g1b_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_g1b_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_muZinv_g1b_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_muZinv_g1b_loose200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_2el",                      {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_2el,                    10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv",                   {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv,                 10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_ht200",             {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_ht200,           10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_blnotagmt2",        {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_blnotagmt2,      10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_bl",                {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_bl,              10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_0b",                {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_0b,              10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_0b_ht200",          {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_0b_ht200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_0b_blnotagmt2,   10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_elZinv_0b_blnotag",        {dcData_DoubleEG_nb, dcMC_nb},   {1, 2}, s_elZinv_0b_blnotag,      10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmu",                   {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmu,                   10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv",               {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv,               10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_ht200",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_ht200,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_loose200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_blnotagmt2,    10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_bl",            {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_bl,            10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b",            {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b,            10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_ht200",      {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_ht200,      10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_loose200,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_0b_blnotag,    10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_g1b_loose0,          10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_g1b_loose50,         10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_g1b_loose100,        10, 0, 10,   true, false,  "Nb",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nb_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nb, dcMC_nb}, {1, 2}, s_elmuZinv_g1b_loose200,        10, 0, 10,   true, false,  "Nb",  ""));

    // nj
    vh.push_back(PHS("DataMC_SingleMuon_nj_nosel",                  {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_nosel,                  20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_2mu",                    {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_2mu,                    20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv",                 {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv,                 20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_ht200",           {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_ht200,           20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_loose200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_blnotagmt2",      {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_blnotagmt2,      20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_bl",              {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_bl,              20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b",              {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b,              20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_ht200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_ht200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_loose200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_blnotagmt2,   20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_0b_blnotag",      {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_0b_blnotag,      20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_g1b_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_g1b_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_g1b_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_g1b_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_g1b_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_g1b_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_muZinv_g1b_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_muZinv_g1b_loose200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_2el",                      {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_2el,                    20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv",                   {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv,                 20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_ht200",             {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_ht200,           20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_blnotagmt2",        {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_blnotagmt2,      20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_bl",                {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_bl,              20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_0b",                {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_0b,              20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_0b_ht200",          {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_0b_ht200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_0b_blnotagmt2,   20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_elZinv_0b_blnotag",        {dcData_DoubleEG_nj, dcMC_nj},   {1, 2}, s_elZinv_0b_blnotag,      20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmu",                   {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmu,                   20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv",               {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv,               20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_ht200",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_ht200,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_loose200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_blnotagmt2",    {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_blnotagmt2,    20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_bl",            {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_bl,            20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b",            {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b,            20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_ht200",      {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_ht200,      20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_loose200,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_0b_blnotag",    {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_0b_blnotag,    20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_g1b_loose0",          {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_g1b_loose0,          20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_g1b_loose50",         {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_g1b_loose50,         20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_g1b_loose100",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_g1b_loose100,        20, 0, 20,   true, false,  "Nj",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_nj_elmuZinv_g1b_loose200",        {dcData_SingleMuon_nj, dcMC_nj}, {1, 2}, s_elmuZinv_g1b_loose200,        20, 0, 20,   true, false,  "Nj",  ""));

    // jpt
    vh.push_back(PHS("DataMC_SingleMuon_jpt_nosel",                  {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_2mu",                    {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv",                 {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_ht200",           {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_blnotagmt2",      {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_bl",              {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b",              {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_ht200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_0b_blnotag",      {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_g1b_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_g1b_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_g1b_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_muZinv_g1b_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_2el",                      {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv",                   {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_ht200",             {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_blnotagmt2",        {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_bl",                {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_0b",                {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_0b_ht200",          {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_elZinv_0b_blnotag",        {dcData_DoubleEG_jpt, dcMC_jpt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmu",                   {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv",               {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_ht200",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_bl",            {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b",            {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_ht200",      {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_jpt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_jpt, dcMC_jpt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet pt",  ""));

    // j1pt
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_nosel",                  {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_2mu",                    {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv",                 {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_ht200",           {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_g1b_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_g1b_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_g1b_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_g1b_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_blnotagmt2",      {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_bl",              {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b",              {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_ht200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_muZinv_0b_blnotag",      {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_2el",                      {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv",                   {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_ht200",             {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_blnotagmt2",        {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_bl",                {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_0b",                {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_0b_ht200",          {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_elZinv_0b_blnotag",        {dcData_DoubleEG_j1pt, dcMC_j1pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmu",                   {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv",               {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_ht200",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_bl",            {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b",            {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "jet1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j1pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_j1pt, dcMC_j1pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "jet1 pt",  ""));

    // j2pt
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_nosel",                  {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_2mu",                    {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv",                 {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_ht200",           {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_g1b_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_g1b_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_g1b_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_g1b_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_blnotagmt2",      {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_bl",              {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b",              {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_ht200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_muZinv_0b_blnotag",      {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_2el",                      {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv",                   {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_ht200",             {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_blnotagmt2",        {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_bl",                {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_0b",                {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_0b_ht200",          {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_elZinv_0b_blnotag",        {dcData_DoubleEG_j2pt, dcMC_j2pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmu",                   {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv",               {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_ht200",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_bl",            {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b",            {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "jet2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j2pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_j2pt, dcMC_j2pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "jet2 pt",  ""));

    // j3pt
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_nosel",                  {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_nosel,                  50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_2mu",                    {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_2mu,                    50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv",                 {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv,                 50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_ht200",           {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_g1b_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_g1b_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_g1b_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_g1b_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_blnotagmt2",      {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_bl",              {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_bl,              50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b",              {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b,              50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_ht200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_muZinv_0b_blnotag",      {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_2el",                      {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_2el,                    50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv",                   {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv,                 50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_ht200",             {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_blnotagmt2",        {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_bl",                {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_bl,              50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_0b",                {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_0b,              50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_0b_ht200",          {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_elZinv_0b_blnotag",        {dcData_DoubleEG_j3pt, dcMC_j3pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmu",                   {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmu,                   50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv",               {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv,               50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_ht200",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_bl",            {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b",            {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1500,   true, false,  "jet3 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_j3pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_j3pt, dcMC_j3pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1500,   true, false,  "jet3 pt",  ""));

    // mupt
    vh.push_back(PHS("DataMC_SingleMuon_mupt_nosel",                  {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_2mu",                    {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv",                 {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_ht200",           {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_g1b_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_g1b_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_g1b_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_g1b_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_blnotagmt2",      {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_bl",              {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b",              {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_ht200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_muZinv_0b_blnotag",      {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_2el",                      {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv",                   {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_ht200",             {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_blnotagmt2",        {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_bl",                {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_0b",                {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_0b_ht200",          {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_elZinv_0b_blnotag",        {dcData_DoubleEG_mupt, dcMC_mupt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmu",                   {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv",               {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_ht200",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_bl",            {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b",            {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_ht200",      {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "mu pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mupt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mupt, dcMC_mupt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "mu pt",  ""));

    // mu1pt
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_nosel",                  {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_2mu",                    {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv",                 {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_ht200",           {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_g1b_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_g1b_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_g1b_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_g1b_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_blnotagmt2",      {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_bl",              {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b",              {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_ht200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_muZinv_0b_blnotag",      {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_2el",                      {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv",                   {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_ht200",             {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_blnotagmt2",        {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_bl",                {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_0b",                {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_0b_ht200",          {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_elZinv_0b_blnotag",        {dcData_DoubleEG_mu1pt, dcMC_mu1pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmu",                   {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv",               {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_ht200",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_bl",            {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b",            {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "mu1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu1pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "mu1 pt",  ""));

    // mu2pt
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_nosel",                  {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_2mu",                    {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv",                 {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_ht200",           {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_g1b_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_g1b_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_g1b_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_g1b_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_blnotagmt2",      {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_bl",              {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b",              {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_ht200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_muZinv_0b_blnotag",      {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_2el",                      {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv",                   {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_ht200",             {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_blnotagmt2",        {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_bl",                {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_0b",                {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_0b_ht200",          {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_elZinv_0b_blnotag",        {dcData_DoubleEG_mu2pt, dcMC_mu2pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmu",                   {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv",               {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_ht200",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_bl",            {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b",            {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "mu2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mu2pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "mu2 pt",  ""));

    // elpt
    vh.push_back(PHS("DataMC_SingleMuon_elpt_nosel",                  {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_2mu",                    {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv",                 {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_ht200",           {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_g1b_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_g1b_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_g1b_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_g1b_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_blnotagmt2",      {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_bl",              {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b",              {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_ht200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_muZinv_0b_blnotag",      {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_2el",                      {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv",                   {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_ht200",             {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_blnotagmt2",        {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_bl",                {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_0b",                {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_0b_ht200",          {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_elZinv_0b_blnotag",        {dcData_DoubleEG_elpt, dcMC_elpt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmu",                   {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv",               {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_ht200",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_bl",            {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b",            {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_ht200",      {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "el pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_elpt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_elpt, dcMC_elpt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "el pt",  ""));

    // el1pt
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_nosel",                  {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_2mu",                    {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv",                 {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_ht200",           {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_g1b_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_g1b_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_g1b_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_g1b_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_blnotagmt2",      {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_bl",              {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b",              {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_ht200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_muZinv_0b_blnotag",      {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_2el",                      {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv",                   {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_ht200",             {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_blnotagmt2",        {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_bl",                {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_0b",                {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_0b_ht200",          {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_elZinv_0b_blnotag",        {dcData_DoubleEG_el1pt, dcMC_el1pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmu",                   {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv",               {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_ht200",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_bl",            {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b",            {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "el1 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el1pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "el1 pt",  ""));

    // el2pt
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_nosel",                  {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_nosel,                  50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_2mu",                    {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_2mu,                    50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv",                 {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv,                 50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_ht200",           {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_ht200,           50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_g1b_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_g1b_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_g1b_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_g1b_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_blnotagmt2",      {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_bl",              {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_bl,              50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b",              {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b,              50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_ht200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_ht200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_muZinv_0b_blnotag",      {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_muZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_2el",                      {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_2el,                    50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv",                   {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv,                 50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_ht200",             {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_ht200,           50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_blnotagmt2",        {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_blnotagmt2,      50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_bl",                {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_bl,              50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_0b",                {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_0b,              50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_0b_ht200",          {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_0b_ht200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_0b_blnotagmt2,   50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_elZinv_0b_blnotag",        {dcData_DoubleEG_el2pt, dcMC_el2pt},   {1, 2}, s_elZinv_0b_blnotag,      50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmu",                   {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmu,                   50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv",               {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv,               50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_ht200",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_ht200,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_g1b_loose0",          {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_g1b_loose0,          50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_g1b_loose50",         {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_g1b_loose50,         50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_g1b_loose100",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_g1b_loose100,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_g1b_loose200",        {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_g1b_loose200,        50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_blnotagmt2",    {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_blnotagmt2,    50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_bl",            {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_bl,            50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b",            {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b,            50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_ht200",      {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_ht200,      50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 50, 0, 1000,   true, false,  "el2 pt",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_el2pt_elmuZinv_0b_blnotag",    {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, s_elmuZinv_0b_blnotag,    50, 0, 1000,   true, false,  "el2 pt",  ""));

    // mll
    vh.push_back(PHS("DataMC_SingleMuon_mll_nosel",                  {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_nosel,                  40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_2mu",                    {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_2mu,                    40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv",                 {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv,                 40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_ht200",           {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_ht200,           40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_g1b_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_g1b_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_g1b_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_g1b_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_g1b_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_g1b_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_g1b_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_g1b_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_blnotagmt2",      {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_blnotagmt2,      40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_bl",              {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_bl,              40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b",              {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b,              40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_ht200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_ht200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_blnotagmt2",   {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_blnotagmt2,   40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_muZinv_0b_blnotag",      {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_muZinv_0b_blnotag,      40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_2el",                      {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_2el,                    40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv",                   {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv,                 40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_ht200",             {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_ht200,           40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_blnotagmt2",        {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_blnotagmt2,      40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_bl",                {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_bl,              40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_0b",                {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_0b,              40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_0b_ht200",          {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_0b_ht200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_0b_blnotagmt2",     {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_0b_blnotagmt2,   40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_elZinv_0b_blnotag",        {dcData_DoubleEG_mll, dcMC_mll},   {1, 2}, s_elZinv_0b_blnotag,      40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmu",                   {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmu,                   40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv",               {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv,               40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_ht200",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_ht200,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_g1b_loose0",          {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_g1b_loose0,          40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_g1b_loose50",         {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_g1b_loose50,         40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_g1b_loose100",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_g1b_loose100,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_g1b_loose200",        {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_g1b_loose200,        40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_blnotagmt2",    {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_blnotagmt2,    40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_bl",            {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_bl,            40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b",            {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b,            40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_ht200",      {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_ht200,      40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_blnotagmt2", {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_blnotagmt2, 40, 0, 200,   true, false,  "mll",  ""));
    vh.push_back(PHS("DataMC_SingleMuon_mll_elmuZinv_0b_blnotag",    {dcData_SingleMuon_mll, dcMC_mll}, {1, 2}, s_elmuZinv_0b_blnotag,    40, 0, 200,   true, false,  "mll",  ""));



    //Generate cutflows
    vector<string> cfsZ = {"",
                           "",
                           "passLeptVeto",
                           "passLeptVeto",
                           "passLeptVeto;passnJetsZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
                           "passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
                           "passLeptVeto;passBaselineZinv"};

    vector<string> cfsDYmm = {"",
                              "",
                              "passMuZinvSel",
                              "passMuZinvSel",
                              "passMuZinvSel;passnJetsZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
                              "passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
                              "passMuZinvSel;passBaselineZinv"};

    vector<string> cfsDatamm = {"",
                                "passNoiseEventFilterZinv;",
                                "passNoiseEventFilterZinv;passMuZinvSel",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
                                "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passBaselineZinv"};

    vector<string> cfsDataem = {"",
                                "passNoiseEventFilterZinv;",
                                "passNoiseEventFilterZinv;passElMuZinvSel",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
                                "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passBaselineZinv"};

    vector<Plotter::CutFlowSummary> cutFlowSummaries;
    
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("ZtoNuNu",           PDC("", "", {dsDY_nunu}),           cfsZ));
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("DYtoMuMu",          PDC("", "", {dsDY_ll_zAcc_scaled}), cfsDYmm));
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("DYtoMuMu_Unscaled", PDC("", "", {dsDY_ll_scaled}),      cfsDYmm));
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Data_MuMu",         PDC("", "", {dsData_SingleMuon}),   cfsDatamm));
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Data_ElMu",         PDC("", "", {dsData_SingleMuon}),   cfsDataem));
    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("TTBar_ElMu",        PDC("", "", {dstt2l}),              cfsDataem));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple();

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setCutFlows(cutFlowSummaries);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave)  plotter.saveHists();
    if(doPlots) plotter.plot();
}

//  LocalWords:  mll
