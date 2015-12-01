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
        //sampleloc = "";
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

    // -------------------------------
    // - Data/MC plots: 2mu, 2e, emu
    // -------------------------------

    // DatasetSummary for each cut level and process
    // --> No selection apart from noise filters -- mainly for basic debugging
    Plotter::DatasetSummary dsData_2015C_nosel("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsDY_nosel(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dstt2l_nosel(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dstW_nosel(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsttZ_nosel(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsVV_nosel(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv", "");

    // 2mu datasetsummary
    // --> Loose region: only presence of 2 muons
    Plotter::DatasetSummary dsData_2015C_2mu("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passDiMuSel", "");
    Plotter::DatasetSummary dsDY_2mu(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passDiMuSel", "");
    Plotter::DatasetSummary dstt2l_2mu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passDiMuSel", "");
    Plotter::DatasetSummary dstW_2mu(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passDiMuSel", "");
    Plotter::DatasetSummary dsttZ_2mu(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passDiMuSel", "");
    Plotter::DatasetSummary dsVV_2mu(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passDiMuSel", "");
    // --> Loose region: only passMuZinvSel
    Plotter::DatasetSummary dsData_2015C_muZinv("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_muZinv(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_muZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstW_muZinv(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_muZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsVV_muZinv(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passMuZinvSel", "");
    // --> Loose region: only passMuZinvSel + HT>200
    Plotter::DatasetSummary dsData_2015C_ht200("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsDY_ht200(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dstt2l_ht200(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dstW_ht200(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsttZ_ht200(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsVV_ht200(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passMuZinvSel;ht>200", "");
    // --> baseline without btag, top tag and mt2 cuts
    Plotter::DatasetSummary dsData_2015C_blnotag("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_blnotag(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_blnotag(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstW_blnotag(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_blnotag(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsVV_blnotag(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    // --> full baseline
    Plotter::DatasetSummary dsData_2015C_bl("Data",       fileMap["Data_SingleMuon"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_bl(        "DY",         fileMap["DYJetsToLL"],      "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_bl(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstW_bl(        "single top", fileMap["tW"],              "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_bl(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsVV_bl(        "Rare",       fileMap["Rare"],            "passBaselineZinv;passMuZinvSel", "");
    // --> 0b + muZinv
    Plotter::DatasetSummary dsData_2015C_0bmuZinv("Data",       fileMap["Data_SingleMuon"], "passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0bmuZinv(        "DY",         fileMap["DYJetsToLL"],      "passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0bmuZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0bmuZinv(        "single top", fileMap["tW"],              "passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0bmuZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0bmuZinv(        "Rare",       fileMap["Rare"],            "passMuZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0bnomt2(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0bnomt2(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0bnomt2(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0bnomt2(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0bnomt2(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0b(        "single top", fileMap["tW"],              "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");

    // Weighted according to 1 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0bmuZinv_w1b("Data",       fileMap["Data_SingleMuon"], "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0bmuZinv_w1b(        "DY",         fileMap["DYJetsToLL"],      "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0bmuZinv_w1b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0bmuZinv_w1b(        "single top", fileMap["tW"],              "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0bmuZinv_w1b(       "t#bar{t}Z",  fileMap["TTZ"],             "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0bmuZinv_w1b(        "Rare",       fileMap["Rare"],            "passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w1b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w1b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w1b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w1b(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w1b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w1b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w1b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0b_w1b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w1b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0b_w1b(        "single top", fileMap["tW"],              "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w1b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0b_w1b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv1b;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");

    // Weighted according to 2 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0bmuZinv_w2b("Data",       fileMap["Data_SingleMuon"], "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0bmuZinv_w2b(        "DY",         fileMap["DYJetsToLL"],      "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0bmuZinv_w2b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0bmuZinv_w2b(        "single top", fileMap["tW"],              "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0bmuZinv_w2b(       "t#bar{t}Z",  fileMap["TTZ"],             "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0bmuZinv_w2b(        "Rare",       fileMap["Rare"],            "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w2b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w2b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w2b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w2b(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w2b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w2b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w2b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0b_w2b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w2b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0b_w2b(        "single top", fileMap["tW"],              "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w2b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0b_w2b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv2b;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");

    // Weighted according to 3 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0bmuZinv_w3b("Data",       fileMap["Data_SingleMuon"], "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0bmuZinv_w3b(        "DY",         fileMap["DYJetsToLL"],      "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0bmuZinv_w3b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0bmuZinv_w3b(        "single top", fileMap["tW"],              "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0bmuZinv_w3b(       "t#bar{t}Z",  fileMap["TTZ"],             "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0bmuZinv_w3b(        "Rare",       fileMap["Rare"],            "passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w3b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w3b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w3b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w3b(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w3b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w3b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w3b("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0b_w3b(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w3b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0b_w3b(        "single top", fileMap["tW"],              "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w3b(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0b_w3b(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv3b;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");


    // 2e datasetsummary
    // --> Loose region: only presence of 2 electrons
    Plotter::DatasetSummary dsData_2015C_2el("Data",       fileMap["Data_DoubleEG"],   "passNoiseEventFilterZinv;passDiElecSel", "");
    Plotter::DatasetSummary dsDY_2el(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passDiElecSel", "");
    Plotter::DatasetSummary dstt2l_2el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passDiElecSel", "");
    Plotter::DatasetSummary dstW_2el(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passDiElecSel", "");
    Plotter::DatasetSummary dsttZ_2el(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passDiElecSel", "");
    Plotter::DatasetSummary dsVV_2el(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passDiElecSel", "");
    // --> Loose region: only passElZinvSel
    Plotter::DatasetSummary dsData_2015C_elZinv("Data",       fileMap["Data_DoubleEG"],   "passNoiseEventFilterZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsDY_elZinv(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstt2l_elZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstW_elZinv(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsttZ_elZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsVV_elZinv(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passElecZinvSel", "");
    // --> Loose region: only passElecZinvSel + HT>200
    Plotter::DatasetSummary dsData_2015C_ht200_el("Data",       fileMap["Data_DoubleEG"],   "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    Plotter::DatasetSummary dsDY_ht200_el(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    Plotter::DatasetSummary dstt2l_ht200_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    Plotter::DatasetSummary dstW_ht200_el(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    Plotter::DatasetSummary dsttZ_ht200_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    Plotter::DatasetSummary dsVV_ht200_el(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passElecZinvSel;ht>200", "");
    // --> baseline without btag, top tag and mt2 cuts
    Plotter::DatasetSummary dsData_2015C_blnotag_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsDY_blnotag_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstt2l_blnotag_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstW_blnotag_el(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsttZ_blnotag_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsVV_blnotag_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passElecZinvSel", "");
    // --> full baseline
    Plotter::DatasetSummary dsData_2015C_bl_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsDY_bl_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstt2l_bl_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dstW_bl_el(        "single top", fileMap["tW"],              "passBaselineZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsttZ_bl_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineZinv;passElecZinvSel", "");
    Plotter::DatasetSummary dsVV_bl_el(        "Rare",       fileMap["Rare"],            "passBaselineZinv;passElecZinvSel", "");
    // --> 0b + muZinv
    Plotter::DatasetSummary dsData_2015C_0belZinv("Data",       fileMap["Data_DoubleEG"],   "passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0belZinv(        "DY",         fileMap["DYJetsToLL"],      "passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0belZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0belZinv(        "single top", fileMap["tW"],              "passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0belZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0belZinv(        "Rare",       fileMap["Rare"],            "passElecZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0bnomt2_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0bnomt2_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0bnomt2_el(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0bnomt2_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0bnomt2_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passElecZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0b_el(        "single top", fileMap["tW"],              "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv;passElecZinvSel;cntCSVSZinv=0", "");

    // Weighted according to 1 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belZinv_w1b("Data",       fileMap["Data_DoubleEG"],   "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0belZinv_w1b(        "DY",         fileMap["DYJetsToLL"],      "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0belZinv_w1b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0belZinv_w1b(        "single top", fileMap["tW"],              "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0belZinv_w1b(       "t#bar{t}Z",  fileMap["TTZ"],             "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0belZinv_w1b(        "Rare",       fileMap["Rare"],            "passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w1b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w1b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w1b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w1b_el(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w1b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w1b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w1b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0b_w1b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w1b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0b_w1b_el(        "single top", fileMap["tW"],              "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w1b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0b_w1b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv1b;passElecZinvSel;cntCSVSZinv=0", "weight1fakebComb");

    // Weighted according to 2 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belZinv_w2b("Data",       fileMap["Data_DoubleEG"],   "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0belZinv_w2b(        "DY",         fileMap["DYJetsToLL"],      "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0belZinv_w2b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0belZinv_w2b(        "single top", fileMap["tW"],              "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0belZinv_w2b(       "t#bar{t}Z",  fileMap["TTZ"],             "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0belZinv_w2b(        "Rare",       fileMap["Rare"],            "passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w2b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w2b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w2b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w2b_el(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w2b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w2b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w2b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0b_w2b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w2b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0b_w2b_el(        "single top", fileMap["tW"],              "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w2b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0b_w2b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv2b;passElecZinvSel;cntCSVSZinv=0", "weight2fakebComb");

    // Weighted according to 3 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belZinv_w3b("Data",       fileMap["Data_DoubleEG"],   "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0belZinv_w3b(        "DY",         fileMap["DYJetsToLL"],      "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0belZinv_w3b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0belZinv_w3b(        "single top", fileMap["tW"],              "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0belZinv_w3b(       "t#bar{t}Z",  fileMap["TTZ"],             "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0belZinv_w3b(        "Rare",       fileMap["Rare"],            "passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w3b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w3b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w3b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w3b_el(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w3b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w3b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w3b_el("Data",       fileMap["Data_DoubleEG"],   "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0b_w3b_el(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w3b_el(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0b_w3b_el(        "single top", fileMap["tW"],              "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w3b_el(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0b_w3b_el(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv3b;passElecZinvSel;cntCSVSZinv=0", "weight3fakebComb");


    // emu datasetsummary for ttbar subtraction
    // --> Loose region: only passElMuSel (no zmass cut)
    Plotter::DatasetSummary dsData_2015C_elmu("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passElMuSel", "");
    Plotter::DatasetSummary dsDY_elmu(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passElMuSel", "");
    Plotter::DatasetSummary dstt2l_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passElMuSel", "");
    Plotter::DatasetSummary dstW_elmu(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passElMuSel", "");
    Plotter::DatasetSummary dsttZ_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passElMuSel", "");
    Plotter::DatasetSummary dsVV_elmu(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passElMuSel", "");
    // --> Loose region: only passElMuZinvSel
    Plotter::DatasetSummary dsData_2015C_elmuZinv("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsDY_elmuZinv(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_elmuZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstW_elmuZinv(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_elmuZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsVV_elmuZinv(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passElMuZinvSel", "");
    // --> Loose region: only passElMuZinvSel + HT>200
    Plotter::DatasetSummary dsData_2015C_ht200_elmu("Data",       fileMap["Data_SingleMuon"], "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsDY_ht200_elmu(        "DY",         fileMap["DYJetsToLL"],      "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dstt2l_ht200_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dstW_ht200_elmu(        "single top", fileMap["tW"],              "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsttZ_ht200_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    Plotter::DatasetSummary dsVV_ht200_elmu(        "Rare",       fileMap["Rare"],            "passNoiseEventFilterZinv;passElMuZinvSel;ht>200", "");
    // --> baseline without btag, top tag and mt2 cuts
    Plotter::DatasetSummary dsData_2015C_blnotag_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsDY_blnotag_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_blnotag_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstW_blnotag_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_blnotag_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsVV_blnotag_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passElMuZinvSel", "");
    // --> full baseline
    Plotter::DatasetSummary dsData_2015C_bl_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsDY_bl_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_bl_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dstW_bl_elmu(        "single top", fileMap["tW"],              "passBaselineZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_bl_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineZinv;passElMuZinvSel", "");
    Plotter::DatasetSummary dsVV_bl_elmu(        "Rare",       fileMap["Rare"],            "passBaselineZinv;passElMuZinvSel", "");
    // --> 0b + muZinv
    Plotter::DatasetSummary dsData_2015C_0belmuZinv("Data",       fileMap["Data_SingleMuon"], "passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0belmuZinv(        "DY",         fileMap["DYJetsToLL"],      "passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0belmuZinv(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0belmuZinv(        "single top", fileMap["tW"],              "passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0belmuZinv(       "t#bar{t}Z",  fileMap["TTZ"],             "passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0belmuZinv(        "Rare",       fileMap["Rare"],            "passElMuZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0bnomt2_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0bnomt2_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0bnomt2_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0bnomt2_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0bnomt2_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv;passElMuZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv;passElMuZinvSel;cntCSVSZinv=0", "");

    // Weighted according to 1 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belmuZinv_w1b("Data",       fileMap["Data_SingleMuon"], "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0belmuZinv_w1b(        "DY",         fileMap["DYJetsToLL"],      "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0belmuZinv_w1b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0belmuZinv_w1b(        "single top", fileMap["tW"],              "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0belmuZinv_w1b(       "t#bar{t}Z",  fileMap["TTZ"],             "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0belmuZinv_w1b(        "Rare",       fileMap["Rare"],            "passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w1b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w1b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w1b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w1b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w1b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w1b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w1b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0b_w1b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w1b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0b_w1b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w1b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0b_w1b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv1b;passElMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");

    // Weighted according to 2 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belmuZinv_w2b("Data",       fileMap["Data_SingleMuon"], "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0belmuZinv_w2b(        "DY",         fileMap["DYJetsToLL"],      "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0belmuZinv_w2b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0belmuZinv_w2b(        "single top", fileMap["tW"],              "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0belmuZinv_w2b(       "t#bar{t}Z",  fileMap["TTZ"],             "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0belmuZinv_w2b(        "Rare",       fileMap["Rare"],            "passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w2b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w2b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w2b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w2b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w2b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w2b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w2b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0b_w2b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w2b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0b_w2b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w2b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0b_w2b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv2b;passElMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");

    // Weighted according to 3 fake b
    // --> 0b control region, muZinv
    Plotter::DatasetSummary dsData_2015C_0belmuZinv_w3b("Data",       fileMap["Data_SingleMuon"], "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0belmuZinv_w3b(        "DY",         fileMap["DYJetsToLL"],      "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0belmuZinv_w3b(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0belmuZinv_w3b(        "single top", fileMap["tW"],              "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0belmuZinv_w3b(       "t#bar{t}Z",  fileMap["TTZ"],             "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0belmuZinv_w3b(        "Rare",       fileMap["Rare"],            "passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w3b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w3b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w3b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w3b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w3b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w3b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagMT2Zinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w3b_elmu("Data",       fileMap["Data_SingleMuon"], "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0b_w3b_elmu(        "DY",         fileMap["DYJetsToLL"],      "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w3b_elmu(      "t#bar{t}",   fileMap["TTbarNoHad"],      "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0b_w3b_elmu(        "single top", fileMap["tW"],              "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w3b_elmu(       "t#bar{t}Z",  fileMap["TTZ"],             "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0b_w3b_elmu(        "Rare",       fileMap["Rare"],            "passBaselineNoTagZinv3b;passElMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");



    // Define the collections, i.e. variables per selection, and group samples in a stack as needed
    // all with muons
    // --> met
    Plotter::DataCollection dcData_2015CD_met_nosel(   "data",   "cleanMetPt", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_met_nosel(            "stack",  "cleanMetPt", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_met_2mu(     "data",   "cleanMetPt", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_met_2mu(              "stack",  "cleanMetPt", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_met_muZinv(  "data",   "cleanMetPt", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_met_muZinv(           "stack",  "cleanMetPt", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_met_ht200(  "data",   "cleanMetPt", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_met_ht200(           "stack",  "cleanMetPt", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_met_blnotag( "data",   "cleanMetPt", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_met_blnotag(          "stack",  "cleanMetPt", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_met_bl(      "data",   "cleanMetPt", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_met_bl(               "stack",  "cleanMetPt", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_met_0bmuZinv("data",   "cleanMetPt", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_met_0bmuZinv(         "stack",  "cleanMetPt", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_met_0bnomt2( "data",   "cleanMetPt", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_met_0bnomt2(          "stack",  "cleanMetPt", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_met_0b(      "data",   "cleanMetPt", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_met_0b(               "stack",  "cleanMetPt", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> ntops
    Plotter::DataCollection dcData_2015CD_nt_nosel("data",  "nTopCandSortedCntZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_nt_nosel("stack",  "nTopCandSortedCntZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nt_2mu("data",  "nTopCandSortedCntZinv", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_nt_2mu("stack",  "nTopCandSortedCntZinv", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_nt_muZinv("data",  "nTopCandSortedCntZinv", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_nt_muZinv("stack",  "nTopCandSortedCntZinv", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_nt_ht200("data",  "nTopCandSortedCntZinv", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_nt_ht200("stack",  "nTopCandSortedCntZinv", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_nt_blnotag("data",  "nTopCandSortedCntZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_nt_blnotag("stack",  "nTopCandSortedCntZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nt_bl("data",  "nTopCandSortedCntZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_nt_bl("stack",  "nTopCandSortedCntZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nt_0bmuZinv("data",  "nTopCandSortedCntZinv", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_nt_0bmuZinv("stack",  "nTopCandSortedCntZinv", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_nt_0bnomt2("data",  "nTopCandSortedCntZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_nt_0bnomt2("stack",  "nTopCandSortedCntZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nt_0b("data",  "nTopCandSortedCntZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_nt_0b("stack",  "nTopCandSortedCntZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> MT2
    Plotter::DataCollection dcData_2015CD_mt2_nosel("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_mt2_nosel("stack",  "best_had_brJet_MT2Zinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mt2_2mu("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mt2_2mu("stack",  "best_had_brJet_MT2Zinv", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mt2_muZinv("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mt2_muZinv("stack",  "best_had_brJet_MT2Zinv", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mt2_ht200("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mt2_ht200("stack",  "best_had_brJet_MT2Zinv", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mt2_blnotag("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mt2_blnotag("stack",  "best_had_brJet_MT2Zinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mt2_bl("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mt2_bl("stack",  "best_had_brJet_MT2Zinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mt2_0bmuZinv("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mt2_0bmuZinv("stack",  "best_had_brJet_MT2Zinv", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mt2_0bnomt2("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mt2_0bnomt2("stack",  "best_had_brJet_MT2Zinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mt2_0b("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mt2_0b("stack",  "best_had_brJet_MT2Zinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> ntops with bfaking
    Plotter::DataCollection dcData_2015CD_nt1b_0bmuZinv("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0bmuZinv_w1b});
    Plotter::DataCollection dcMC_nt1b_0bmuZinv("stack",  "nTopCandSortedCntZinv1b", {dsDY_0bmuZinv_w1b, dstt2l_0bmuZinv_w1b, dstW_0bmuZinv_w1b, dsttZ_0bmuZinv_w1b, dsVV_0bmuZinv_w1b});
    Plotter::DataCollection dcData_2015CD_nt1b_0bnomt2("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0bnomt2_w1b});
    Plotter::DataCollection dcMC_nt1b_0bnomt2("stack",  "nTopCandSortedCntZinv1b", {dsDY_0bnomt2_w1b, dstt2l_0bnomt2_w1b, dstW_0bnomt2_w1b, dsttZ_0bnomt2_w1b, dsVV_0bnomt2_w1b});
    Plotter::DataCollection dcData_2015CD_nt1b_0b("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0b_w1b});
    Plotter::DataCollection dcMC_nt1b_0b("stack",  "nTopCandSortedCntZinv1b", {dsDY_0b_w1b, dstt2l_0b_w1b, dstW_0b_w1b, dsttZ_0b_w1b, dsVV_0b_w1b});
    Plotter::DataCollection dcData_2015CD_nt2b_0bmuZinv("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0bmuZinv_w2b});
    Plotter::DataCollection dcMC_nt2b_0bmuZinv("stack",  "nTopCandSortedCntZinv2b", {dsDY_0bmuZinv_w2b, dstt2l_0bmuZinv_w2b, dstW_0bmuZinv_w2b, dsttZ_0bmuZinv_w2b, dsVV_0bmuZinv_w2b});
    Plotter::DataCollection dcData_2015CD_nt2b_0bnomt2("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0bnomt2_w2b});
    Plotter::DataCollection dcMC_nt2b_0bnomt2("stack",  "nTopCandSortedCntZinv2b", {dsDY_0bnomt2_w2b, dstt2l_0bnomt2_w2b, dstW_0bnomt2_w2b, dsttZ_0bnomt2_w2b, dsVV_0bnomt2_w2b});
    Plotter::DataCollection dcData_2015CD_nt2b_0b("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0b_w2b});
    Plotter::DataCollection dcMC_nt2b_0b("stack",  "nTopCandSortedCntZinv2b", {dsDY_0b_w2b, dstt2l_0b_w2b, dstW_0b_w2b, dsttZ_0b_w2b, dsVV_0b_w2b});
    Plotter::DataCollection dcData_2015CD_nt3b_0bmuZinv("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0bmuZinv_w3b});
    Plotter::DataCollection dcMC_nt3b_0bmuZinv("stack",  "nTopCandSortedCntZinv3b", {dsDY_0bmuZinv_w3b, dstt2l_0bmuZinv_w3b, dstW_0bmuZinv_w3b, dsttZ_0bmuZinv_w3b, dsVV_0bmuZinv_w3b});
    Plotter::DataCollection dcData_2015CD_nt3b_0bnomt2("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0bnomt2_w3b});
    Plotter::DataCollection dcMC_nt3b_0bnomt2("stack",  "nTopCandSortedCntZinv3b", {dsDY_0bnomt2_w3b, dstt2l_0bnomt2_w3b, dstW_0bnomt2_w3b, dsttZ_0bnomt2_w3b, dsVV_0bnomt2_w3b});
    Plotter::DataCollection dcData_2015CD_nt3b_0b("data", "nTopCandSortedCntZinv3b", {dsData_2015C_0b_w3b});
    Plotter::DataCollection dcMC_nt3b_0b("stack",  "nTopCandSortedCntZinv3b", {dsDY_0b_w3b, dstt2l_0b_w3b, dstW_0b_w3b, dsttZ_0b_w3b, dsVV_0b_w3b});
    // --> MT2 with bfaking
    Plotter::DataCollection dcData_2015CD_mt21b_0bmuZinv("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0bmuZinv_w1b});
    Plotter::DataCollection dcMC_mt21b_0bmuZinv("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0bmuZinv_w1b, dstt2l_0bmuZinv_w1b, dstW_0bmuZinv_w1b, dsttZ_0bmuZinv_w1b, dsVV_0bmuZinv_w1b});
    Plotter::DataCollection dcData_2015CD_mt21b_0bnomt2("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0bnomt2_w1b});
    Plotter::DataCollection dcMC_mt21b_0bnomt2("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0bnomt2_w1b, dstt2l_0bnomt2_w1b, dstW_0bnomt2_w1b, dsttZ_0bnomt2_w1b, dsVV_0bnomt2_w1b});
    Plotter::DataCollection dcData_2015CD_mt21b_0b("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0b_w1b});
    Plotter::DataCollection dcMC_mt21b_0b("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0b_w1b, dstt2l_0b_w1b, dstW_0b_w1b, dsttZ_0b_w1b, dsVV_0b_w1b});
    Plotter::DataCollection dcData_2015CD_mt22b_0bmuZinv("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0bmuZinv_w2b});
    Plotter::DataCollection dcMC_mt22b_0bmuZinv("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0bmuZinv_w2b, dstt2l_0bmuZinv_w2b, dstW_0bmuZinv_w2b, dsttZ_0bmuZinv_w2b, dsVV_0bmuZinv_w2b});
    Plotter::DataCollection dcData_2015CD_mt22b_0bnomt2("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0bnomt2_w2b});
    Plotter::DataCollection dcMC_mt22b_0bnomt2("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0bnomt2_w2b, dstt2l_0bnomt2_w2b, dstW_0bnomt2_w2b, dsttZ_0bnomt2_w2b, dsVV_0bnomt2_w2b});
    Plotter::DataCollection dcData_2015CD_mt22b_0b("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0b_w2b});
    Plotter::DataCollection dcMC_mt22b_0b("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0b_w2b, dstt2l_0b_w2b, dstW_0b_w2b, dsttZ_0b_w2b, dsVV_0b_w2b});
    Plotter::DataCollection dcData_2015CD_mt23b_0bmuZinv("data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0bmuZinv_w3b});
    Plotter::DataCollection dcMC_mt23b_0bmuZinv("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0bmuZinv_w3b, dstt2l_0bmuZinv_w3b, dstW_0bmuZinv_w3b, dsttZ_0bmuZinv_w3b, dsVV_0bmuZinv_w3b});
    Plotter::DataCollection dcData_2015CD_mt23b_0bnomt2("data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0bnomt2_w3b});
    Plotter::DataCollection dcMC_mt23b_0bnomt2("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0bnomt2_w3b, dstt2l_0bnomt2_w3b, dstW_0bnomt2_w3b, dsttZ_0bnomt2_w3b, dsVV_0bnomt2_w3b});
    Plotter::DataCollection dcData_2015CD_mt23b_0b("data", "best_had_brJet_MT2Zinv3b", {dsData_2015C_0b_w3b});
    Plotter::DataCollection dcMC_mt23b_0b("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0b_w3b, dstt2l_0b_w3b, dstW_0b_w3b, dsttZ_0b_w3b, dsVV_0b_w3b});
    // --> nbs
    Plotter::DataCollection dcData_2015CD_nb_nosel("data",  "cntCSVSZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_nb_nosel("stack",  "cntCSVSZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nb_2mu("data",  "cntCSVSZinv", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_nb_2mu("stack",  "cntCSVSZinv", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_nb_muZinv("data",  "cntCSVSZinv", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_nb_muZinv("stack",  "cntCSVSZinv", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_nb_ht200("data",  "cntCSVSZinv", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_nb_ht200("stack",  "cntCSVSZinv", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_nb_blnotag("data",  "cntCSVSZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_nb_blnotag("stack",  "cntCSVSZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nb_bl("data",  "cntCSVSZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_nb_bl("stack",  "cntCSVSZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nb_0bnomt2("data",  "cntCSVSZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_nb_0bnomt2("stack",  "cntCSVSZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nb_0b("data",  "cntCSVSZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_nb_0b("stack",  "cntCSVSZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> njets
    Plotter::DataCollection dcData_2015CD_nj_nosel("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_nj_nosel("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nj_2mu("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_nj_2mu("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_nj_muZinv("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_nj_muZinv("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_nj_ht200("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_nj_ht200("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_nj_blnotag("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_nj_blnotag("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nj_bl("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_nj_bl("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nj_0bmuZinv("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_nj_0bmuZinv("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_nj_0bnomt2("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_nj_0bnomt2("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nj_0b("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_nj_0b("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> HT
    Plotter::DataCollection dcData_2015CD_ht_nosel("data",  "HTZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_ht_nosel("stack",  "HTZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_ht_2mu("data",  "HTZinv", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_ht_2mu("stack",  "HTZinv", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_ht_muZinv("data",  "HTZinv", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_ht_muZinv("stack",  "HTZinv", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_ht_ht200("data",  "HTZinv", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_ht_ht200("stack",  "HTZinv", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_ht_blnotag("data",  "HTZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_ht_blnotag("stack",  "HTZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_ht_bl("data",  "HTZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_ht_bl("stack",  "HTZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_ht_0bmuZinv("data",  "HTZinv", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_ht_0bmuZinv("stack",  "HTZinv", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_ht_0bnomt2("data",  "HTZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_ht_0bnomt2("stack",  "HTZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_ht_0b("data",  "HTZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_ht_0b("stack",  "HTZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> MHT
    Plotter::DataCollection dcData_2015CD_mht_nosel("data",  "cleanMHt", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_mht_nosel("stack",  "cleanMHt", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mht_2mu("data",  "cleanMHt", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mht_2mu("stack",  "cleanMHt", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mht_muZinv("data",  "cleanMHt", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mht_muZinv("stack",  "cleanMHt", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mht_ht200("data",  "cleanMHt", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mht_ht200("stack",  "cleanMHt", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mht_blnotag("data",  "cleanMHt", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mht_blnotag("stack",  "cleanMHt", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mht_bl("data",  "cleanMHt", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mht_bl("stack",  "cleanMHt", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mht_0bmuZinv("data",  "cleanMHt", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mht_0bmuZinv("stack",  "cleanMHt", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mht_0bnomt2("data",  "cleanMHt", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mht_0bnomt2("stack",  "cleanMHt", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mht_0b("data",  "cleanMHt", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mht_0b("stack",  "cleanMHt", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet pt
    Plotter::DataCollection dcData_2015CD_jpt_nosel("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_jpt_nosel("stack",  "jetsLVecLepCleaned(pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_jpt_2mu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_jpt_2mu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_jpt_muZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_jpt_muZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_jpt_ht200("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_jpt_ht200("stack",  "jetsLVecLepCleaned(pt)", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_jpt_blnotag("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_jpt_blnotag("stack",  "jetsLVecLepCleaned(pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_jpt_bl("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_jpt_bl("stack",  "jetsLVecLepCleaned(pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_jpt_0bmuZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_jpt_0bmuZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_jpt_0bnomt2("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_jpt_0bnomt2("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_jpt_0b("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_jpt_0b("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet1 pt 
    Plotter::DataCollection dcData_2015CD_j1pt_nosel("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_j1pt_nosel("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j1pt_2mu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_j1pt_2mu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_j1pt_muZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_j1pt_muZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_ht200("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_j1pt_ht200("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_j1pt_blnotag("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_j1pt_blnotag("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j1pt_bl("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_j1pt_bl("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j1pt_0bmuZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_j1pt_0bmuZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_0bnomt2("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_j1pt_0bnomt2("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j1pt_0b("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_j1pt_0b("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet2 pt
    Plotter::DataCollection dcData_2015CD_j2pt_nosel("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_j2pt_nosel("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j2pt_2mu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_j2pt_2mu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_j2pt_muZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_j2pt_muZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_ht200("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_j2pt_ht200("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_j2pt_blnotag("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_j2pt_blnotag("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j2pt_bl("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_j2pt_bl("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j2pt_0bmuZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_j2pt_0bmuZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_0bnomt2("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_j2pt_0bnomt2("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j2pt_0b("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_j2pt_0b("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet3 pt
    Plotter::DataCollection dcData_2015CD_j3pt_nosel("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_j3pt_nosel("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j3pt_2mu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_j3pt_2mu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_j3pt_muZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_j3pt_muZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_ht200("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_j3pt_ht200("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_j3pt_blnotag("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_j3pt_blnotag("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j3pt_bl("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_j3pt_bl("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j3pt_0bmuZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_j3pt_0bmuZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_0bnomt2("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_j3pt_0bnomt2("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j3pt_0b("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_j3pt_0b("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> mu pt
    Plotter::DataCollection dcData_2015CD_mupt_nosel("data",  "cutMuVec(pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_mupt_nosel("stack",  "cutMuVec(pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mupt_2mu("data",  "cutMuVec(pt)", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mupt_2mu("stack",  "cutMuVec(pt)", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mupt_muZinv("data",  "cutMuVec(pt)", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mupt_muZinv("stack",  "cutMuVec(pt)", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mupt_ht200("data",  "cutMuVec(pt)", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mupt_ht200("stack",  "cutMuVec(pt)", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mupt_blnotag("data",  "cutMuVec(pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mupt_blnotag("stack",  "cutMuVec(pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mupt_bl("data",  "cutMuVec(pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mupt_bl("stack",  "cutMuVec(pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mupt_0bmuZinv("data",  "cutMuVec(pt)", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mupt_0bmuZinv("stack",  "cutMuVec(pt)", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mupt_0bnomt2("data",  "cutMuVec(pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mupt_0bnomt2("stack",  "cutMuVec(pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mupt_0b("data",  "cutMuVec(pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mupt_0b("stack",  "cutMuVec(pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> mu pt
    Plotter::DataCollection dcData_2015CD_mu1pt_nosel("data",  "cutMuPt1", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_mu1pt_nosel("stack",  "cutMuPt1", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mu1pt_2mu("data",  "cutMuPt1", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mu1pt_2mu("stack",  "cutMuPt1", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mu1pt_muZinv("data",  "cutMuPt1", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mu1pt_muZinv("stack",  "cutMuPt1", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mu1pt_ht200("data",  "cutMuPt1", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mu1pt_ht200("stack",  "cutMuPt1", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mu1pt_blnotag("data",  "cutMuPt1", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mu1pt_blnotag("stack",  "cutMuPt1", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mu1pt_bl("data",  "cutMuPt1", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mu1pt_bl("stack",  "cutMuPt1", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mu1pt_0bmuZinv("data",  "cutMuPt1", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mu1pt_0bmuZinv("stack",  "cutMuPt1", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mu1pt_0bnomt2("data",  "cutMuPt1", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mu1pt_0bnomt2("stack",  "cutMuPt1", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mu1pt_0b("data",  "cutMuPt1", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mu1pt_0b("stack",  "cutMuPt1", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> mu pt
    Plotter::DataCollection dcData_2015CD_mu2pt_nosel("data",  "cutMuPt2", {dsData_2015C_nosel});
    Plotter::DataCollection dcMC_mu2pt_nosel("stack",  "cutMuPt2", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mu2pt_2mu("data",  "cutMuPt2", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mu2pt_2mu("stack",  "cutMuPt2", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mu2pt_muZinv("data",  "cutMuPt2", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mu2pt_muZinv("stack",  "cutMuPt2", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mu2pt_ht200("data",  "cutMuPt2", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mu2pt_ht200("stack",  "cutMuPt2", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mu2pt_blnotag("data",  "cutMuPt2", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mu2pt_blnotag("stack",  "cutMuPt2", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mu2pt_bl("data",  "cutMuPt2", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mu2pt_bl("stack",  "cutMuPt2", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mu2pt_0bmuZinv("data",  "cutMuPt2", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mu2pt_0bmuZinv("stack",  "cutMuPt2", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mu2pt_0bnomt2("data",  "cutMuPt2", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mu2pt_0bnomt2("stack",  "cutMuPt2", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mu2pt_0b("data",  "cutMuPt2", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mu2pt_0b("stack",  "cutMuPt2", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> mll
    Plotter::DataCollection dcData_2015CD_mll_2mu("data",  "bestRecoZM", {dsData_2015C_2mu});
    Plotter::DataCollection dcMC_mll_2mu("stack",  "bestRecoZM", {dsDY_2mu, dstt2l_2mu, dstW_2mu, dsttZ_2mu, dsVV_2mu});
    Plotter::DataCollection dcData_2015CD_mll_muZinv("data",  "bestRecoZM", {dsData_2015C_muZinv});
    Plotter::DataCollection dcMC_mll_muZinv("stack",  "bestRecoZM", {dsDY_muZinv, dstt2l_muZinv, dstW_muZinv, dsttZ_muZinv, dsVV_muZinv});
    Plotter::DataCollection dcData_2015CD_mll_ht200("data",  "bestRecoZM", {dsData_2015C_ht200});
    Plotter::DataCollection dcMC_mll_ht200("stack",  "bestRecoZM", {dsDY_ht200, dstt2l_ht200, dstW_ht200, dsttZ_ht200, dsVV_ht200});
    Plotter::DataCollection dcData_2015CD_mll_blnotag("data",  "bestRecoZM", {dsData_2015C_blnotag});
    Plotter::DataCollection dcMC_mll_blnotag("stack",  "bestRecoZM", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mll_bl("data",  "bestRecoZM", {dsData_2015C_bl});
    Plotter::DataCollection dcMC_mll_bl("stack",  "bestRecoZM", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mll_0bmuZinv("data",  "bestRecoZM", {dsData_2015C_0bmuZinv});
    Plotter::DataCollection dcMC_mll_0bmuZinv("stack",  "bestRecoZM", {dsDY_0bmuZinv, dstt2l_0bmuZinv, dstW_0bmuZinv, dsttZ_0bmuZinv, dsVV_0bmuZinv});
    Plotter::DataCollection dcData_2015CD_mll_0bnomt2("data",  "bestRecoZM", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcMC_mll_0bnomt2("stack",  "bestRecoZM", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mll_0b("data",  "bestRecoZM", {dsData_2015C_0b});
    Plotter::DataCollection dcMC_mll_0b("stack",  "bestRecoZM", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});


    // all with electrons
    // --> met
    Plotter::DataCollection dcData_2015CD_met_2el(        "data",   "cleanMetPt", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_met_2el(                 "stack",  "cleanMetPt", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_met_elZinv(     "data",   "cleanMetPt", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_met_elZinv(              "stack",  "cleanMetPt", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_met_ht200_el(   "data",   "cleanMetPt", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_met_ht200_el(            "stack",  "cleanMetPt", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_met_blnotag_el( "data",   "cleanMetPt", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_met_blnotag_el(          "stack",  "cleanMetPt", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_met_bl_el(      "data",   "cleanMetPt", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_met_bl_el(               "stack",  "cleanMetPt", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_met_0belZinv(   "data",   "cleanMetPt", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_met_0belZinv(            "stack",  "cleanMetPt", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_met_0bnomt2_el( "data",   "cleanMetPt", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_met_0bnomt2_el(          "stack",  "cleanMetPt", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_met_0b_el(      "data",   "cleanMetPt", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_met_0b_el(               "stack",  "cleanMetPt", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> ntops
    Plotter::DataCollection dcData_2015CD_nt_2el("data",  "nTopCandSortedCntZinv", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_nt_2el("stack",  "nTopCandSortedCntZinv", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_nt_elZinv("data",  "nTopCandSortedCntZinv", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_nt_elZinv("stack",  "nTopCandSortedCntZinv", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_nt_ht200_el("data",  "nTopCandSortedCntZinv", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_nt_ht200_el("stack",  "nTopCandSortedCntZinv", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_nt_blnotag_el("data",  "nTopCandSortedCntZinv", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_nt_blnotag_el("stack",  "nTopCandSortedCntZinv", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_nt_bl_el("data",  "nTopCandSortedCntZinv", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_nt_bl_el("stack",  "nTopCandSortedCntZinv", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_nt_0belZinv("data",  "nTopCandSortedCntZinv", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_nt_0belZinv("stack",  "nTopCandSortedCntZinv", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_nt_0bnomt2_el("data",  "nTopCandSortedCntZinv", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_nt_0bnomt2_el("stack",  "nTopCandSortedCntZinv", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_nt_0b_el("data",  "nTopCandSortedCntZinv", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_nt_0b_el("stack",  "nTopCandSortedCntZinv", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> MT2
    Plotter::DataCollection dcData_2015CD_mt2_2el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_mt2_2el("stack",  "best_had_brJet_MT2Zinv", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_mt2_elZinv("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_mt2_elZinv("stack",  "best_had_brJet_MT2Zinv", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_mt2_ht200_el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_mt2_ht200_el("stack",  "best_had_brJet_MT2Zinv", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_mt2_blnotag_el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_mt2_blnotag_el("stack",  "best_had_brJet_MT2Zinv", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_mt2_bl_el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_mt2_bl_el("stack",  "best_had_brJet_MT2Zinv", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_mt2_0belZinv("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_mt2_0belZinv("stack",  "best_had_brJet_MT2Zinv", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_mt2_0bnomt2_el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_mt2_0bnomt2_el("stack",  "best_had_brJet_MT2Zinv", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_mt2_0b_el("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_mt2_0b_el("stack",  "best_had_brJet_MT2Zinv", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> ntops with bfaking
    Plotter::DataCollection dcData_2015CD_nt1b_0belZinv("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0belZinv_w1b});
    Plotter::DataCollection dcMC_nt1b_0belZinv("stack",  "nTopCandSortedCntZinv1b", {dsDY_0belZinv_w1b, dstt2l_0belZinv_w1b, dstW_0belZinv_w1b, dsttZ_0belZinv_w1b, dsVV_0belZinv_w1b});
    Plotter::DataCollection dcData_2015CD_nt1b_0bnomt2_el("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0bnomt2_w1b_el});
    Plotter::DataCollection dcMC_nt1b_0bnomt2_el("stack",  "nTopCandSortedCntZinv1b", {dsDY_0bnomt2_w1b_el, dstt2l_0bnomt2_w1b_el, dstW_0bnomt2_w1b_el, dsttZ_0bnomt2_w1b_el, dsVV_0bnomt2_w1b_el});
    Plotter::DataCollection dcData_2015CD_nt1b_0b_el("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0b_w1b_el});
    Plotter::DataCollection dcMC_nt1b_0b_el("stack",  "nTopCandSortedCntZinv1b", {dsDY_0b_w1b_el, dstt2l_0b_w1b_el, dstW_0b_w1b_el, dsttZ_0b_w1b_el, dsVV_0b_w1b_el});
    Plotter::DataCollection dcData_2015CD_nt2b_0belZinv("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0belZinv_w2b});
    Plotter::DataCollection dcMC_nt2b_0belZinv("stack",  "nTopCandSortedCntZinv2b", {dsDY_0belZinv_w2b, dstt2l_0belZinv_w2b, dstW_0belZinv_w2b, dsttZ_0belZinv_w2b, dsVV_0belZinv_w2b});
    Plotter::DataCollection dcData_2015CD_nt2b_0bnomt2_el("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0bnomt2_w2b_el});
    Plotter::DataCollection dcMC_nt2b_0bnomt2_el("stack",  "nTopCandSortedCntZinv2b", {dsDY_0bnomt2_w2b_el, dstt2l_0bnomt2_w2b_el, dstW_0bnomt2_w2b_el, dsttZ_0bnomt2_w2b_el, dsVV_0bnomt2_w2b_el});
    Plotter::DataCollection dcData_2015CD_nt2b_0b_el("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0b_w2b_el});
    Plotter::DataCollection dcMC_nt2b_0b_el("stack",  "nTopCandSortedCntZinv2b", {dsDY_0b_w2b_el, dstt2l_0b_w2b_el, dstW_0b_w2b_el, dsttZ_0b_w2b_el, dsVV_0b_w2b_el});
    Plotter::DataCollection dcData_2015CD_nt3b_0belZinv("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0belZinv_w3b});
    Plotter::DataCollection dcMC_nt3b_0belZinv("stack",  "nTopCandSortedCntZinv3b", {dsDY_0belZinv_w3b, dstt2l_0belZinv_w3b, dstW_0belZinv_w3b, dsttZ_0belZinv_w3b, dsVV_0belZinv_w3b});
    Plotter::DataCollection dcData_2015CD_nt3b_0bnomt2_el("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0bnomt2_w3b_el});
    Plotter::DataCollection dcMC_nt3b_0bnomt2_el("stack",  "nTopCandSortedCntZinv3b", {dsDY_0bnomt2_w3b_el, dstt2l_0bnomt2_w3b_el, dstW_0bnomt2_w3b_el, dsttZ_0bnomt2_w3b_el, dsVV_0bnomt2_w3b_el});
    Plotter::DataCollection dcData_2015CD_nt3b_0b_el("data", "nTopCandSortedCntZinv3b", {dsData_2015C_0b_w3b_el});
    Plotter::DataCollection dcMC_nt3b_0b_el("stack",  "nTopCandSortedCntZinv3b", {dsDY_0b_w3b_el, dstt2l_0b_w3b_el, dstW_0b_w3b_el, dsttZ_0b_w3b_el, dsVV_0b_w3b_el});
    // --> MT2 with bfaking
    Plotter::DataCollection dcData_2015CD_mt21b_0belZinv("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0belZinv_w1b});
    Plotter::DataCollection dcMC_mt21b_0belZinv("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0belZinv_w1b, dstt2l_0belZinv_w1b, dstW_0belZinv_w1b, dsttZ_0belZinv_w1b, dsVV_0belZinv_w1b});
    Plotter::DataCollection dcData_2015CD_mt21b_0bnomt2_el("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0bnomt2_w1b_el});
    Plotter::DataCollection dcMC_mt21b_0bnomt2_el("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0bnomt2_w1b_el, dstt2l_0bnomt2_w1b_el, dstW_0bnomt2_w1b_el, dsttZ_0bnomt2_w1b_el, dsVV_0bnomt2_w1b_el});
    Plotter::DataCollection dcData_2015CD_mt21b_0b_el("data", "best_had_brJet_MT2Zinv1b", {dsData_2015C_0b_w1b_el});
    Plotter::DataCollection dcMC_mt21b_0b_el("stack",  "best_had_brJet_MT2Zinv1b", {dsDY_0b_w1b_el, dstt2l_0b_w1b_el, dstW_0b_w1b_el, dsttZ_0b_w1b_el, dsVV_0b_w1b_el});
    Plotter::DataCollection dcData_2015CD_mt22b_0belZinv("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0belZinv_w2b});
    Plotter::DataCollection dcMC_mt22b_0belZinv("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0belZinv_w2b, dstt2l_0belZinv_w2b, dstW_0belZinv_w2b, dsttZ_0belZinv_w2b, dsVV_0belZinv_w2b});
    Plotter::DataCollection dcData_2015CD_mt22b_0bnomt2_el("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0bnomt2_w2b_el});
    Plotter::DataCollection dcMC_mt22b_0bnomt2_el("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0bnomt2_w2b_el, dstt2l_0bnomt2_w2b_el, dstW_0bnomt2_w2b_el, dsttZ_0bnomt2_w2b_el, dsVV_0bnomt2_w2b_el});
    Plotter::DataCollection dcData_2015CD_mt22b_0b_el("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0b_w2b_el});
    Plotter::DataCollection dcMC_mt22b_0b_el("stack",  "best_had_brJet_MT2Zinv2b", {dsDY_0b_w2b_el, dstt2l_0b_w2b_el, dstW_0b_w2b_el, dsttZ_0b_w2b_el, dsVV_0b_w2b_el});
    Plotter::DataCollection dcData_2015CD_mt23b_0belZinv("data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0belZinv_w3b});
    Plotter::DataCollection dcMC_mt23b_0belZinv("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0belZinv_w3b, dstt2l_0belZinv_w3b, dstW_0belZinv_w3b, dsttZ_0belZinv_w3b, dsVV_0belZinv_w3b});
    Plotter::DataCollection dcData_2015CD_mt23b_0bnomt2_el("data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0bnomt2_w3b_el});
    Plotter::DataCollection dcMC_mt23b_0bnomt2_el("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0bnomt2_w3b_el, dstt2l_0bnomt2_w3b_el, dstW_0bnomt2_w3b_el, dsttZ_0bnomt2_w3b_el, dsVV_0bnomt2_w3b_el});
    Plotter::DataCollection dcData_2015CD_mt23b_0b_el("data", "best_had_brJet_MT2Zinv3b", {dsData_2015C_0b_w3b_el});
    Plotter::DataCollection dcMC_mt23b_0b_el("stack",  "best_had_brJet_MT2Zinv3b", {dsDY_0b_w3b_el, dstt2l_0b_w3b_el, dstW_0b_w3b_el, dsttZ_0b_w3b_el, dsVV_0b_w3b_el});
    // --> nbs
    Plotter::DataCollection dcData_2015CD_nb_2el("data",  "cntCSVSZinv", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_nb_2el("stack",  "cntCSVSZinv", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_nb_elZinv("data",  "cntCSVSZinv", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_nb_elZinv("stack",  "cntCSVSZinv", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_nb_ht200_el("data",  "cntCSVSZinv", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_nb_ht200_el("stack",  "cntCSVSZinv", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_nb_blnotag_el("data",  "cntCSVSZinv", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_nb_blnotag_el("stack",  "cntCSVSZinv", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_nb_bl_el("data",  "cntCSVSZinv", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_nb_bl_el("stack",  "cntCSVSZinv", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_nb_0bnomt2_el("data",  "cntCSVSZinv", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_nb_0bnomt2_el("stack",  "cntCSVSZinv", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_nb_0b_el("data",  "cntCSVSZinv", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_nb_0b_el("stack",  "cntCSVSZinv", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> njets
    Plotter::DataCollection dcData_2015CD_nj_2el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_nj_2el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_nj_elZinv("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_nj_elZinv("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_nj_ht200_el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_nj_ht200_el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_nj_blnotag_el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_nj_blnotag_el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_nj_bl_el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_nj_bl_el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_nj_0belZinv("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_nj_0belZinv("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_nj_0bnomt2_el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_nj_0bnomt2_el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_nj_0b_el("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_nj_0b_el("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> HT
    Plotter::DataCollection dcData_2015CD_ht_2el("data",  "HTZinv", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_ht_2el("stack",  "HTZinv", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_ht_elZinv("data",  "HTZinv", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_ht_elZinv("stack",  "HTZinv", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_ht_ht200_el("data",  "HTZinv", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_ht_ht200_el("stack",  "HTZinv", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_ht_blnotag_el("data",  "HTZinv", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_ht_blnotag_el("stack",  "HTZinv", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_ht_bl_el("data",  "HTZinv", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_ht_bl_el("stack",  "HTZinv", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_ht_0belZinv("data",  "HTZinv", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_ht_0belZinv("stack",  "HTZinv", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_ht_0bnomt2_el("data",  "HTZinv", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_ht_0bnomt2_el("stack",  "HTZinv", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_ht_0b_el("data",  "HTZinv", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_ht_0b_el("stack",  "HTZinv", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> MHT
    Plotter::DataCollection dcData_2015CD_mht_2el("data",  "cleanMHt", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_mht_2el("stack",  "cleanMHt", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_mht_elZinv("data",  "cleanMHt", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_mht_elZinv("stack",  "cleanMHt", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_mht_ht200_el("data",  "cleanMHt", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_mht_ht200_el("stack",  "cleanMHt", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_mht_blnotag_el("data",  "cleanMHt", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_mht_blnotag_el("stack",  "cleanMHt", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_mht_bl_el("data",  "cleanMHt", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_mht_bl_el("stack",  "cleanMHt", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_mht_0belZinv("data",  "cleanMHt", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_mht_0belZinv("stack",  "cleanMHt", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_mht_0bnomt2_el("data",  "cleanMHt", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_mht_0bnomt2_el("stack",  "cleanMHt", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_mht_0b_el("data",  "cleanMHt", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_mht_0b_el("stack",  "cleanMHt", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> jet pt
    Plotter::DataCollection dcData_2015CD_jpt_2el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_jpt_2el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_jpt_elZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_jpt_elZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_jpt_ht200_el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_jpt_ht200_el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_jpt_blnotag_el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_jpt_blnotag_el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_jpt_bl_el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_jpt_bl_el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_jpt_0belZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_jpt_0belZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_jpt_0bnomt2_el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_jpt_0bnomt2_el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_jpt_0b_el("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_jpt_0b_el("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> jet1 pt 
    Plotter::DataCollection dcData_2015CD_j1pt_2el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_j1pt_2el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_j1pt_elZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_j1pt_elZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_ht200_el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_j1pt_ht200_el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_j1pt_blnotag_el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_j1pt_blnotag_el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_j1pt_bl_el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_j1pt_bl_el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_j1pt_0belZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_j1pt_0belZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_0bnomt2_el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_j1pt_0bnomt2_el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_j1pt_0b_el("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_j1pt_0b_el("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> jet2 pt
    Plotter::DataCollection dcData_2015CD_j2pt_2el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_j2pt_2el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_j2pt_elZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_j2pt_elZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_ht200_el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_j2pt_ht200_el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_j2pt_blnotag_el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_j2pt_blnotag_el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_j2pt_bl_el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_j2pt_bl_el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_j2pt_0belZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_j2pt_0belZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_0bnomt2_el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_j2pt_0bnomt2_el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_j2pt_0b_el("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_j2pt_0b_el("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> jet3 pt
    Plotter::DataCollection dcData_2015CD_j3pt_2el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_j3pt_2el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_j3pt_elZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_j3pt_elZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_ht200_el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_j3pt_ht200_el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_j3pt_blnotag_el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_j3pt_blnotag_el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_j3pt_bl_el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_j3pt_bl_el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_j3pt_0belZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_j3pt_0belZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_0bnomt2_el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_j3pt_0bnomt2_el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_j3pt_0b_el("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_j3pt_0b_el("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> el pt
    Plotter::DataCollection dcData_2015CD_elpt_2el("data",  "cutElecVec(pt)", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_elpt_2el("stack",  "cutElecVec(pt)", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_elpt_elZinv("data",  "cutElecVec(pt)", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_elpt_elZinv("stack",  "cutElecVec(pt)", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_elpt_ht200("data",  "cutElecVec(pt)", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_elpt_ht200("stack",  "cutElecVec(pt)", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_elpt_blnotag("data",  "cutElecVec(pt)", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_elpt_blnotag("stack",  "cutElecVec(pt)", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_elpt_bl("data",  "cutElecVec(pt)", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_elpt_bl("stack",  "cutElecVec(pt)", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_elpt_0belZinv("data",  "cutElecVec(pt)", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_elpt_0belZinv("stack",  "cutElecVec(pt)", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_elpt_0bnomt2("data",  "cutElecVec(pt)", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_elpt_0bnomt2("stack",  "cutElecVec(pt)", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_elpt_0b("data",  "cutElecVec(pt)", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_elpt_0b("stack",  "cutElecVec(pt)", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> el pt
    Plotter::DataCollection dcData_2015CD_el1pt_2el("data",  "cutElecPt1", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_el1pt_2el("stack",  "cutElecPt1", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_el1pt_elZinv("data",  "cutElecPt1", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_el1pt_elZinv("stack",  "cutElecPt1", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_el1pt_ht200("data",  "cutElecPt1", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_el1pt_ht200("stack",  "cutElecPt1", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_el1pt_blnotag("data",  "cutElecPt1", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_el1pt_blnotag("stack",  "cutElecPt1", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_el1pt_bl("data",  "cutElecPt1", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_el1pt_bl("stack",  "cutElecPt1", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_el1pt_0belZinv("data",  "cutElecPt1", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_el1pt_0belZinv("stack",  "cutElecPt1", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_el1pt_0bnomt2("data",  "cutElecPt1", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_el1pt_0bnomt2("stack",  "cutElecPt1", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_el1pt_0b("data",  "cutElecPt1", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_el1pt_0b("stack",  "cutElecPt1", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> el pt
    Plotter::DataCollection dcData_2015CD_el2pt_2el("data",  "cutElecPt2", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_el2pt_2el("stack",  "cutElecPt2", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_el2pt_elZinv("data",  "cutElecPt2", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_el2pt_elZinv("stack",  "cutElecPt2", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_el2pt_ht200("data",  "cutElecPt2", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_el2pt_ht200("stack",  "cutElecPt2", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_el2pt_blnotag("data",  "cutElecPt2", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_el2pt_blnotag("stack",  "cutElecPt2", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_el2pt_bl("data",  "cutElecPt2", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_el2pt_bl("stack",  "cutElecPt2", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_el2pt_0belZinv("data",  "cutElecPt2", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_el2pt_0belZinv("stack",  "cutElecPt2", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_el2pt_0bnomt2("data",  "cutElecPt2", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_el2pt_0bnomt2("stack",  "cutElecPt2", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_el2pt_0b("data",  "cutElecPt2", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_el2pt_0b("stack",  "cutElecPt2", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});
    // --> mll
    Plotter::DataCollection dcData_2015CD_mll_2el("data",  "bestRecoZM", {dsData_2015C_2el});
    Plotter::DataCollection dcMC_mll_2el("stack",  "bestRecoZM", {dsDY_2el, dstt2l_2el, dstW_2el, dsttZ_2el, dsVV_2el});
    Plotter::DataCollection dcData_2015CD_mll_elZinv("data",  "bestRecoZM", {dsData_2015C_elZinv});
    Plotter::DataCollection dcMC_mll_elZinv("stack",  "bestRecoZM", {dsDY_elZinv, dstt2l_elZinv, dstW_elZinv, dsttZ_elZinv, dsVV_elZinv});
    Plotter::DataCollection dcData_2015CD_mll_ht200_el("data",  "bestRecoZM", {dsData_2015C_ht200_el});
    Plotter::DataCollection dcMC_mll_ht200_el("stack",  "bestRecoZM", {dsDY_ht200_el, dstt2l_ht200_el, dstW_ht200_el, dsttZ_ht200_el, dsVV_ht200_el});
    Plotter::DataCollection dcData_2015CD_mll_blnotag_el("data",  "bestRecoZM", {dsData_2015C_blnotag_el});
    Plotter::DataCollection dcMC_mll_blnotag_el("stack",  "bestRecoZM", {dsDY_blnotag_el, dstt2l_blnotag_el, dstW_blnotag_el, dsttZ_blnotag_el, dsVV_blnotag_el});
    Plotter::DataCollection dcData_2015CD_mll_bl_el("data",  "bestRecoZM", {dsData_2015C_bl_el});
    Plotter::DataCollection dcMC_mll_bl_el("stack",  "bestRecoZM", {dsDY_bl_el, dstt2l_bl_el, dstW_bl_el, dsttZ_bl_el, dsVV_bl_el});
    Plotter::DataCollection dcData_2015CD_mll_0belZinv("data",  "bestRecoZM", {dsData_2015C_0belZinv});
    Plotter::DataCollection dcMC_mll_0belZinv("stack",  "bestRecoZM", {dsDY_0belZinv, dstt2l_0belZinv, dstW_0belZinv, dsttZ_0belZinv, dsVV_0belZinv});
    Plotter::DataCollection dcData_2015CD_mll_0bnomt2_el("data",  "bestRecoZM", {dsData_2015C_0bnomt2_el});
    Plotter::DataCollection dcMC_mll_0bnomt2_el("stack",  "bestRecoZM", {dsDY_0bnomt2_el, dstt2l_0bnomt2_el, dstW_0bnomt2_el, dsttZ_0bnomt2_el, dsVV_0bnomt2_el});
    Plotter::DataCollection dcData_2015CD_mll_0b_el("data",  "bestRecoZM", {dsData_2015C_0b_el});
    Plotter::DataCollection dcMC_mll_0b_el("stack",  "bestRecoZM", {dsDY_0b_el, dstt2l_0b_el, dstW_0b_el, dsttZ_0b_el, dsVV_0b_el});


    // all with e+mu
    // --> met
    Plotter::DataCollection dcData_2015CD_met_elmu(        "data",  "cleanMetPt", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_met_elmu(                 "stack", "cleanMetPt", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_met_elmuZinv(    "data",  "cleanMetPt", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_met_elmuZinv(             "stack", "cleanMetPt", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_met_ht200_elmu(  "data",  "cleanMetPt", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_met_ht200_elmu(           "stack", "cleanMetPt", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_met_blnotag_elmu("data",  "cleanMetPt", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_met_blnotag_elmu(         "stack", "cleanMetPt", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_met_bl_elmu(     "data",  "cleanMetPt", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_met_bl_elmu(              "stack", "cleanMetPt", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_met_0belmuZinv(  "data",  "cleanMetPt", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_met_0belmuZinv(           "stack", "cleanMetPt", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_met_0bnomt2_elmu("data",  "cleanMetPt", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_met_0bnomt2_elmu(         "stack", "cleanMetPt", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_met_0b_elmu(     "data",  "cleanMetPt", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_met_0b_elmu(              "stack", "cleanMetPt", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> ntops
    Plotter::DataCollection dcData_2015CD_nt_elmu(        "data",  "nTopCandSortedCntZinv", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_nt_elmu(                 "stack", "nTopCandSortedCntZinv", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_nt_elmuZinv(    "data",  "nTopCandSortedCntZinv", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_nt_elmuZinv(             "stack", "nTopCandSortedCntZinv", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_nt_ht200_elmu(  "data",  "nTopCandSortedCntZinv", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_nt_ht200_elmu(           "stack", "nTopCandSortedCntZinv", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_nt_blnotag_elmu("data",  "nTopCandSortedCntZinv", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_nt_blnotag_elmu(         "stack", "nTopCandSortedCntZinv", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_nt_bl_elmu(     "data",  "nTopCandSortedCntZinv", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_nt_bl_elmu(              "stack", "nTopCandSortedCntZinv", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_nt_0belmuZinv(  "data",  "nTopCandSortedCntZinv", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_nt_0belmuZinv(           "stack", "nTopCandSortedCntZinv", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_nt_0bnomt2_elmu("data",  "nTopCandSortedCntZinv", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_nt_0bnomt2_elmu(         "stack", "nTopCandSortedCntZinv", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_nt_0b_elmu(     "data",  "nTopCandSortedCntZinv", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_nt_0b_elmu(              "stack", "nTopCandSortedCntZinv", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> MT2
    Plotter::DataCollection dcData_2015CD_mt2_elmu(        "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_mt2_elmu(                 "stack", "best_had_brJet_MT2Zinv", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_mt2_elmuZinv(    "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_mt2_elmuZinv(             "stack", "best_had_brJet_MT2Zinv", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_mt2_ht200_elmu(  "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_mt2_ht200_elmu(           "stack", "best_had_brJet_MT2Zinv", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_mt2_blnotag_elmu("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_mt2_blnotag_elmu(         "stack", "best_had_brJet_MT2Zinv", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_mt2_bl_elmu(     "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_mt2_bl_elmu(              "stack", "best_had_brJet_MT2Zinv", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_mt2_0belmuZinv(  "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_mt2_0belmuZinv(           "stack", "best_had_brJet_MT2Zinv", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_mt2_0bnomt2_elmu("data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_mt2_0bnomt2_elmu(         "stack", "best_had_brJet_MT2Zinv", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_mt2_0b_elmu(     "data",  "best_had_brJet_MT2Zinv", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_mt2_0b_elmu(              "stack", "best_had_brJet_MT2Zinv", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> ntops with bfaking
    Plotter::DataCollection dcData_2015CD_nt1b_0belmuZinv(  "data",  "nTopCandSortedCntZinv1b", {dsData_2015C_0belmuZinv_w1b});
    Plotter::DataCollection dcMC_nt1b_0belmuZinv(           "stack", "nTopCandSortedCntZinv1b", {dsDY_0belmuZinv_w1b, dstt2l_0belmuZinv_w1b, dstW_0belmuZinv_w1b, dsttZ_0belmuZinv_w1b, dsVV_0belmuZinv_w1b});
    Plotter::DataCollection dcData_2015CD_nt1b_0bnomt2_elmu("data",  "nTopCandSortedCntZinv1b", {dsData_2015C_0bnomt2_w1b_elmu});
    Plotter::DataCollection dcMC_nt1b_0bnomt2_elmu(         "stack", "nTopCandSortedCntZinv1b", {dsDY_0bnomt2_w1b_elmu, dstt2l_0bnomt2_w1b_elmu, dstW_0bnomt2_w1b_elmu, dsttZ_0bnomt2_w1b_elmu, dsVV_0bnomt2_w1b_elmu});
    Plotter::DataCollection dcData_2015CD_nt1b_0b_elmu(     "data",  "nTopCandSortedCntZinv1b", {dsData_2015C_0b_w1b_elmu});
    Plotter::DataCollection dcMC_nt1b_0b_elmu(              "stack", "nTopCandSortedCntZinv1b", {dsDY_0b_w1b_elmu, dstt2l_0b_w1b_elmu, dstW_0b_w1b_elmu, dsttZ_0b_w1b_elmu, dsVV_0b_w1b_elmu});
    Plotter::DataCollection dcData_2015CD_nt2b_0belmuZinv(  "data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0belmuZinv_w2b});
    Plotter::DataCollection dcMC_nt2b_0belmuZinv(           "stack", "nTopCandSortedCntZinv2b", {dsDY_0belmuZinv_w2b, dstt2l_0belmuZinv_w2b, dstW_0belmuZinv_w2b, dsttZ_0belmuZinv_w2b, dsVV_0belmuZinv_w2b});
    Plotter::DataCollection dcData_2015CD_nt2b_0bnomt2_elmu("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0bnomt2_w2b_elmu});
    Plotter::DataCollection dcMC_nt2b_0bnomt2_elmu(         "stack", "nTopCandSortedCntZinv2b", {dsDY_0bnomt2_w2b_elmu, dstt2l_0bnomt2_w2b_elmu, dstW_0bnomt2_w2b_elmu, dsttZ_0bnomt2_w2b_elmu, dsVV_0bnomt2_w2b_elmu});
    Plotter::DataCollection dcData_2015CD_nt2b_0b_elmu(     "data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0b_w2b_elmu});
    Plotter::DataCollection dcMC_nt2b_0b_elmu(              "stack", "nTopCandSortedCntZinv2b", {dsDY_0b_w2b_elmu, dstt2l_0b_w2b_elmu, dstW_0b_w2b_elmu, dsttZ_0b_w2b_elmu, dsVV_0b_w2b_elmu});
    Plotter::DataCollection dcData_2015CD_nt3b_0belmuZinv(  "data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0belmuZinv_w3b});
    Plotter::DataCollection dcMC_nt3b_0belmuZinv(           "stack", "nTopCandSortedCntZinv3b", {dsDY_0belmuZinv_w3b, dstt2l_0belmuZinv_w3b, dstW_0belmuZinv_w3b, dsttZ_0belmuZinv_w3b, dsVV_0belmuZinv_w3b});
    Plotter::DataCollection dcData_2015CD_nt3b_0bnomt2_elmu("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0bnomt2_w3b_elmu});
    Plotter::DataCollection dcMC_nt3b_0bnomt2_elmu(         "stack", "nTopCandSortedCntZinv3b", {dsDY_0bnomt2_w3b_elmu, dstt2l_0bnomt2_w3b_elmu, dstW_0bnomt2_w3b_elmu, dsttZ_0bnomt2_w3b_elmu, dsVV_0bnomt2_w3b_elmu});
    Plotter::DataCollection dcData_2015CD_nt3b_0b_elmu(     "data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0b_w3b_elmu});
    Plotter::DataCollection dcMC_nt3b_0b_elmu(              "stack", "nTopCandSortedCntZinv3b", {dsDY_0b_w3b_elmu, dstt2l_0b_w3b_elmu, dstW_0b_w3b_elmu, dsttZ_0b_w3b_elmu, dsVV_0b_w3b_elmu});
    // --> MT2 with bfaking
    Plotter::DataCollection dcData_2015CD_mt21b_0belmuZinv(  "data",  "best_had_brJet_MT2Zinv1b", {dsData_2015C_0belmuZinv_w1b});
    Plotter::DataCollection dcMC_mt21b_0belmuZinv(           "stack", "best_had_brJet_MT2Zinv1b", {dsDY_0belmuZinv_w1b, dstt2l_0belmuZinv_w1b, dstW_0belmuZinv_w1b, dsttZ_0belmuZinv_w1b, dsVV_0belmuZinv_w1b});
    Plotter::DataCollection dcData_2015CD_mt21b_0bnomt2_elmu("data",  "best_had_brJet_MT2Zinv1b", {dsData_2015C_0bnomt2_w1b_elmu});
    Plotter::DataCollection dcMC_mt21b_0bnomt2_elmu(         "stack", "best_had_brJet_MT2Zinv1b", {dsDY_0bnomt2_w1b_elmu, dstt2l_0bnomt2_w1b_elmu, dstW_0bnomt2_w1b_elmu, dsttZ_0bnomt2_w1b_elmu, dsVV_0bnomt2_w1b_elmu});
    Plotter::DataCollection dcData_2015CD_mt21b_0b_elmu(     "data",  "best_had_brJet_MT2Zinv1b", {dsData_2015C_0b_w1b_elmu});
    Plotter::DataCollection dcMC_mt21b_0b_elmu(              "stack", "best_had_brJet_MT2Zinv1b", {dsDY_0b_w1b_elmu, dstt2l_0b_w1b_elmu, dstW_0b_w1b_elmu, dsttZ_0b_w1b_elmu, dsVV_0b_w1b_elmu});
    Plotter::DataCollection dcData_2015CD_mt22b_0belmuZinv(  "data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0belmuZinv_w2b});
    Plotter::DataCollection dcMC_mt22b_0belmuZinv(           "stack", "best_had_brJet_MT2Zinv2b", {dsDY_0belmuZinv_w2b, dstt2l_0belmuZinv_w2b, dstW_0belmuZinv_w2b, dsttZ_0belmuZinv_w2b, dsVV_0belmuZinv_w2b});
    Plotter::DataCollection dcData_2015CD_mt22b_0bnomt2_elmu("data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0bnomt2_w2b_elmu});
    Plotter::DataCollection dcMC_mt22b_0bnomt2_elmu(         "stack", "best_had_brJet_MT2Zinv2b", {dsDY_0bnomt2_w2b_elmu, dstt2l_0bnomt2_w2b_elmu, dstW_0bnomt2_w2b_elmu, dsttZ_0bnomt2_w2b_elmu, dsVV_0bnomt2_w2b_elmu});
    Plotter::DataCollection dcData_2015CD_mt22b_0b_elmu(     "data",  "best_had_brJet_MT2Zinv2b", {dsData_2015C_0b_w2b_elmu});
    Plotter::DataCollection dcMC_mt22b_0b_elmu(              "stack", "best_had_brJet_MT2Zinv2b", {dsDY_0b_w2b_elmu, dstt2l_0b_w2b_elmu, dstW_0b_w2b_elmu, dsttZ_0b_w2b_elmu, dsVV_0b_w2b_elmu});
    Plotter::DataCollection dcData_2015CD_mt23b_0belmuZinv(  "data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0belmuZinv_w3b});
    Plotter::DataCollection dcMC_mt23b_0belmuZinv(           "stack", "best_had_brJet_MT2Zinv3b", {dsDY_0belmuZinv_w3b, dstt2l_0belmuZinv_w3b, dstW_0belmuZinv_w3b, dsttZ_0belmuZinv_w3b, dsVV_0belmuZinv_w3b});
    Plotter::DataCollection dcData_2015CD_mt23b_0bnomt2_elmu("data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0bnomt2_w3b_elmu});
    Plotter::DataCollection dcMC_mt23b_0bnomt2_elmu(         "stack", "best_had_brJet_MT2Zinv3b", {dsDY_0bnomt2_w3b_elmu, dstt2l_0bnomt2_w3b_elmu, dstW_0bnomt2_w3b_elmu, dsttZ_0bnomt2_w3b_elmu, dsVV_0bnomt2_w3b_elmu});
    Plotter::DataCollection dcData_2015CD_mt23b_0b_elmu(     "data",  "best_had_brJet_MT2Zinv3b", {dsData_2015C_0b_w3b_elmu});
    Plotter::DataCollection dcMC_mt23b_0b_elmu(              "stack", "best_had_brJet_MT2Zinv3b", {dsDY_0b_w3b_elmu, dstt2l_0b_w3b_elmu, dstW_0b_w3b_elmu, dsttZ_0b_w3b_elmu, dsVV_0b_w3b_elmu});
    // --> nbs
    Plotter::DataCollection dcData_2015CD_nb_elmu(        "data",  "cntCSVSZinv", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_nb_elmu(                 "stack", "cntCSVSZinv", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_nb_elmuZinv(    "data",  "cntCSVSZinv", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_nb_elmuZinv(             "stack", "cntCSVSZinv", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_nb_ht200_elmu(  "data",  "cntCSVSZinv", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_nb_ht200_elmu(           "stack", "cntCSVSZinv", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_nb_blnotag_elmu("data",  "cntCSVSZinv", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_nb_blnotag_elmu(         "stack", "cntCSVSZinv", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_nb_bl_elmu(     "data",  "cntCSVSZinv", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_nb_bl_elmu(              "stack", "cntCSVSZinv", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_nb_0bnomt2_elmu("data",  "cntCSVSZinv", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_nb_0bnomt2_elmu(         "stack", "cntCSVSZinv", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_nb_0b_elmu(     "data",  "cntCSVSZinv", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_nb_0b_elmu(              "stack", "cntCSVSZinv", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> njets
    Plotter::DataCollection dcData_2015CD_nj_elmu(        "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_nj_elmu(                 "stack", "cntNJetsPt30Eta24Zinv", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_nj_elmuZinv(    "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_nj_elmuZinv(             "stack", "cntNJetsPt30Eta24Zinv", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_nj_ht200_elmu(  "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_nj_ht200_elmu(           "stack", "cntNJetsPt30Eta24Zinv", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_nj_blnotag_elmu("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_nj_blnotag_elmu(         "stack", "cntNJetsPt30Eta24Zinv", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_nj_bl_elmu(     "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_nj_bl_elmu(              "stack", "cntNJetsPt30Eta24Zinv", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_nj_0belmuZinv(  "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_nj_0belmuZinv(           "stack", "cntNJetsPt30Eta24Zinv", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_nj_0bnomt2_elmu("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_nj_0bnomt2_elmu(         "stack", "cntNJetsPt30Eta24Zinv", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_nj_0b_elmu(     "data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_nj_0b_elmu(              "stack", "cntNJetsPt30Eta24Zinv", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> HT
    Plotter::DataCollection dcData_2015CD_ht_elmu("data",  "HTZinv", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_ht_elmu("stack",  "HTZinv", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_ht_elmuZinv("data",  "HTZinv", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_ht_elmuZinv("stack",  "HTZinv", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_ht_ht200_elmu("data",  "HTZinv", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_ht_ht200_elmu("stack",  "HTZinv", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_ht_blnotag_elmu("data",  "HTZinv", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_ht_blnotag_elmu("stack",  "HTZinv", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_ht_bl_elmu("data",  "HTZinv", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_ht_bl_elmu("stack",  "HTZinv", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_ht_0belmuZinv("data",  "HTZinv", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_ht_0belmuZinv("stack",  "HTZinv", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_ht_0bnomt2_elmu("data",  "HTZinv", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_ht_0bnomt2_elmu("stack",  "HTZinv", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_ht_0b_elmu("data",  "HTZinv", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_ht_0b_elmu("stack",  "HTZinv", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> MHT
    Plotter::DataCollection dcData_2015CD_mht_elmu("data",  "cleanMHt", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_mht_elmu("stack",  "cleanMHt", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_mht_elmuZinv("data",  "cleanMHt", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_mht_elmuZinv("stack",  "cleanMHt", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_mht_ht200_elmu("data",  "cleanMHt", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_mht_ht200_elmu("stack",  "cleanMHt", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_mht_blnotag_elmu("data",  "cleanMHt", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_mht_blnotag_elmu("stack",  "cleanMHt", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_mht_bl_elmu("data",  "cleanMHt", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_mht_bl_elmu("stack",  "cleanMHt", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_mht_0belmuZinv("data",  "cleanMHt", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_mht_0belmuZinv("stack",  "cleanMHt", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_mht_0bnomt2_elmu("data",  "cleanMHt", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_mht_0bnomt2_elmu("stack",  "cleanMHt", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_mht_0b_elmu("data",  "cleanMHt", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_mht_0b_elmu("stack",  "cleanMHt", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> jet pt
    Plotter::DataCollection dcData_2015CD_jpt_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_jpt_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_jpt_elmuZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_jpt_elmuZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_jpt_ht200_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_jpt_ht200_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_jpt_blnotag_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_jpt_blnotag_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_jpt_bl_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_jpt_bl_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_jpt_0belmuZinv("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_jpt_0belmuZinv("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_jpt_0bnomt2_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_jpt_0bnomt2_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_jpt_0b_elmu("data",  "jetsLVecLepCleaned(pt)", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_jpt_0b_elmu("stack",  "jetsLVecLepCleaned(pt)", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> jet1 pt 
    Plotter::DataCollection dcData_2015CD_j1pt_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_j1pt_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_j1pt_elmuZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_j1pt_elmuZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_ht200_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_j1pt_ht200_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_j1pt_blnotag_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_j1pt_blnotag_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_j1pt_bl_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_j1pt_bl_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_j1pt_0belmuZinv("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_j1pt_0belmuZinv("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_j1pt_0bnomt2_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_j1pt_0bnomt2_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_j1pt_0b_elmu("data",  "jetsLVecLepCleaned[0](pt)", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_j1pt_0b_elmu("stack",  "jetsLVecLepCleaned[0](pt)", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> jet2 pt
    Plotter::DataCollection dcData_2015CD_j2pt_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_j2pt_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_j2pt_elmuZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_j2pt_elmuZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_ht200_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_j2pt_ht200_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_j2pt_blnotag_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_j2pt_blnotag_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_j2pt_bl_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_j2pt_bl_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_j2pt_0belmuZinv("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_j2pt_0belmuZinv("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_j2pt_0bnomt2_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_j2pt_0bnomt2_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_j2pt_0b_elmu("data",  "jetsLVecLepCleaned[1](pt)", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_j2pt_0b_elmu("stack",  "jetsLVecLepCleaned[1](pt)", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> jet3 pt
    Plotter::DataCollection dcData_2015CD_j3pt_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_j3pt_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_j3pt_elmuZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_j3pt_elmuZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_ht200_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_j3pt_ht200_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_j3pt_blnotag_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_j3pt_blnotag_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_j3pt_bl_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_j3pt_bl_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_j3pt_0belmuZinv("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_j3pt_0belmuZinv("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_j3pt_0bnomt2_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_j3pt_0bnomt2_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_j3pt_0b_elmu("data",  "jetsLVecLepCleaned[2](pt)", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_j3pt_0b_elmu("stack",  "jetsLVecLepCleaned[2](pt)", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> mu pt
    Plotter::DataCollection dcData_2015CD_mu1pt_elmu("data",  "cutMuPt1", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_mu1pt_elmu("stack",  "cutMuPt1", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_mu1pt_elmuZinv("data",  "cutMuPt1", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_mu1pt_elmuZinv("stack",  "cutMuPt1", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_mu1pt_ht200_elmu("data",  "cutMuPt1", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_mu1pt_ht200_elmu("stack",  "cutMuPt1", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_mu1pt_blnotag_elmu("data",  "cutMuPt1", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_mu1pt_blnotag_elmu("stack",  "cutMuPt1", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_mu1pt_bl_elmu("data",  "cutMuPt1", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_mu1pt_bl_elmu("stack",  "cutMuPt1", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_mu1pt_0belmuZinv("data",  "cutMuPt1", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_mu1pt_0belmuZinv("stack",  "cutMuPt1", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_mu1pt_0bnomt2_elmu("data",  "cutMuPt1", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_mu1pt_0bnomt2_elmu("stack",  "cutMuPt1", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_mu1pt_0b_elmu("data",  "cutMuPt1", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_mu1pt_0b_elmu("stack",  "cutMuPt1", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> el pt
    Plotter::DataCollection dcData_2015CD_el1pt_elmu("data",  "cutElecPt1", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_el1pt_elmu("stack",  "cutElecPt1", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_el1pt_elmuZinv("data",  "cutElecPt1", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_el1pt_elmuZinv("stack",  "cutElecPt1", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_el1pt_ht200_elmu("data",  "cutElecPt1", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_el1pt_ht200_elmu("stack",  "cutElecPt1", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_el1pt_blnotag_elmu("data",  "cutElecPt1", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_el1pt_blnotag_elmu("stack",  "cutElecPt1", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_el1pt_bl_elmu("data",  "cutElecPt1", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_el1pt_bl_elmu("stack",  "cutElecPt1", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_el1pt_0belmuZinv("data",  "cutElecPt1", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_el1pt_0belmuZinv("stack",  "cutElecPt1", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_el1pt_0bnomt2_elmu("data",  "cutElecPt1", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_el1pt_0bnomt2_elmu("stack",  "cutElecPt1", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_el1pt_0b_elmu("data",  "cutElecPt1", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_el1pt_0b_elmu("stack",  "cutElecPt1", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});
    // --> mll
    Plotter::DataCollection dcData_2015CD_mll_elmu("data",  "bestRecoZM", {dsData_2015C_elmu});
    Plotter::DataCollection dcMC_mll_elmu("stack",  "bestRecoZM", {dsDY_elmu, dstt2l_elmu, dstW_elmu, dsttZ_elmu, dsVV_elmu});
    Plotter::DataCollection dcData_2015CD_mll_elmuZinv("data",  "bestRecoZM", {dsData_2015C_elmuZinv});
    Plotter::DataCollection dcMC_mll_elmuZinv("stack",  "bestRecoZM", {dsDY_elmuZinv, dstt2l_elmuZinv, dstW_elmuZinv, dsttZ_elmuZinv, dsVV_elmuZinv});
    Plotter::DataCollection dcData_2015CD_mll_ht200_elmu("data",  "bestRecoZM", {dsData_2015C_ht200_elmu});
    Plotter::DataCollection dcMC_mll_ht200_elmu("stack",  "bestRecoZM", {dsDY_ht200_elmu, dstt2l_ht200_elmu, dstW_ht200_elmu, dsttZ_ht200_elmu, dsVV_ht200_elmu});
    Plotter::DataCollection dcData_2015CD_mll_blnotag_elmu("data",  "bestRecoZM", {dsData_2015C_blnotag_elmu});
    Plotter::DataCollection dcMC_mll_blnotag_elmu("stack",  "bestRecoZM", {dsDY_blnotag_elmu, dstt2l_blnotag_elmu, dstW_blnotag_elmu, dsttZ_blnotag_elmu, dsVV_blnotag_elmu});
    Plotter::DataCollection dcData_2015CD_mll_bl_elmu("data",  "bestRecoZM", {dsData_2015C_bl_elmu});
    Plotter::DataCollection dcMC_mll_bl_elmu("stack",  "bestRecoZM", {dsDY_bl_elmu, dstt2l_bl_elmu, dstW_bl_elmu, dsttZ_bl_elmu, dsVV_bl_elmu});
    Plotter::DataCollection dcData_2015CD_mll_0belmuZinv("data",  "bestRecoZM", {dsData_2015C_0belmuZinv});
    Plotter::DataCollection dcMC_mll_0belmuZinv("stack",  "bestRecoZM", {dsDY_0belmuZinv, dstt2l_0belmuZinv, dstW_0belmuZinv, dsttZ_0belmuZinv, dsVV_0belmuZinv});
    Plotter::DataCollection dcData_2015CD_mll_0bnomt2_elmu("data",  "bestRecoZM", {dsData_2015C_0bnomt2_elmu});
    Plotter::DataCollection dcMC_mll_0bnomt2_elmu("stack",  "bestRecoZM", {dsDY_0bnomt2_elmu, dstt2l_0bnomt2_elmu, dstW_0bnomt2_elmu, dsttZ_0bnomt2_elmu, dsVV_0bnomt2_elmu});
    Plotter::DataCollection dcData_2015CD_mll_0b_elmu("data",  "bestRecoZM", {dsData_2015C_0b_elmu});
    Plotter::DataCollection dcMC_mll_0b_elmu("stack",  "bestRecoZM", {dsDY_0b_elmu, dstt2l_0b_elmu, dstW_0b_elmu, dsttZ_0b_elmu, dsVV_0b_elmu});




    // Push back the histograms into the histo vector
    // --> nosel
    vh.push_back(PHS("DataMC_2015CD_met_nosel",   {dcData_2015CD_met_nosel, dcMC_met_nosel},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_nosel",    {dcData_2015CD_ht_nosel, dcMC_ht_nosel},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_nosel",    {dcData_2015CD_nt_nosel, dcMC_nt_nosel},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_nosel",   {dcData_2015CD_mt2_nosel, dcMC_mt2_nosel},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_nosel",    {dcData_2015CD_nb_nosel, dcMC_nb_nosel},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_nosel",    {dcData_2015CD_nj_nosel, dcMC_nj_nosel},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_nosel",   {dcData_2015CD_mht_nosel, dcMC_mht_nosel},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_nosel",   {dcData_2015CD_jpt_nosel, dcMC_jpt_nosel},      {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_nosel",  {dcData_2015CD_j1pt_nosel, dcMC_j1pt_nosel},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_nosel",  {dcData_2015CD_j2pt_nosel, dcMC_j2pt_nosel},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_nosel",  {dcData_2015CD_j3pt_nosel, dcMC_j3pt_nosel},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_nosel",  {dcData_2015CD_mupt_nosel, dcMC_mupt_nosel},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_nosel", {dcData_2015CD_mu1pt_nosel, dcMC_mu1pt_nosel}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_nosel", {dcData_2015CD_mu2pt_nosel, dcMC_mu2pt_nosel}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",        ""));
    // --> 2mu
    vh.push_back(PHS("DataMC_2015CD_met_2mu",   {dcData_2015CD_met_2mu, dcMC_met_2mu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_2mu",    {dcData_2015CD_ht_2mu, dcMC_ht_2mu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_2mu",    {dcData_2015CD_nt_2mu, dcMC_nt_2mu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_2mu",   {dcData_2015CD_mt2_2mu, dcMC_mt2_2mu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_2mu",    {dcData_2015CD_nb_2mu, dcMC_nb_2mu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_2mu",    {dcData_2015CD_nj_2mu, dcMC_nj_2mu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_2mu",   {dcData_2015CD_mht_2mu, dcMC_mht_2mu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_2mu",   {dcData_2015CD_jpt_2mu, dcMC_jpt_2mu},      {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_2mu",  {dcData_2015CD_j1pt_2mu, dcMC_j1pt_2mu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_2mu",  {dcData_2015CD_j2pt_2mu, dcMC_j2pt_2mu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_2mu",  {dcData_2015CD_j3pt_2mu, dcMC_j3pt_2mu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_2mu",  {dcData_2015CD_mupt_2mu, dcMC_mupt_2mu},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_2mu", {dcData_2015CD_mu1pt_2mu, dcMC_mu1pt_2mu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_2mu", {dcData_2015CD_mu2pt_2mu, dcMC_mu2pt_2mu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_2mu",   {dcData_2015CD_mll_2mu, dcMC_mll_2mu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> 2mu
    vh.push_back(PHS("DataMC_2015CD_met_muZinv",   {dcData_2015CD_met_muZinv, dcMC_met_muZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_muZinv",    {dcData_2015CD_ht_muZinv, dcMC_ht_muZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_muZinv",    {dcData_2015CD_nt_muZinv, dcMC_nt_muZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_muZinv",   {dcData_2015CD_mt2_muZinv, dcMC_mt2_muZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_muZinv",    {dcData_2015CD_nb_muZinv, dcMC_nb_muZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_muZinv",    {dcData_2015CD_nj_muZinv, dcMC_nj_muZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_muZinv",   {dcData_2015CD_mht_muZinv, dcMC_mht_muZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_muZinv",   {dcData_2015CD_jpt_muZinv, dcMC_jpt_muZinv},      {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_muZinv",  {dcData_2015CD_j1pt_muZinv, dcMC_j1pt_muZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_muZinv",  {dcData_2015CD_j2pt_muZinv, dcMC_j2pt_muZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_muZinv",  {dcData_2015CD_j3pt_muZinv, dcMC_j3pt_muZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_muZinv",  {dcData_2015CD_mupt_muZinv, dcMC_mupt_muZinv},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_muZinv", {dcData_2015CD_mu1pt_muZinv, dcMC_mu1pt_muZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_muZinv", {dcData_2015CD_mu2pt_muZinv, dcMC_mu2pt_muZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_muZinv",   {dcData_2015CD_mll_muZinv, dcMC_mll_muZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> 2mu
    vh.push_back(PHS("DataMC_2015CD_met_ht200",   {dcData_2015CD_met_ht200, dcMC_met_ht200},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_ht200",    {dcData_2015CD_ht_ht200, dcMC_ht_ht200},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_ht200",    {dcData_2015CD_nt_ht200, dcMC_nt_ht200},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_ht200",   {dcData_2015CD_mt2_ht200, dcMC_mt2_ht200},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_ht200",    {dcData_2015CD_nb_ht200, dcMC_nb_ht200},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_ht200",    {dcData_2015CD_nj_ht200, dcMC_nj_ht200},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_ht200",   {dcData_2015CD_mht_ht200, dcMC_mht_ht200},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_ht200",   {dcData_2015CD_jpt_ht200, dcMC_jpt_ht200},      {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_ht200",  {dcData_2015CD_j1pt_ht200, dcMC_j1pt_ht200},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_ht200",  {dcData_2015CD_j2pt_ht200, dcMC_j2pt_ht200},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_ht200",  {dcData_2015CD_j3pt_ht200, dcMC_j3pt_ht200},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_ht200",  {dcData_2015CD_mupt_ht200, dcMC_mupt_ht200},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_ht200", {dcData_2015CD_mu1pt_ht200, dcMC_mu1pt_ht200}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_ht200", {dcData_2015CD_mu2pt_ht200, dcMC_mu2pt_ht200}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_ht200",   {dcData_2015CD_mll_ht200, dcMC_mll_ht200},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> baselineNoTag
    vh.push_back(PHS("DataMC_2015CD_met_baselineNoTag",   {dcData_2015CD_met_blnotag, dcMC_met_blnotag},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baselineNoTag",    {dcData_2015CD_ht_blnotag, dcMC_ht_blnotag},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baselineNoTag",    {dcData_2015CD_nt_blnotag, dcMC_nt_blnotag},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baselineNoTag",   {dcData_2015CD_mt2_blnotag, dcMC_mt2_blnotag},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baselineNoTag",    {dcData_2015CD_nb_blnotag, dcMC_nb_blnotag},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baselineNoTag",    {dcData_2015CD_nj_blnotag, dcMC_nj_blnotag},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baselineNoTag",   {dcData_2015CD_mht_blnotag, dcMC_mht_blnotag},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baselineNoTag",   {dcData_2015CD_jpt_blnotag, dcMC_jpt_blnotag},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baselineNoTag",  {dcData_2015CD_j1pt_blnotag, dcMC_j1pt_blnotag},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baselineNoTag",  {dcData_2015CD_j2pt_blnotag, dcMC_j2pt_blnotag},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baselineNoTag",  {dcData_2015CD_j3pt_blnotag, dcMC_j3pt_blnotag},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baselineNoTag",  {dcData_2015CD_mupt_blnotag, dcMC_mupt_blnotag},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baselineNoTag", {dcData_2015CD_mu1pt_blnotag, dcMC_mu1pt_blnotag}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baselineNoTag", {dcData_2015CD_mu2pt_blnotag, dcMC_mu2pt_blnotag}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baselineNoTag",   {dcData_2015CD_mll_blnotag, dcMC_mll_blnotag},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> baseline
    vh.push_back(PHS("DataMC_2015CD_met_baseline",   {dcData_2015CD_met_bl, dcMC_met_bl},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",            ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baseline",    {dcData_2015CD_ht_bl, dcMC_ht_bl},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baseline",    {dcData_2015CD_nt_bl, dcMC_nt_bl},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",           ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baseline",   {dcData_2015CD_mt2_bl, dcMC_mt2_bl},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",            ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baseline",    {dcData_2015CD_nb_bl, dcMC_nb_bl},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",             ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baseline",    {dcData_2015CD_nj_bl, dcMC_nj_bl},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",             ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baseline",   {dcData_2015CD_mht_bl, dcMC_mht_bl},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",            ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baseline",   {dcData_2015CD_jpt_bl, dcMC_jpt_bl},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baseline",  {dcData_2015CD_j1pt_bl, dcMC_j1pt_bl},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baseline",  {dcData_2015CD_j2pt_bl, dcMC_j2pt_bl},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baseline",  {dcData_2015CD_j3pt_bl, dcMC_j3pt_bl},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baseline",  {dcData_2015CD_mupt_bl, dcMC_mupt_bl},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baseline", {dcData_2015CD_mu1pt_bl, dcMC_mu1pt_bl}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baseline", {dcData_2015CD_mu2pt_bl, dcMC_mu2pt_bl}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baseline",   {dcData_2015CD_mll_bl, dcMC_mll_bl},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",            ""));
    // --> 0b muZinv
    vh.push_back(PHS("DataMC_2015CD_met_0bmuZinv",   {dcData_2015CD_met_0bmuZinv, dcMC_met_0bmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0bmuZinv",    {dcData_2015CD_ht_0bmuZinv, dcMC_ht_0bmuZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0bmuZinv",    {dcData_2015CD_nt_0bmuZinv, dcMC_nt_0bmuZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0bmuZinv",  {dcData_2015CD_nt1b_0bmuZinv, dcMC_nt1b_0bmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0bmuZinv",  {dcData_2015CD_nt2b_0bmuZinv, dcMC_nt2b_0bmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0bmuZinv",  {dcData_2015CD_nt3b_0bmuZinv, dcMC_nt3b_0bmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0bmuZinv",   {dcData_2015CD_mt2_0bmuZinv, dcMC_mt2_0bmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0bmuZinv", {dcData_2015CD_mt21b_0bmuZinv, dcMC_mt21b_0bmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0bmuZinv", {dcData_2015CD_mt22b_0bmuZinv, dcMC_mt22b_0bmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0bmuZinv", {dcData_2015CD_mt23b_0bmuZinv, dcMC_mt23b_0bmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    //vh.push_back(PHS("DataMC_2015CD_nb_0bmuZinv",    {dcData_2015CD_nb_0bmuZinv, dcMC_nb_0bmuZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0bmuZinv",    {dcData_2015CD_nj_0bmuZinv, dcMC_nj_0bmuZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0bmuZinv",   {dcData_2015CD_mht_0bmuZinv, dcMC_mht_0bmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0bmuZinv",   {dcData_2015CD_jpt_0bmuZinv, dcMC_jpt_0bmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0bmuZinv",  {dcData_2015CD_j1pt_0bmuZinv, dcMC_j1pt_0bmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0bmuZinv",  {dcData_2015CD_j2pt_0bmuZinv, dcMC_j2pt_0bmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0bmuZinv",  {dcData_2015CD_j3pt_0bmuZinv, dcMC_j3pt_0bmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_0bmuZinv",  {dcData_2015CD_mupt_0bmuZinv, dcMC_mupt_0bmuZinv},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",            ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0bmuZinv", {dcData_2015CD_mu1pt_0bmuZinv, dcMC_mu1pt_0bmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_0bmuZinv", {dcData_2015CD_mu2pt_0bmuZinv, dcMC_mu2pt_0bmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0bmuZinv",   {dcData_2015CD_mll_0bmuZinv, dcMC_mll_0bmuZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));
    // --> 0b no mt2
    vh.push_back(PHS("DataMC_2015CD_met_0bnomt2",   {dcData_2015CD_met_0bnomt2, dcMC_met_0bnomt2},     {1, 2}, "", 50, 0, 1500,  true, false,  "met",             ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0bnomt2",    {dcData_2015CD_ht_0bnomt2, dcMC_ht_0bnomt2},       {1, 2}, "", 50, 0, 1500,  true, false,  "ht",              ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0bnomt2",    {dcData_2015CD_nt_0bnomt2, dcMC_nt_0bnomt2},       {1, 2}, "", 5, 0, 5,      true, false,  "ntop",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0bnomt2",  {dcData_2015CD_nt1b_0bnomt2, dcMC_nt1b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0bnomt2",  {dcData_2015CD_nt2b_0bnomt2, dcMC_nt2b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0bnomt2",  {dcData_2015CD_nt3b_0bnomt2, dcMC_nt3b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0bnomt2",   {dcData_2015CD_mt2_0bnomt2, dcMC_mt2_0bnomt2},     {1, 2}, "", 50, 0, 1500,  true, false,  "mt2",             ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0bnomt2", {dcData_2015CD_mt21b_0bnomt2, dcMC_mt21b_0bnomt2}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0bnomt2", {dcData_2015CD_mt22b_0bnomt2, dcMC_mt22b_0bnomt2}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0bnomt2", {dcData_2015CD_mt23b_0bnomt2, dcMC_mt23b_0bnomt2}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0bnomt2",    {dcData_2015CD_nb_0bnomt2, dcMC_nb_0bnomt2},       {1, 2}, "", 10, 0, 10,    true, false,  "nb",              ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0bnomt2",    {dcData_2015CD_nj_0bnomt2, dcMC_nj_0bnomt2},       {1, 2}, "", 20, 0, 20,    true, false,  "nj",              ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0bnomt2",   {dcData_2015CD_mht_0bnomt2, dcMC_mht_0bnomt2},     {1, 2}, "", 50, 0, 1500,  true, false,  "mht",             ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0bnomt2",   {dcData_2015CD_jpt_0bnomt2, dcMC_jpt_0bnomt2},     {1, 2}, "", 50, 0, 1500,  true, false,  "jet pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0bnomt2",  {dcData_2015CD_j1pt_0bnomt2, dcMC_j1pt_0bnomt2},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0bnomt2",  {dcData_2015CD_j2pt_0bnomt2, dcMC_j2pt_0bnomt2},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0bnomt2",  {dcData_2015CD_j3pt_0bnomt2, dcMC_j3pt_0bnomt2},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_0bnomt2",  {dcData_2015CD_mupt_0bnomt2, dcMC_mupt_0bnomt2},   {1, 2}, "", 50, 0, 1000,  true, false,  "mu pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0bnomt2", {dcData_2015CD_mu1pt_0bnomt2, dcMC_mu1pt_0bnomt2}, {1, 2}, "", 50, 0, 1000,  true, false,  "mu1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_0bnomt2", {dcData_2015CD_mu2pt_0bnomt2, dcMC_mu2pt_0bnomt2}, {1, 2}, "", 50, 0, 1000,  true, false,  "mu2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0bnomt2",   {dcData_2015CD_mll_0bnomt2, dcMC_mll_0bnomt2},     {1, 2}, "", 50, 0, 1000,  true, false,  "mll",             ""));
    // --> 0b
    vh.push_back(PHS("DataMC_2015CD_met_0b",   {dcData_2015CD_met_0b, dcMC_met_0b},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0b",    {dcData_2015CD_ht_0b, dcMC_ht_0b},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0b",    {dcData_2015CD_nt_0b, dcMC_nt_0b},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0b",  {dcData_2015CD_nt1b_0b, dcMC_nt1b_0b},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0b",  {dcData_2015CD_nt2b_0b, dcMC_nt2b_0b},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0b",  {dcData_2015CD_nt3b_0b, dcMC_nt3b_0b},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0b",   {dcData_2015CD_mt2_0b, dcMC_mt2_0b},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0b", {dcData_2015CD_mt21b_0b, dcMC_mt21b_0b}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0b", {dcData_2015CD_mt22b_0b, dcMC_mt22b_0b}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0b", {dcData_2015CD_mt23b_0b, dcMC_mt23b_0b}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0b",    {dcData_2015CD_nb_0b, dcMC_nb_0b},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0b",    {dcData_2015CD_nj_0b, dcMC_nj_0b},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0b",   {dcData_2015CD_mht_0b, dcMC_mht_0b},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0b",   {dcData_2015CD_jpt_0b, dcMC_jpt_0b},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0b",  {dcData_2015CD_j1pt_0b, dcMC_j1pt_0b},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0b",  {dcData_2015CD_j2pt_0b, dcMC_j2pt_0b},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0b",  {dcData_2015CD_j3pt_0b, dcMC_j3pt_0b},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_0b",  {dcData_2015CD_mupt_0b, dcMC_mupt_0b},   {1, 2}, "", 50, 0, 1000,   true, false,  "mu pt",            ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0b", {dcData_2015CD_mu1pt_0b, dcMC_mu1pt_0b}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_0b", {dcData_2015CD_mu2pt_0b, dcMC_mu2pt_0b}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu2 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0b",   {dcData_2015CD_mll_0b, dcMC_mll_0b},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));


    // electrons
    // --> 2el
    vh.push_back(PHS("DataMC_2015CD_met_2el",   {dcData_2015CD_met_2el, dcMC_met_2el},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_2el",    {dcData_2015CD_ht_2el, dcMC_ht_2el},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_2el",    {dcData_2015CD_nt_2el, dcMC_nt_2el},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_2el",   {dcData_2015CD_mt2_2el, dcMC_mt2_2el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_2el",    {dcData_2015CD_nb_2el, dcMC_nb_2el},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_2el",    {dcData_2015CD_nj_2el, dcMC_nj_2el},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_2el",   {dcData_2015CD_mht_2el, dcMC_mht_2el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_2el",   {dcData_2015CD_jpt_2el, dcMC_jpt_2el},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_2el",  {dcData_2015CD_j1pt_2el, dcMC_j1pt_2el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_2el",  {dcData_2015CD_j2pt_2el, dcMC_j2pt_2el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_2el",  {dcData_2015CD_j3pt_2el, dcMC_j3pt_2el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_2el",  {dcData_2015CD_elpt_2el, dcMC_elpt_2el},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_2el", {dcData_2015CD_el1pt_2el, dcMC_el1pt_2el}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_2el", {dcData_2015CD_el2pt_2el, dcMC_el2pt_2el}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_2el",   {dcData_2015CD_mll_2el, dcMC_mll_2el},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> 2mu
    vh.push_back(PHS("DataMC_2015CD_met_elZinv",   {dcData_2015CD_met_elZinv, dcMC_met_elZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_elZinv",    {dcData_2015CD_ht_elZinv, dcMC_ht_elZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_elZinv",    {dcData_2015CD_nt_elZinv, dcMC_nt_elZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_elZinv",   {dcData_2015CD_mt2_elZinv, dcMC_mt2_elZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_elZinv",    {dcData_2015CD_nb_elZinv, dcMC_nb_elZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_elZinv",    {dcData_2015CD_nj_elZinv, dcMC_nj_elZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_elZinv",   {dcData_2015CD_mht_elZinv, dcMC_mht_elZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_elZinv",   {dcData_2015CD_jpt_elZinv, dcMC_jpt_elZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_elZinv",  {dcData_2015CD_j1pt_elZinv, dcMC_j1pt_elZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_elZinv",  {dcData_2015CD_j2pt_elZinv, dcMC_j2pt_elZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_elZinv",  {dcData_2015CD_j3pt_elZinv, dcMC_j3pt_elZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_elZinv",  {dcData_2015CD_elpt_elZinv, dcMC_elpt_elZinv},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_elZinv", {dcData_2015CD_el1pt_elZinv, dcMC_el1pt_elZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_elZinv", {dcData_2015CD_el2pt_elZinv, dcMC_el2pt_elZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_elZinv",   {dcData_2015CD_mll_elZinv, dcMC_mll_elZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> 2mu
    vh.push_back(PHS("DataMC_2015CD_met_ht200_el",   {dcData_2015CD_met_ht200_el, dcMC_met_ht200_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_ht200_el",    {dcData_2015CD_ht_ht200_el, dcMC_ht_ht200_el},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_ht200_el",    {dcData_2015CD_nt_ht200_el, dcMC_nt_ht200_el},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_ht200_el",   {dcData_2015CD_mt2_ht200_el, dcMC_mt2_ht200_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_ht200_el",    {dcData_2015CD_nb_ht200_el, dcMC_nb_ht200_el},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_ht200_el",    {dcData_2015CD_nj_ht200_el, dcMC_nj_ht200_el},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_ht200_el",   {dcData_2015CD_mht_ht200_el, dcMC_mht_ht200_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_ht200_el",   {dcData_2015CD_jpt_ht200_el, dcMC_jpt_ht200_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_ht200_el",  {dcData_2015CD_j1pt_ht200_el, dcMC_j1pt_ht200_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_ht200_el",  {dcData_2015CD_j2pt_ht200_el, dcMC_j2pt_ht200_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_ht200_el",  {dcData_2015CD_j3pt_ht200_el, dcMC_j3pt_ht200_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_ht200",     {dcData_2015CD_elpt_ht200, dcMC_elpt_ht200},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_ht200",    {dcData_2015CD_el1pt_ht200, dcMC_el1pt_ht200}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_ht200",    {dcData_2015CD_el2pt_ht200, dcMC_el2pt_ht200}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_ht200_el",   {dcData_2015CD_mll_ht200_el, dcMC_mll_ht200_el},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> baselineNoTag
    vh.push_back(PHS("DataMC_2015CD_met_baselineNoTag_el",   {dcData_2015CD_met_blnotag_el, dcMC_met_blnotag_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",        ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baselineNoTag_el",    {dcData_2015CD_ht_blnotag_el, dcMC_ht_blnotag_el},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baselineNoTag_el",    {dcData_2015CD_nt_blnotag_el, dcMC_nt_blnotag_el},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",       ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baselineNoTag_el",   {dcData_2015CD_mt2_blnotag_el, dcMC_mt2_blnotag_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",        ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baselineNoTag_el",    {dcData_2015CD_nb_blnotag_el, dcMC_nb_blnotag_el},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baselineNoTag_el",    {dcData_2015CD_nj_blnotag_el, dcMC_nj_blnotag_el},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baselineNoTag_el",   {dcData_2015CD_mht_blnotag_el, dcMC_mht_blnotag_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",        ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baselineNoTag_el",   {dcData_2015CD_jpt_blnotag_el, dcMC_jpt_blnotag_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baselineNoTag_el",  {dcData_2015CD_j1pt_blnotag_el, dcMC_j1pt_blnotag_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baselineNoTag_el",  {dcData_2015CD_j2pt_blnotag_el, dcMC_j2pt_blnotag_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baselineNoTag_el",  {dcData_2015CD_j3pt_blnotag_el, dcMC_j3pt_blnotag_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_baselineNoTag",     {dcData_2015CD_elpt_blnotag, dcMC_elpt_blnotag},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",      ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_baselineNoTag",    {dcData_2015CD_el1pt_blnotag, dcMC_el1pt_blnotag}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_baselineNoTag",    {dcData_2015CD_el2pt_blnotag, dcMC_el2pt_blnotag}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baselineNoTag_el",   {dcData_2015CD_mll_blnotag_el, dcMC_mll_blnotag_el},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",        ""));
    // --> baseline
    vh.push_back(PHS("DataMC_2015CD_met_baseline_el",   {dcData_2015CD_met_bl_el, dcMC_met_bl_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",            ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baseline_el",    {dcData_2015CD_ht_bl_el, dcMC_ht_bl_el},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baseline_el",    {dcData_2015CD_nt_bl_el, dcMC_nt_bl_el},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",           ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baseline_el",   {dcData_2015CD_mt2_bl_el, dcMC_mt2_bl_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",            ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baseline_el",    {dcData_2015CD_nb_bl_el, dcMC_nb_bl_el},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",             ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baseline_el",    {dcData_2015CD_nj_bl_el, dcMC_nj_bl_el},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",             ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baseline_el",   {dcData_2015CD_mht_bl_el, dcMC_mht_bl_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",            ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baseline_el",   {dcData_2015CD_jpt_bl_el, dcMC_jpt_bl_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baseline_el",  {dcData_2015CD_j1pt_bl_el, dcMC_j1pt_bl_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baseline_el",  {dcData_2015CD_j2pt_bl_el, dcMC_j2pt_bl_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baseline_el",  {dcData_2015CD_j3pt_bl_el, dcMC_j3pt_bl_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_baseline",  {dcData_2015CD_elpt_bl, dcMC_elpt_bl},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_baseline", {dcData_2015CD_el1pt_bl, dcMC_el1pt_bl}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_baseline", {dcData_2015CD_el2pt_bl, dcMC_el2pt_bl}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baseline_el",   {dcData_2015CD_mll_bl_el, dcMC_mll_bl_el},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",            ""));
    // --> 0b muZinv
    vh.push_back(PHS("DataMC_2015CD_met_0belZinv",   {dcData_2015CD_met_0belZinv, dcMC_met_0belZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0belZinv",    {dcData_2015CD_ht_0belZinv, dcMC_ht_0belZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0belZinv",    {dcData_2015CD_nt_0belZinv, dcMC_nt_0belZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0belZinv",  {dcData_2015CD_nt1b_0belZinv, dcMC_nt1b_0belZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0belZinv",  {dcData_2015CD_nt2b_0belZinv, dcMC_nt2b_0belZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0belZinv",  {dcData_2015CD_nt3b_0belZinv, dcMC_nt3b_0belZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0belZinv",   {dcData_2015CD_mt2_0belZinv, dcMC_mt2_0belZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0belZinv", {dcData_2015CD_mt21b_0belZinv, dcMC_mt21b_0belZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0belZinv", {dcData_2015CD_mt22b_0belZinv, dcMC_mt22b_0belZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0belZinv", {dcData_2015CD_mt23b_0belZinv, dcMC_mt23b_0belZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    //vh.push_back(PHS("DataMC_2015CD_nb_0belZinv",    {dcData_2015CD_nb_0belZinv, dcMC_nb_0belZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0belZinv",    {dcData_2015CD_nj_0belZinv, dcMC_nj_0belZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0belZinv",   {dcData_2015CD_mht_0belZinv, dcMC_mht_0belZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0belZinv",   {dcData_2015CD_jpt_0belZinv, dcMC_jpt_0belZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0belZinv",  {dcData_2015CD_j1pt_0belZinv, dcMC_j1pt_0belZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0belZinv",  {dcData_2015CD_j2pt_0belZinv, dcMC_j2pt_0belZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0belZinv",  {dcData_2015CD_j3pt_0belZinv, dcMC_j3pt_0belZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_0belZinv",  {dcData_2015CD_elpt_0belZinv, dcMC_elpt_0belZinv},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",            ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0belZinv", {dcData_2015CD_el1pt_0belZinv, dcMC_el1pt_0belZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_0belZinv", {dcData_2015CD_el2pt_0belZinv, dcMC_el2pt_0belZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0belZinv",   {dcData_2015CD_mll_0belZinv, dcMC_mll_0belZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));
    // --> 0b no mt2
    vh.push_back(PHS("DataMC_2015CD_met_0bnomt2_el",   {dcData_2015CD_met_0bnomt2_el, dcMC_met_0bnomt2_el},     {1, 2}, "", 50, 0, 1500,  true, false,  "met",             ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0bnomt2_el",    {dcData_2015CD_ht_0bnomt2_el, dcMC_ht_0bnomt2_el},       {1, 2}, "", 50, 0, 1500,  true, false,  "ht",              ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0bnomt2_el",    {dcData_2015CD_nt_0bnomt2_el, dcMC_nt_0bnomt2_el},       {1, 2}, "", 5, 0, 5,      true, false,  "ntop",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0bnomt2_el",  {dcData_2015CD_nt1b_0bnomt2_el, dcMC_nt1b_0bnomt2_el},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0bnomt2_el",  {dcData_2015CD_nt2b_0bnomt2_el, dcMC_nt2b_0bnomt2_el},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0bnomt2_el",  {dcData_2015CD_nt3b_0bnomt2_el, dcMC_nt3b_0bnomt2_el},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0bnomt2_el",   {dcData_2015CD_mt2_0bnomt2_el, dcMC_mt2_0bnomt2_el},     {1, 2}, "", 50, 0, 1500,  true, false,  "mt2",             ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0bnomt2_el", {dcData_2015CD_mt21b_0bnomt2_el, dcMC_mt21b_0bnomt2_el}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0bnomt2_el", {dcData_2015CD_mt22b_0bnomt2_el, dcMC_mt22b_0bnomt2_el}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0bnomt2_el", {dcData_2015CD_mt23b_0bnomt2_el, dcMC_mt23b_0bnomt2_el}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0bnomt2_el",    {dcData_2015CD_nb_0bnomt2_el, dcMC_nb_0bnomt2_el},       {1, 2}, "", 10, 0, 10,    true, false,  "nb",              ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0bnomt2_el",    {dcData_2015CD_nj_0bnomt2_el, dcMC_nj_0bnomt2_el},       {1, 2}, "", 20, 0, 20,    true, false,  "nj",              ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0bnomt2_el",   {dcData_2015CD_mht_0bnomt2_el, dcMC_mht_0bnomt2_el},     {1, 2}, "", 50, 0, 1500,  true, false,  "mht",             ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0bnomt2_el",   {dcData_2015CD_jpt_0bnomt2_el, dcMC_jpt_0bnomt2_el},     {1, 2}, "", 50, 0, 1500,  true, false,  "jet pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0bnomt2_el",  {dcData_2015CD_j1pt_0bnomt2_el, dcMC_j1pt_0bnomt2_el},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0bnomt2_el",  {dcData_2015CD_j2pt_0bnomt2_el, dcMC_j2pt_0bnomt2_el},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0bnomt2_el",  {dcData_2015CD_j3pt_0bnomt2_el, dcMC_j3pt_0bnomt2_el},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_0bnomt2",     {dcData_2015CD_elpt_0bnomt2, dcMC_elpt_0bnomt2},         {1, 2}, "", 50, 0, 1000,  true, false,  "el pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0bnomt2",    {dcData_2015CD_el1pt_0bnomt2, dcMC_el1pt_0bnomt2},       {1, 2}, "", 50, 0, 1000,  true, false,  "el1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_0bnomt2",    {dcData_2015CD_el2pt_0bnomt2, dcMC_el2pt_0bnomt2},       {1, 2}, "", 50, 0, 1000,  true, false,  "el2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0bnomt2_el",   {dcData_2015CD_mll_0bnomt2_el, dcMC_mll_0bnomt2_el},     {1, 2}, "", 50, 0, 1000,  true, false,  "mll",             ""));
    // --> 0b
    vh.push_back(PHS("DataMC_2015CD_met_0b_el",   {dcData_2015CD_met_0b_el, dcMC_met_0b_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0b_el",    {dcData_2015CD_ht_0b_el, dcMC_ht_0b_el},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0b_el",    {dcData_2015CD_nt_0b_el, dcMC_nt_0b_el},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0b_el",  {dcData_2015CD_nt1b_0b_el, dcMC_nt1b_0b_el},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0b_el",  {dcData_2015CD_nt2b_0b_el, dcMC_nt2b_0b_el},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0b_el",  {dcData_2015CD_nt3b_0b_el, dcMC_nt3b_0b_el},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0b_el",   {dcData_2015CD_mt2_0b_el, dcMC_mt2_0b_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0b_el", {dcData_2015CD_mt21b_0b_el, dcMC_mt21b_0b_el}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0b_el", {dcData_2015CD_mt22b_0b_el, dcMC_mt22b_0b_el}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0b_el", {dcData_2015CD_mt23b_0b_el, dcMC_mt23b_0b_el}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0b_el",    {dcData_2015CD_nb_0b_el, dcMC_nb_0b_el},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0b_el",    {dcData_2015CD_nj_0b_el, dcMC_nj_0b_el},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0b_el",   {dcData_2015CD_mht_0b_el, dcMC_mht_0b_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0b_el",   {dcData_2015CD_jpt_0b_el, dcMC_jpt_0b_el},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0b_el",  {dcData_2015CD_j1pt_0b_el, dcMC_j1pt_0b_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0b_el",  {dcData_2015CD_j2pt_0b_el, dcMC_j2pt_0b_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0b_el",  {dcData_2015CD_j3pt_0b_el, dcMC_j3pt_0b_el},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_elpt_0b",     {dcData_2015CD_elpt_0b, dcMC_elpt_0b},   {1, 2}, "", 50, 0, 1000,   true, false,  "el pt",            ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0b",    {dcData_2015CD_el1pt_0b, dcMC_el1pt_0b}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_el2pt_0b",    {dcData_2015CD_el2pt_0b, dcMC_el2pt_0b}, {1, 2}, "", 50, 0, 1000,   true, false,  "el2 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0b_el",   {dcData_2015CD_mll_0b_el, dcMC_mll_0b_el},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));

    // el + mu
    // --> elmu
    vh.push_back(PHS("DataMC_2015CD_met_elmu",   {dcData_2015CD_met_elmu, dcMC_met_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_elmu",    {dcData_2015CD_ht_elmu, dcMC_ht_elmu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_elmu",    {dcData_2015CD_nt_elmu, dcMC_nt_elmu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_elmu",   {dcData_2015CD_mt2_elmu, dcMC_mt2_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_elmu",    {dcData_2015CD_nb_elmu, dcMC_nb_elmu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_elmu",    {dcData_2015CD_nj_elmu, dcMC_nj_elmu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_elmu",   {dcData_2015CD_mht_elmu, dcMC_mht_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_elmu",   {dcData_2015CD_jpt_elmu, dcMC_jpt_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_elmu",  {dcData_2015CD_j1pt_elmu, dcMC_j1pt_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_elmu",  {dcData_2015CD_j2pt_elmu, dcMC_j2pt_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_elmu",  {dcData_2015CD_j3pt_elmu, dcMC_j3pt_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_elmu", {dcData_2015CD_mu1pt_elmu, dcMC_mu1pt_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_elmu", {dcData_2015CD_el1pt_elmu, dcMC_el1pt_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_elmu",   {dcData_2015CD_mll_elmu, dcMC_mll_elmu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> elmu Zinv
    vh.push_back(PHS("DataMC_2015CD_met_elmuZinv",   {dcData_2015CD_met_elmuZinv, dcMC_met_elmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_elmuZinv",    {dcData_2015CD_ht_elmuZinv, dcMC_ht_elmuZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_elmuZinv",    {dcData_2015CD_nt_elmuZinv, dcMC_nt_elmuZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_elmuZinv",   {dcData_2015CD_mt2_elmuZinv, dcMC_mt2_elmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_elmuZinv",    {dcData_2015CD_nb_elmuZinv, dcMC_nb_elmuZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_elmuZinv",    {dcData_2015CD_nj_elmuZinv, dcMC_nj_elmuZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_elmuZinv",   {dcData_2015CD_mht_elmuZinv, dcMC_mht_elmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_elmuZinv",   {dcData_2015CD_jpt_elmuZinv, dcMC_jpt_elmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_elmuZinv",  {dcData_2015CD_j1pt_elmuZinv, dcMC_j1pt_elmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_elmuZinv",  {dcData_2015CD_j2pt_elmuZinv, dcMC_j2pt_elmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_elmuZinv",  {dcData_2015CD_j3pt_elmuZinv, dcMC_j3pt_elmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_elmuZinv", {dcData_2015CD_mu1pt_elmuZinv, dcMC_mu1pt_elmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_elmuZinv", {dcData_2015CD_el1pt_elmuZinv, dcMC_el1pt_elmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_elmuZinv",   {dcData_2015CD_mll_elmuZinv, dcMC_mll_elmuZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> ht200
    vh.push_back(PHS("DataMC_2015CD_met_ht200_elmu",   {dcData_2015CD_met_ht200_elmu, dcMC_met_ht200_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",           ""));
    vh.push_back(PHS("DataMC_2015CD_ht_ht200_elmu",    {dcData_2015CD_ht_ht200_elmu, dcMC_ht_ht200_elmu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt_ht200_elmu",    {dcData_2015CD_nt_ht200_elmu, dcMC_nt_ht200_elmu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",          ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_ht200_elmu",   {dcData_2015CD_mt2_ht200_elmu, dcMC_mt2_ht200_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",           ""));
    vh.push_back(PHS("DataMC_2015CD_nb_ht200_elmu",    {dcData_2015CD_nb_ht200_elmu, dcMC_nb_ht200_elmu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",            ""));
    vh.push_back(PHS("DataMC_2015CD_nj_ht200_elmu",    {dcData_2015CD_nj_ht200_elmu, dcMC_nj_ht200_elmu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",            ""));
    vh.push_back(PHS("DataMC_2015CD_mht_ht200_elmu",   {dcData_2015CD_mht_ht200_elmu, dcMC_mht_ht200_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",           ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_ht200_elmu",   {dcData_2015CD_jpt_ht200_elmu, dcMC_jpt_ht200_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_ht200_elmu",  {dcData_2015CD_j1pt_ht200_elmu, dcMC_j1pt_ht200_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_ht200_elmu",  {dcData_2015CD_j2pt_ht200_elmu, dcMC_j2pt_ht200_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_ht200_elmu",  {dcData_2015CD_j3pt_ht200_elmu, dcMC_j3pt_ht200_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",       ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_ht200_elmu", {dcData_2015CD_mu1pt_ht200_elmu, dcMC_mu1pt_ht200_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_ht200_elmu", {dcData_2015CD_el1pt_ht200_elmu, dcMC_el1pt_ht200_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mll_ht200_elmu",   {dcData_2015CD_mll_ht200_elmu, dcMC_mll_ht200_elmu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",           ""));
    // --> baselineNoTag
    vh.push_back(PHS("DataMC_2015CD_met_baselineNoTag_elmu",   {dcData_2015CD_met_blnotag_elmu, dcMC_met_blnotag_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",        ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baselineNoTag_elmu",    {dcData_2015CD_ht_blnotag_elmu, dcMC_ht_blnotag_elmu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baselineNoTag_elmu",    {dcData_2015CD_nt_blnotag_elmu, dcMC_nt_blnotag_elmu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",       ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baselineNoTag_elmu",   {dcData_2015CD_mt2_blnotag_elmu, dcMC_mt2_blnotag_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",        ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baselineNoTag_elmu",    {dcData_2015CD_nb_blnotag_elmu, dcMC_nb_blnotag_elmu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baselineNoTag_elmu",    {dcData_2015CD_nj_blnotag_elmu, dcMC_nj_blnotag_elmu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baselineNoTag_elmu",   {dcData_2015CD_mht_blnotag_elmu, dcMC_mht_blnotag_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",        ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baselineNoTag_elmu",   {dcData_2015CD_jpt_blnotag_elmu, dcMC_jpt_blnotag_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baselineNoTag_elmu",  {dcData_2015CD_j1pt_blnotag_elmu, dcMC_j1pt_blnotag_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baselineNoTag_elmu",  {dcData_2015CD_j2pt_blnotag_elmu, dcMC_j2pt_blnotag_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baselineNoTag_elmu",  {dcData_2015CD_j3pt_blnotag_elmu, dcMC_j3pt_blnotag_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",    ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baselineNoTag_elmu", {dcData_2015CD_mu1pt_blnotag_elmu, dcMC_mu1pt_blnotag_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_baselineNoTag_elmu", {dcData_2015CD_el1pt_blnotag_elmu, dcMC_el1pt_blnotag_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",     ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baselineNoTag_elmu",   {dcData_2015CD_mll_blnotag_elmu, dcMC_mll_blnotag_elmu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",        ""));
    // --> baseline
    vh.push_back(PHS("DataMC_2015CD_met_baseline_elmu",   {dcData_2015CD_met_bl_elmu, dcMC_met_bl_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",            ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baseline_elmu",    {dcData_2015CD_ht_bl_elmu, dcMC_ht_bl_elmu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baseline_elmu",    {dcData_2015CD_nt_bl_elmu, dcMC_nt_bl_elmu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",           ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baseline_elmu",   {dcData_2015CD_mt2_bl_elmu, dcMC_mt2_bl_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",            ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baseline_elmu",    {dcData_2015CD_nb_bl_elmu, dcMC_nb_bl_elmu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",             ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baseline_elmu",    {dcData_2015CD_nj_bl_elmu, dcMC_nj_bl_elmu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",             ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baseline_elmu",   {dcData_2015CD_mht_bl_elmu, dcMC_mht_bl_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",            ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baseline_elmu",   {dcData_2015CD_jpt_bl_elmu, dcMC_jpt_bl_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baseline_elmu",  {dcData_2015CD_j1pt_bl_elmu, dcMC_j1pt_bl_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baseline_elmu",  {dcData_2015CD_j2pt_bl_elmu, dcMC_j2pt_bl_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baseline_elmu",  {dcData_2015CD_j3pt_bl_elmu, dcMC_j3pt_bl_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",        ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baseline_elmu", {dcData_2015CD_mu1pt_bl_elmu, dcMC_mu1pt_bl_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_baseline_elmu", {dcData_2015CD_el1pt_bl_elmu, dcMC_el1pt_bl_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mll_baseline_elmu",   {dcData_2015CD_mll_bl_elmu, dcMC_mll_bl_elmu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",            ""));
    // --> 0b muZinv
    vh.push_back(PHS("DataMC_2015CD_met_0belmuZinv",   {dcData_2015CD_met_0belmuZinv, dcMC_met_0belmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0belmuZinv",    {dcData_2015CD_ht_0belmuZinv, dcMC_ht_0belmuZinv},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0belmuZinv",    {dcData_2015CD_nt_0belmuZinv, dcMC_nt_0belmuZinv},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0belmuZinv",  {dcData_2015CD_nt1b_0belmuZinv, dcMC_nt1b_0belmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0belmuZinv",  {dcData_2015CD_nt2b_0belmuZinv, dcMC_nt2b_0belmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0belmuZinv",  {dcData_2015CD_nt3b_0belmuZinv, dcMC_nt3b_0belmuZinv},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0belmuZinv",   {dcData_2015CD_mt2_0belmuZinv, dcMC_mt2_0belmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0belmuZinv", {dcData_2015CD_mt21b_0belmuZinv, dcMC_mt21b_0belmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0belmuZinv", {dcData_2015CD_mt22b_0belmuZinv, dcMC_mt22b_0belmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0belmuZinv", {dcData_2015CD_mt23b_0belmuZinv, dcMC_mt23b_0belmuZinv}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    //vh.push_back(PHS("DataMC_2015CD_nb_0belmuZinv",    {dcData_2015CD_nb_0belmuZinv, dcMC_nb_0belmuZinv},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0belmuZinv",    {dcData_2015CD_nj_0belmuZinv, dcMC_nj_0belmuZinv},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0belmuZinv",   {dcData_2015CD_mht_0belmuZinv, dcMC_mht_0belmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0belmuZinv",   {dcData_2015CD_jpt_0belmuZinv, dcMC_jpt_0belmuZinv},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0belmuZinv",  {dcData_2015CD_j1pt_0belmuZinv, dcMC_j1pt_0belmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0belmuZinv",  {dcData_2015CD_j2pt_0belmuZinv, dcMC_j2pt_0belmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0belmuZinv",  {dcData_2015CD_j3pt_0belmuZinv, dcMC_j3pt_0belmuZinv},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0belmuZinv", {dcData_2015CD_mu1pt_0belmuZinv, dcMC_mu1pt_0belmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0belmuZinv", {dcData_2015CD_el1pt_0belmuZinv, dcMC_el1pt_0belmuZinv}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0belmuZinv",   {dcData_2015CD_mll_0belmuZinv, dcMC_mll_0belmuZinv},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));
    // --> 0b no mt2
    vh.push_back(PHS("DataMC_2015CD_met_0bnomt2_elmu",   {dcData_2015CD_met_0bnomt2_elmu, dcMC_met_0bnomt2_elmu},     {1, 2}, "", 50, 0, 1500,  true, false,  "met",             ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0bnomt2_elmu",    {dcData_2015CD_ht_0bnomt2_elmu, dcMC_ht_0bnomt2_elmu},       {1, 2}, "", 50, 0, 1500,  true, false,  "ht",              ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0bnomt2_elmu",    {dcData_2015CD_nt_0bnomt2_elmu, dcMC_nt_0bnomt2_elmu},       {1, 2}, "", 5, 0, 5,      true, false,  "ntop",            ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0bnomt2_elmu",  {dcData_2015CD_nt1b_0bnomt2_elmu, dcMC_nt1b_0bnomt2_elmu},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0bnomt2_elmu",  {dcData_2015CD_nt2b_0bnomt2_elmu, dcMC_nt2b_0bnomt2_elmu},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0bnomt2_elmu",  {dcData_2015CD_nt3b_0bnomt2_elmu, dcMC_nt3b_0bnomt2_elmu},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0bnomt2_elmu",   {dcData_2015CD_mt2_0bnomt2_elmu, dcMC_mt2_0bnomt2_elmu},     {1, 2}, "", 50, 0, 1500,  true, false,  "mt2",             ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0bnomt2_elmu", {dcData_2015CD_mt21b_0bnomt2_elmu, dcMC_mt21b_0bnomt2_elmu}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0bnomt2_elmu", {dcData_2015CD_mt22b_0bnomt2_elmu, dcMC_mt22b_0bnomt2_elmu}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0bnomt2_elmu", {dcData_2015CD_mt23b_0bnomt2_elmu, dcMC_mt23b_0bnomt2_elmu}, {1, 2}, "", 50, 0, 1500,  true, false,  "mt2 (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0bnomt2_elmu",    {dcData_2015CD_nb_0bnomt2_elmu, dcMC_nb_0bnomt2_elmu},       {1, 2}, "", 10, 0, 10,    true, false,  "nb",              ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0bnomt2_elmu",    {dcData_2015CD_nj_0bnomt2_elmu, dcMC_nj_0bnomt2_elmu},       {1, 2}, "", 20, 0, 20,    true, false,  "nj",              ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0bnomt2_elmu",   {dcData_2015CD_mht_0bnomt2_elmu, dcMC_mht_0bnomt2_elmu},     {1, 2}, "", 50, 0, 1500,  true, false,  "mht",             ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0bnomt2_elmu",   {dcData_2015CD_jpt_0bnomt2_elmu, dcMC_jpt_0bnomt2_elmu},     {1, 2}, "", 50, 0, 1500,  true, false,  "jet pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0bnomt2_elmu",  {dcData_2015CD_j1pt_0bnomt2_elmu, dcMC_j1pt_0bnomt2_elmu},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0bnomt2_elmu",  {dcData_2015CD_j2pt_0bnomt2_elmu, dcMC_j2pt_0bnomt2_elmu},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0bnomt2_elmu",  {dcData_2015CD_j3pt_0bnomt2_elmu, dcMC_j3pt_0bnomt2_elmu},   {1, 2}, "", 50, 0, 1500,  true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0bnomt2_elmu", {dcData_2015CD_mu1pt_0bnomt2_elmu, dcMC_mu1pt_0bnomt2_elmu}, {1, 2}, "", 50, 0, 1000,  true, false,  "mu1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0bnomt2_elmu", {dcData_2015CD_el1pt_0bnomt2_elmu, dcMC_el1pt_0bnomt2_elmu}, {1, 2}, "", 50, 0, 1000,  true, false,  "el1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0bnomt2_elmu",   {dcData_2015CD_mll_0bnomt2_elmu, dcMC_mll_0bnomt2_elmu},     {1, 2}, "", 50, 0, 1000,  true, false,  "mll",             ""));
    // --> 0b
    vh.push_back(PHS("DataMC_2015CD_met_0b_elmu",   {dcData_2015CD_met_0b_elmu, dcMC_met_0b_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "met",              ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0b_elmu",    {dcData_2015CD_ht_0b_elmu, dcMC_ht_0b_elmu},       {1, 2}, "", 50, 0, 1500,   true, false,  "ht",               ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0b_elmu",    {dcData_2015CD_nt_0b_elmu, dcMC_nt_0b_elmu},       {1, 2}, "", 5, 0, 5,       true, false,  "ntop",             ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0b_elmu",  {dcData_2015CD_nt1b_0b_elmu, dcMC_nt1b_0b_elmu},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0b_elmu",  {dcData_2015CD_nt2b_0b_elmu, dcMC_nt2b_0b_elmu},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0b_elmu",  {dcData_2015CD_nt3b_0b_elmu, dcMC_nt3b_0b_elmu},   {1, 2}, "", 5, 0, 5,       true, false,  "ntop (3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0b_elmu",   {dcData_2015CD_mt2_0b_elmu, dcMC_mt2_0b_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mt2",              ""));
    vh.push_back(PHS("DataMC_2015CD_mt21b_0b_elmu", {dcData_2015CD_mt21b_0b_elmu, dcMC_mt21b_0b_elmu}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (1b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt22b_0b_elmu", {dcData_2015CD_mt22b_0b_elmu, dcMC_mt22b_0b_elmu}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (2b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_mt23b_0b_elmu", {dcData_2015CD_mt23b_0b_elmu, dcMC_mt23b_0b_elmu}, {1, 2}, "", 50, 0, 1500,   true, false,  "mt2 (3b fake)",    ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0b_elmu",    {dcData_2015CD_nb_0b_elmu, dcMC_nb_0b_elmu},       {1, 2}, "", 10, 0, 10,     true, false,  "nb",               ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0b_elmu",    {dcData_2015CD_nj_0b_elmu, dcMC_nj_0b_elmu},       {1, 2}, "", 20, 0, 20,     true, false,  "nj",               ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0b_elmu",   {dcData_2015CD_mht_0b_elmu, dcMC_mht_0b_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "mht",              ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0b_elmu",   {dcData_2015CD_jpt_0b_elmu, dcMC_jpt_0b_elmu},     {1, 2}, "", 50, 0, 1500,   true, false,  "jet pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0b_elmu",  {dcData_2015CD_j1pt_0b_elmu, dcMC_j1pt_0b_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet1 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0b_elmu",  {dcData_2015CD_j2pt_0b_elmu, dcMC_j2pt_0b_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet2 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0b_elmu",  {dcData_2015CD_j3pt_0b_elmu, dcMC_j3pt_0b_elmu},   {1, 2}, "", 50, 0, 1500,   true, false,  "jet3 pt",          ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0b_elmu", {dcData_2015CD_mu1pt_0b_elmu, dcMC_mu1pt_0b_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "mu1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_el1pt_0b_elmu", {dcData_2015CD_el1pt_0b_elmu, dcMC_el1pt_0b_elmu}, {1, 2}, "", 50, 0, 1000,   true, false,  "el1 pt",           ""));
    vh.push_back(PHS("DataMC_2015CD_mll_0b_elmu",   {dcData_2015CD_mll_0b_elmu, dcMC_mll_0b_elmu},     {1, 2}, "", 50, 0, 1000,   true, false,  "mll",              ""));


    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.read();
    if(doSave)  plotter.saveHists();
    if(doPlots) plotter.plot();
}
