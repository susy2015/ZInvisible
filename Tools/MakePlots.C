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
        if(sc[dataSets] != sc.null())
        {
            fileMap[dataSets] = {sc[dataSets]};
            int i = 0;
            for(const auto& fs : sc[dataSets])
            {
                fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
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
	// If we didn't find such a collection, we pass in the dataset itself
        //if (fileMap.size() == 0)
	//{
	//    fileMap[dataSets] = {ss[dataSets]};
	//	}

    }



	/*
    else if(dataSets.compare("Data_SingleMuon_2015C") == 0)
    {
        //add this to the 25ns filemap
        fileMap[dataSets] = {ss[dataSets]};
	fileMap["Data_SingleMuon25ns"] = {ss[dataSets]};
    }
    else if(dataSets.compare("Data_SingleMuon_2015D") == 0)
    {
        //add this to the 25ns filemap
        fileMap[dataSets] = {ss[dataSets]};
	fileMap["Data_SingleMuon25ns"] = {ss[dataSets]};
    }
	*/



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

    // ------------------------
    // - Data/MC plots
    // ------------------------

    // DatasetSummary for each cut level and process
    // --> No selection apart from noise filters
    Plotter::DatasetSummary dsData_2015C_nosel("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsDY_nosel("DY", fileMap["DYJetsToLL"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dstt2l_nosel("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dstW_nosel("single top", fileMap["tW"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsttZ_nosel("t#bar{t}Z", fileMap["TTZ"], "passNoiseEventFilterZinv", "");
    Plotter::DatasetSummary dsVV_nosel("Diboson", fileMap["Diboson"], "passNoiseEventFilterZinv", "");
    // --> baseline without btag, top tag and mt2 cuts
    Plotter::DatasetSummary dsData_2015C_blnotag("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_blnotag("DY", fileMap["DYJetsToLL"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_blnotag("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstW_blnotag("single top", fileMap["tW"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_blnotag("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsVV_blnotag("Diboson", fileMap["Diboson"], "passBaselineNoTagMT2Zinv;passMuZinvSel", "");
    // --> full baseline
    Plotter::DatasetSummary dsData_2015C_bl("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsDY_bl("DY", fileMap["DYJetsToLL"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstt2l_bl("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dstW_bl("single top", fileMap["tW"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsttZ_bl("t#bar{t}Z", fileMap["TTZ"], "passBaselineZinv;passMuZinvSel", "");
    Plotter::DatasetSummary dsVV_bl("Diboson", fileMap["Diboson"], "passBaselineZinv;passMuZinvSel", "");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0bnomt2("DY", fileMap["DYJetsToLL"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0bnomt2("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0bnomt2("single top", fileMap["tW"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0bnomt2("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0bnomt2("Diboson", fileMap["Diboson"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsDY_0b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstt2l_0b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dstW_0b("single top", fileMap["tW"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsttZ_0b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");
    Plotter::DatasetSummary dsVV_0b("Diboson", fileMap["Diboson"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "");

    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w1b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w1b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w1b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w1b("single top", fileMap["tW"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w1b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w1b("Diboson", fileMap["Diboson"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w1b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsDY_0b_w1b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w1b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dstW_0b_w1b("single top", fileMap["tW"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w1b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    Plotter::DatasetSummary dsVV_0b_w1b("Diboson", fileMap["Diboson"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight1fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w2b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w2b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w2b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w2b("single top", fileMap["tW"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w2b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w2b("Diboson", fileMap["Diboson"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w2b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsDY_0b_w2b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w2b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dstW_0b_w2b("single top", fileMap["tW"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w2b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    Plotter::DatasetSummary dsVV_0b_w2b("Diboson", fileMap["Diboson"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight2fakebComb");
    // --> 0b control region without mt2 cut
    Plotter::DatasetSummary dsData_2015C_0bnomt2_w3b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0bnomt2_w3b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0bnomt2_w3b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0bnomt2_w3b("single top", fileMap["tW"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0bnomt2_w3b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0bnomt2_w3b("Diboson", fileMap["Diboson"], "passBaselineNoTagMT2Zinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    // --> 0b control region 
    Plotter::DatasetSummary dsData_2015C_0b_w3b("Data Run 2015C", fileMap["Data_SingleMuon25ns"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsDY_0b_w3b("DY", fileMap["DYJetsToLL"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstt2l_0b_w3b("t#bar{t} dilepton", fileMap["TTbarDiLep"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dstW_0b_w3b("single top", fileMap["tW"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsttZ_0b_w3b("t#bar{t}Z", fileMap["TTZ"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");
    Plotter::DatasetSummary dsVV_0b_w3b("Diboson", fileMap["Diboson"], "passBaselineNoTagZinv;passMuZinvSel;cntCSVSZinv=0", "weight3fakebComb");


    // Define the collections, i.e. variables per selection
    // --> met
    Plotter::DataCollection dcData_2015CD_met_nosel("data",  "cleanMetPt", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_met_nosel("stack",  "cleanMetPt", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_met_blnotag("data",  "cleanMetPt", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_met_blnotag("stack",  "cleanMetPt", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_met_bl("data",  "cleanMetPt", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_met_bl("stack",  "cleanMetPt", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_met_0bnomt2("data",  "cleanMetPt", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_met_0bnomt2("stack",  "cleanMetPt", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_met_0b("data",  "cleanMetPt", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_met_0b("stack",  "cleanMetPt", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> ntops
    Plotter::DataCollection dcData_2015CD_nt_nosel("data",  "nTopCandSortedCntZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_nt_nosel("stack",  "nTopCandSortedCntZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nt_blnotag("data",  "nTopCandSortedCntZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_nt_blnotag("stack",  "nTopCandSortedCntZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nt_bl("data",  "nTopCandSortedCntZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_nt_bl("stack",  "nTopCandSortedCntZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nt_0bnomt2("data",  "nTopCandSortedCntZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_nt_0bnomt2("stack",  "nTopCandSortedCntZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nt_0b("data",  "nTopCandSortedCntZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_nt_0b("stack",  "nTopCandSortedCntZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> ntops with bfaking
    Plotter::DataCollection dcData_2015CD_nt1b_0bnomt2("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0bnomt2_w1b});
    Plotter::DataCollection dcDY_nt1b_0bnomt2("stack",  "nTopCandSortedCntZinv1b", {dsDY_0bnomt2_w1b, dstt2l_0bnomt2_w1b, dstW_0bnomt2_w1b, dsttZ_0bnomt2_w1b, dsVV_0bnomt2_w1b});
    Plotter::DataCollection dcData_2015CD_nt1b_0b("data", "nTopCandSortedCntZinv1b", {dsData_2015C_0b_w1b});
    Plotter::DataCollection dcDY_nt1b_0b("stack",  "nTopCandSortedCntZinv1b", {dsDY_0b_w1b, dstt2l_0b_w1b, dstW_0b_w1b, dsttZ_0b_w1b, dsVV_0b_w1b});
    Plotter::DataCollection dcData_2015CD_nt2b_0bnomt2("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0bnomt2_w2b});
    Plotter::DataCollection dcDY_nt2b_0bnomt2("stack",  "nTopCandSortedCntZinv2b", {dsDY_0bnomt2_w2b, dstt2l_0bnomt2_w2b, dstW_0bnomt2_w2b, dsttZ_0bnomt2_w2b, dsVV_0bnomt2_w2b});
    Plotter::DataCollection dcData_2015CD_nt2b_0b("data",  "nTopCandSortedCntZinv2b", {dsData_2015C_0b_w2b});
    Plotter::DataCollection dcDY_nt2b_0b("stack",  "nTopCandSortedCntZinv2b", {dsDY_0b_w2b, dstt2l_0b_w2b, dstW_0b_w2b, dsttZ_0b_w2b, dsVV_0b_w2b});
    Plotter::DataCollection dcData_2015CD_nt3b_0bnomt2("data",  "nTopCandSortedCntZinv3b", {dsData_2015C_0bnomt2_w3b});
    Plotter::DataCollection dcDY_nt3b_0bnomt2("stack",  "nTopCandSortedCntZinv3b", {dsDY_0bnomt2_w3b, dstt2l_0bnomt2_w3b, dstW_0bnomt2_w3b, dsttZ_0bnomt2_w3b, dsVV_0bnomt2_w3b});
    Plotter::DataCollection dcData_2015CD_nt3b_0b("data", "nTopCandSortedCntZinv3b", {dsData_2015C_0b_w3b});
    Plotter::DataCollection dcDY_nt3b_0b("stack",  "nTopCandSortedCntZinv3b", {dsDY_0b_w3b, dstt2l_0b_w3b, dstW_0b_w3b, dsttZ_0b_w3b, dsVV_0b_w3b});
    // --> nbs
    Plotter::DataCollection dcData_2015CD_nb_nosel("data",  "cntCSVSZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_nb_nosel("stack",  "cntCSVSZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nb_blnotag("data",  "cntCSVSZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_nb_blnotag("stack",  "cntCSVSZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nb_bl("data",  "cntCSVSZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_nb_bl("stack",  "cntCSVSZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nb_0bnomt2("data",  "cntCSVSZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_nb_0bnomt2("stack",  "cntCSVSZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nb_0b("data",  "cntCSVSZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_nb_0b("stack",  "cntCSVSZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> njets
    Plotter::DataCollection dcData_2015CD_nj_nosel("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_nj_nosel("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_nj_blnotag("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_nj_blnotag("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_nj_bl("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_nj_bl("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_nj_0bnomt2("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_nj_0bnomt2("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_nj_0b("data",  "cntNJetsPt30Eta24Zinv", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_nj_0b("stack",  "cntNJetsPt30Eta24Zinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> HT
    Plotter::DataCollection dcData_2015CD_ht_nosel("data",  "HTZinv", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_ht_nosel("stack",  "HTZinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_ht_blnotag("data",  "HTZinv", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_ht_blnotag("stack",  "HTZinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_ht_bl("data",  "HTZinv", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_ht_bl("stack",  "HTZinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_ht_0bnomt2("data",  "HTZinv", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_ht_0bnomt2("stack",  "HTZinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_ht_0b("data",  "HTZinv", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_ht_0b("stack",  "HTZinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> MHT
    Plotter::DataCollection dcData_2015CD_mht_nosel("data",  "cleanMHt", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_mht_nosel("stack",  "cleanMHt", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_mht_blnotag("data",  "cleanMHt", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_mht_blnotag("stack",  "cleanMHt", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_mht_bl("data",  "cleanMHt", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_mht_bl("stack",  "cleanMHt", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_mht_0bnomt2("data",  "cleanMHt", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_mht_0bnomt2("stack",  "cleanMHt", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_mht_0b("data",  "cleanMHt", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_mht_0b("stack",  "cleanMHt", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet pt
    Plotter::DataCollection dcData_2015CD_jpt_nosel("data",  "cleanJetVec(pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_jpt_nosel("stack",  "cleanJetVec(pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_jpt_blnotag("data",  "cleanJetVec(pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_jpt_blnotag("stack",  "cleanJetVec(pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_jpt_bl("data",  "cleanJetVec(pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_jpt_bl("stack",  "cleanJetVec(pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_jpt_0bnomt2("data",  "cleanJetVec(pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_jpt_0bnomt2("stack",  "cleanJetVec(pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_jpt_0b("data",  "cleanJetVec(pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_jpt_0b("stack",  "cleanJetVec(pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet1 pt
    Plotter::DataCollection dcData_2015CD_j1pt_nosel("data",  "cleanJetVec[0](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_j1pt_nosel("stack",  "cleanJetVec[0](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j1pt_blnotag("data",  "cleanJetVec[0](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_j1pt_blnotag("stack",  "cleanJetVec[0](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j1pt_bl("data",  "cleanJetVec[0](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_j1pt_bl("stack",  "cleanJetVec[0](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j1pt_0bnomt2("data",  "cleanJetVec[0](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_j1pt_0bnomt2("stack",  "cleanJetVec[0](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j1pt_0b("data",  "cleanJetVec[0](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_j1pt_0b("stack",  "cleanJetVec[0](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet2 pt
    Plotter::DataCollection dcData_2015CD_j2pt_nosel("data",  "cleanJetVec[1](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_j2pt_nosel("stack",  "cleanJetVec[1](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j2pt_blnotag("data",  "cleanJetVec[1](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_j2pt_blnotag("stack",  "cleanJetVec[1](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j2pt_bl("data",  "cleanJetVec[1](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_j2pt_bl("stack",  "cleanJetVec[1](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j2pt_0bnomt2("data",  "cleanJetVec[1](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_j2pt_0bnomt2("stack",  "cleanJetVec[1](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j2pt_0b("data",  "cleanJetVec[1](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_j2pt_0b("stack",  "cleanJetVec[1](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});
    // --> jet3 pt
    Plotter::DataCollection dcData_2015CD_j3pt_nosel("data",  "cleanJetVec[2](pt)", {dsData_2015C_nosel});
    Plotter::DataCollection dcDY_j3pt_nosel("stack",  "cleanJetVec[2](pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel});
    Plotter::DataCollection dcData_2015CD_j3pt_blnotag("data",  "cleanJetVec[2](pt)", {dsData_2015C_blnotag});
    Plotter::DataCollection dcDY_j3pt_blnotag("stack",  "cleanJetVec[2](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag});
    Plotter::DataCollection dcData_2015CD_j3pt_bl("data",  "cleanJetVec[2](pt)", {dsData_2015C_bl});
    Plotter::DataCollection dcDY_j3pt_bl("stack",  "cleanJetVec[2](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl});
    Plotter::DataCollection dcData_2015CD_j3pt_0bnomt2("data",  "cleanJetVec[2](pt)", {dsData_2015C_0bnomt2});
    Plotter::DataCollection dcDY_j3pt_0bnomt2("stack",  "cleanJetVec[2](pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2});
    Plotter::DataCollection dcData_2015CD_j3pt_0b("data",  "cleanJetVec[2](pt)", {dsData_2015C_0b});
    Plotter::DataCollection dcDY_j3pt_0b("stack",  "cleanJetVec[2](pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b});


    // Push back the histograms into the histo vector
    // --> nosel
    vh.push_back(PHS("DataMC_2015CD_met_nosel", {dcData_2015CD_met_nosel, dcDY_met_nosel}, {1, 2}, "", 150, 0, 1500, true, false,  "met",    ""));
    vh.push_back(PHS("DataMC_2015CD_ht_nosel",  {dcData_2015CD_ht_nosel, dcDY_ht_nosel},   {1, 2}, "", 150, 0, 1500, true, false,  "ht",     ""));
    vh.push_back(PHS("DataMC_2015CD_nt_nosel",  {dcData_2015CD_nt_nosel, dcDY_nt_nosel},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt1b_nosel",  {dcData_2015CD_nt1b_nosel, dcDY_nt1b_nosel},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt2b_nosel",  {dcData_2015CD_nt2b_nosel, dcDY_nt2b_nosel},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt3b_nosel",  {dcData_2015CD_nt3b_nosel, dcDY_nt3b_nosel},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_nosel",  {dcData_2015CD_nb_nosel, dcDY_nb_nosel},   {1, 2}, "", 10, 0, 10,    true, false,  "nb",     ""));
    vh.push_back(PHS("DataMC_2015CD_nj_nosel",  {dcData_2015CD_nj_nosel, dcDY_nj_nosel},   {1, 2}, "", 20, 0, 20,    true, false,  "nj",     ""));
    vh.push_back(PHS("DataMC_2015CD_mht_nosel", {dcData_2015CD_mht_nosel, dcDY_mht_nosel}, {1, 2}, "", 150, 0, 1500, true, false,  "mht",    ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_nosel", {dcData_2015CD_jpt_nosel, dcDY_ht_nosel},  {1, 2}, "", 150, 0, 1500, true, false,  "jet pt", ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_nosel", {dcData_2015CD_j1pt_nosel, dcDY_j1pt_nosel},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_nosel", {dcData_2015CD_j2pt_nosel, dcDY_j2pt_nosel}, {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_nosel", {dcData_2015CD_j3pt_nosel, dcDY_j3pt_nosel},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_nosel",  
		     {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_nosel}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_nosel",  
		     {PDC("data", "cutMuVec(pt)", {dsData_2015C_nosel}), PDC("stack", "cutMuVec(pt)", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_nosel",  
		     {PDC("data", "cutMuPt1", {dsData_2015C_nosel}), PDC("stack", "cutMuPt1", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_nosel",  
		     {PDC("data", "cutMuPt2", {dsData_2015C_nosel}), PDC("stack", "cutMuPt2", {dsDY_nosel, dstt2l_nosel, dstW_nosel, dsttZ_nosel, dsVV_nosel})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    // --> baselineNoTag
    vh.push_back(PHS("DataMC_2015CD_met_baselineNoTag", {dcData_2015CD_met_blnotag, dcDY_met_blnotag},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baselineNoTag",  {dcData_2015CD_ht_blnotag, dcDY_ht_blnotag},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baselineNoTag",  {dcData_2015CD_nt_blnotag, dcDY_nt_blnotag},  {1, 2}, "", 5, 0, 5,   true, false,  "ntop",         ""));
    //vh.push_back(PHS("DataMC_2015CD_nt1b_baselineNoTag",  {dcData_2015CD_nt1b_baselineNoTag, dcDY_nt1b_baselineNoTag},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt2b_baselineNoTag",  {dcData_2015CD_nt2b_baselineNoTag, dcDY_nt2b_baselineNoTag},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt3b_baselineNoTag",  {dcData_2015CD_nt3b_baselineNoTag, dcDY_nt3b_baselineNoTag},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baselineNoTag",  {dcData_2015CD_nb_blnotag, dcDY_nb_blnotag},  {1, 2}, "", 10, 0, 10,   true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baselineNoTag",  {dcData_2015CD_nj_blnotag, dcDY_nj_blnotag},  {1, 2}, "", 20, 0, 20,   true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baselineNoTag",  
		     {PDC("data", "cleanMHt", {dsData_2015C_blnotag}), PDC("stack", "cleanMHt", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baselineNoTag",  
		     {PDC("data", "cleanJetVec(pt)", {dsData_2015C_blnotag}), PDC("stack", "cleanJetVec(pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baselineNoTag",  
		     {PDC("data", "cleanJetVec[0](pt)", {dsData_2015C_blnotag}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baselineNoTag",  
		     {PDC("data", "cleanJetVec[1](pt)", {dsData_2015C_blnotag}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baselineNoTag",  
		     {PDC("data", "cleanJetVec[2](pt)", {dsData_2015C_blnotag}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baselineNoTag",  
		     {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_blnotag}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baselineNoTag",  
		     {PDC("data", "cutMuVec(pt)", {dsData_2015C_blnotag}), PDC("stack", "cutMuVec(pt)", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baselineNoTag",  
		     {PDC("data", "cutMuPt1", {dsData_2015C_blnotag}), PDC("stack", "cutMuPt1", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baselineNoTag",  
		     {PDC("data", "cutMuPt2", {dsData_2015C_blnotag}), PDC("stack", "cutMuPt2", {dsDY_blnotag, dstt2l_blnotag, dstW_blnotag, dsttZ_blnotag, dsVV_blnotag})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    // --> baseline
    vh.push_back(PHS("DataMC_2015CD_met_baseline", {dcData_2015CD_met_bl, dcDY_met_bl},  {1, 2}, "", 150, 0, 1500,   true, false,  "met",         ""));
    vh.push_back(PHS("DataMC_2015CD_ht_baseline",  {dcData_2015CD_ht_bl, dcDY_ht_bl},  {1, 2}, "", 150, 0, 1500,   true, false,  "ht",         ""));
    vh.push_back(PHS("DataMC_2015CD_nt_baseline",  {dcData_2015CD_nt_bl, dcDY_nt_bl},  {1, 2}, "", 5, 0, 5,   true, false,  "ntop",         ""));
    //vh.push_back(PHS("DataMC_2015CD_nt1b_baseline",  {dcData_2015CD_nt1b_baseline, dcDY_nt1b_baseline},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt2b_baseline",  {dcData_2015CD_nt2b_baseline, dcDY_nt2b_baseline},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    //vh.push_back(PHS("DataMC_2015CD_nt3b_baseline",  {dcData_2015CD_nt3b_baseline, dcDY_nt3b_baseline},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_baseline",  {dcData_2015CD_nb_bl, dcDY_nb_bl},  {1, 2}, "", 10, 0, 10,   true, false,  "nb",         ""));
    vh.push_back(PHS("DataMC_2015CD_nj_baseline",  {dcData_2015CD_nj_bl, dcDY_nj_bl},  {1, 2}, "", 20, 0, 20,   true, false,  "nj",         ""));
    vh.push_back(PHS("DataMC_2015CD_mht_baseline",  
		     {PDC("data", "cleanMHt", {dsData_2015C_bl}), PDC("stack", "cleanMHt", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mht",         ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_baseline",  
		     {PDC("data", "cleanJetVec(pt)", {dsData_2015C_bl}), PDC("stack", "cleanJetVec(pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_baseline",  
		     {PDC("data", "cleanJetVec[0](pt)", {dsData_2015C_bl}), PDC("stack", "cleanJetVec[0](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_baseline",  
		     {PDC("data", "cleanJetVec[1](pt)", {dsData_2015C_bl}), PDC("stack", "cleanJetVec[1](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_baseline",  
		     {PDC("data", "cleanJetVec[2](pt)", {dsData_2015C_bl}), PDC("stack", "cleanJetVec[2](pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_baseline",  
		     {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_bl}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_baseline",  
		     {PDC("data", "cutMuVec(pt)", {dsData_2015C_bl}), PDC("stack", "cutMuVec(pt)", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_baseline",  
		     {PDC("data", "cutMuPt1", {dsData_2015C_bl}), PDC("stack", "cutMuPt1", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_baseline",  
		     {PDC("data", "cutMuPt2", {dsData_2015C_bl}), PDC("stack", "cutMuPt2", {dsDY_bl, dstt2l_bl, dstW_bl, dsttZ_bl, dsVV_bl})},  
    		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    // --> 0b no mt2
    vh.push_back(PHS("DataMC_2015CD_met_0bnomt2", {dcData_2015CD_met_0bnomt2, dcDY_met_0bnomt2}, {1, 2}, "", 150, 0, 1500, true, false,  "met",    ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0bnomt2",  {dcData_2015CD_ht_0bnomt2, dcDY_ht_0bnomt2},   {1, 2}, "", 150, 0, 1500, true, false,  "ht",     ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0bnomt2",  {dcData_2015CD_nt_0bnomt2, dcDY_nt_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0bnomt2",  {dcData_2015CD_nt1b_0bnomt2, dcDY_nt1b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0bnomt2",  {dcData_2015CD_nt2b_0bnomt2, dcDY_nt2b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0bnomt2",  {dcData_2015CD_nt3b_0bnomt2, dcDY_nt3b_0bnomt2},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0bnomt2",  {dcData_2015CD_nb_0bnomt2, dcDY_nb_0bnomt2},   {1, 2}, "", 10, 0, 10,    true, false,  "nb",     ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0bnomt2",  {dcData_2015CD_nj_0bnomt2, dcDY_nj_0bnomt2},   {1, 2}, "", 20, 0, 20,    true, false,  "nj",     ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0bnomt2", {dcData_2015CD_mht_0bnomt2, dcDY_mht_0bnomt2}, {1, 2}, "", 150, 0, 1500, true, false,  "mht",    ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0bnomt2", {dcData_2015CD_jpt_0bnomt2, dcDY_ht_0bnomt2},  {1, 2}, "", 150, 0, 1500, true, false,  "jet pt", ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0bnomt2", {dcData_2015CD_j1pt_0bnomt2, dcDY_j1pt_0bnomt2},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0bnomt2", {dcData_2015CD_j2pt_0bnomt2, dcDY_j2pt_0bnomt2}, {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0bnomt2", {dcData_2015CD_j3pt_0bnomt2, dcDY_j3pt_0bnomt2},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0bnomt2",  
		     {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_0bnomt2}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_0bnomt2",  
		     {PDC("data", "cutMuVec(pt)", {dsData_2015C_0bnomt2}), PDC("stack", "cutMuVec(pt)", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0bnomt2",  
		     {PDC("data", "cutMuPt1", {dsData_2015C_0bnomt2}), PDC("stack", "cutMuPt1", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_0bnomt2",  
		     {PDC("data", "cutMuPt2", {dsData_2015C_0bnomt2}), PDC("stack", "cutMuPt2", {dsDY_0bnomt2, dstt2l_0bnomt2, dstW_0bnomt2, dsttZ_0bnomt2, dsVV_0bnomt2})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));
    // --> 0b 
    vh.push_back(PHS("DataMC_2015CD_met_0b", {dcData_2015CD_met_0b, dcDY_met_0b}, {1, 2}, "", 150, 0, 1500, true, false,  "met",    ""));
    vh.push_back(PHS("DataMC_2015CD_ht_0b",  {dcData_2015CD_ht_0b, dcDY_ht_0b},   {1, 2}, "", 150, 0, 1500, true, false,  "ht",     ""));
    vh.push_back(PHS("DataMC_2015CD_nt_0b",  {dcData_2015CD_nt_0b, dcDY_nt_0b},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt1b_0b",  {dcData_2015CD_nt1b_0b, dcDY_nt1b_0b},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(1b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt2b_0b",  {dcData_2015CD_nt2b_0b, dcDY_nt2b_0b},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(2b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nt3b_0b",  {dcData_2015CD_nt3b_0b, dcDY_nt3b_0b},   {1, 2}, "", 5, 0, 5,      true, false,  "ntop(3b fake)",   ""));
    vh.push_back(PHS("DataMC_2015CD_nb_0b",  {dcData_2015CD_nb_0b, dcDY_nb_0b},   {1, 2}, "", 10, 0, 10,    true, false,  "nb",     ""));
    vh.push_back(PHS("DataMC_2015CD_nj_0b",  {dcData_2015CD_nj_0b, dcDY_nj_0b},   {1, 2}, "", 20, 0, 20,    true, false,  "nj",     ""));
    vh.push_back(PHS("DataMC_2015CD_mht_0b", {dcData_2015CD_mht_0b, dcDY_mht_0b}, {1, 2}, "", 150, 0, 1500, true, false,  "mht",    ""));
    vh.push_back(PHS("DataMC_2015CD_jpt_0b", {dcData_2015CD_jpt_0b, dcDY_ht_0b},  {1, 2}, "", 150, 0, 1500, true, false,  "jet pt", ""));
    vh.push_back(PHS("DataMC_2015CD_j1pt_0b", {dcData_2015CD_j1pt_0b, dcDY_j1pt_0b},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j2pt_0b", {dcData_2015CD_j2pt_0b, dcDY_j2pt_0b}, {1, 2}, "", 150, 0, 1500,   true, false,  "jet2 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_j3pt_0b", {dcData_2015CD_j3pt_0b, dcDY_j3pt_0b},  {1, 2}, "", 150, 0, 1500,   true, false,  "jet3 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mt2_0b",  
		     {PDC("data", "best_had_brJet_MT2Zinv", {dsData_2015C_0b}), PDC("stack", "best_had_brJet_MT2Zinv", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b})},  
		     {1, 2}, "", 150, 0, 1500,   true, false,  "mt2",         ""));
    vh.push_back(PHS("DataMC_2015CD_mupt_0b",  
		     {PDC("data", "cutMuVec(pt)", {dsData_2015C_0b}), PDC("stack", "cutMuVec(pt)", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu1pt_0b",  
		     {PDC("data", "cutMuPt1", {dsData_2015C_0b}), PDC("stack", "cutMuPt1", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu1 pt",         ""));
    vh.push_back(PHS("DataMC_2015CD_mu2pt_0b",  
		     {PDC("data", "cutMuPt2", {dsData_2015C_0b}), PDC("stack", "cutMuPt2", {dsDY_0b, dstt2l_0b, dstW_0b, dsttZ_0b, dsVV_0b})},  
		     {1, 2}, "", 100, 0, 1000,   true, false,  "mu2 pt",         ""));


    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setPlotDir(plotDir);
    if(doSave)  plotter.saveHists();
    if(doPlots) plotter.plot();
}
