#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/searchBins.h"
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
        //{"savetuple",        no_argument, 0, 't'},
        {"fromFile",         no_argument, 0, 'f'},
        {"condor",           no_argument, 0, 'c'},
        {"histFile",   required_argument, 0, 'H'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
	{"luminosity", required_argument, 0, 'L'},
	{"sbEra",      required_argument, 0, 'S'}
    };

    bool doPlots = true, doSave = true, doTuple = true, fromTuple = true, runOnCondor = false;
    string histFile = "topStudyOutput.root", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "top_plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi = AnaSamples::luminosity;
    std::string sbEra = "SB_37_2015";

    while((opt = getopt_long(argc, argv, "psfcH:D:N:M:E:P:L:S:", long_options, &option_index)) != -1)
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

//        case 't':
//            if(doTuple) doPlots = doSave = false;
//            else        doTuple  = true;
//            break;

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

        case 'S':
            sbEra = optarg;
            break;
        }
    }

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "topStudyOutput_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        //doSave = true;
        //doPlots = false;
        fromTuple = true;
        sampleloc = "condor";
    }

    AnaSamples::SampleSet        ss(sampleloc, lumi);
    AnaSamples::SampleCollection sc(ss);

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_200to400"]};//{ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarInc"] = {ss["TTbarInc"]};
        fileMap["TTbar"] = sc["TTbar"];
        fileMap["Signal_T2tt_mStop500_mLSP325"] = sc["Signal_T2tt_mStop500_mLSP325"];
        fileMap["Signal_T1tttt_mGluino1200_mLSP800"] = sc["Signal_T1tttt_mGluino1200_mLSP800"];
        fileMap["Signal_T1tttt_mGluino1500_mLSP100"] = sc["Signal_T1tttt_mGluino1500_mLSP100"];
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

    Plotter::DatasetSummary ds_Znunu(   "Z#rightarrow#nu#nu", fileMap["ZJetsToNuNu"],                       "", "");
    Plotter::DatasetSummary ds_ttbar(   "t#bar{t}",           fileMap["TTbar"],                             "", "");
    Plotter::DatasetSummary ds_T2tt(    "T2tt (500, 325)",    fileMap["Signal_T2tt_mStop500_mLSP325"],      "", "");
    Plotter::DatasetSummary ds_T1tttt(  "T1tttt (1200, 800)", fileMap["Signal_T1tttt_mGluino1200_mLSP800"], "", "");
    Plotter::DatasetSummary ds_T1tttt_2("T1tttt (1500, 100)", fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "", "");
    
    auto PDCMaker = [&](std::string var) { return Plotter::DataCollection( "single", var, {ds_Znunu, ds_ttbar, ds_T2tt, ds_T1tttt, ds_T1tttt_2}); };

    auto PDCMaker_nTops = [&](Plotter::DatasetSummary& pds) { return Plotter::DataCollection( "single", {{"nTops", pds}, {"nTopsNew", pds}}); };
                                          
    vh.push_back(PHS("nTops",    {PDCMaker("nTops")},    {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTopsNew", {PDCMaker("nTopsNew")}, {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));

    vh.push_back(PHS("nTops_Znunu",    {PDCMaker_nTops(ds_Znunu)},    {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTops_ttbar",    {PDCMaker_nTops(ds_ttbar)},    {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTops_T2tt",     {PDCMaker_nTops(ds_T2tt)},     {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTops_T1tttt_1", {PDCMaker_nTops(ds_T1tttt)},   {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTops_T1tttt_2", {PDCMaker_nTops(ds_T1tttt_2)}, {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));

    vh.push_back(PHS("dTopsPt",    {PDCMaker("vdTops(pt)")},    {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  100, -25, 25,  true,  true,  "#Deltap_{T}(t)", "Norm Events"));
    vh.push_back(PHS("dTopsEta",   {PDCMaker("vdTops(eta)")},   {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  100, -5,   5,  true,  true,  "#Delta#eta(t)", "Norm Events"));
    vh.push_back(PHS("dTopsPhi",   {PDCMaker("vdTops(phi)")},   {1, 2}, "met>200;nTaggerJets>3;cntCSVS>0.5",  100, -5,   5,  true,  true,  "#Delta#phi(t)", "Norm Events"));

    vh.push_back(PHS("met", {PDCMaker("met")}, {1, 2}, "",  100, 0, 1000,  true,  true,  "MET [GeV]", "Norm Events"));

    vh.push_back(PHS("topMass1", {PDCMaker("vTops[0](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass2", {PDCMaker("vTops[1](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass3", {PDCMaker("vTops[2](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass4", {PDCMaker("vTops[3](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));

    vh.push_back(PHS("topMass1New", {PDCMaker("vTopsNew[0](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass2New", {PDCMaker("vTopsNew[1](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass3New", {PDCMaker("vTopsNew[2](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass4New", {PDCMaker("vTopsNew[3](M)")}, {1, 2}, "met>200",  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));

    vh.push_back(PHS("t3j_dr12", {PDCMaker("v3j_dr12[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));
    vh.push_back(PHS("t3j_dr23", {PDCMaker("v3j_dr23[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));
    vh.push_back(PHS("t3j_dr13", {PDCMaker("v3j_dr13[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));

    vh.push_back(PHS("t3j_cm_dr12", {PDCMaker("v3j_cm_dr12[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));
    vh.push_back(PHS("t3j_cm_dr23", {PDCMaker("v3j_cm_dr13[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));
    vh.push_back(PHS("t3j_dcm_r13", {PDCMaker("v3j_cm_dr23[0]")}, {1, 2}, "met>200",  100, 0, 10,  true,  true,  "#DeltaR", "Norm Events"));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsTopStudy();

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(false);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
