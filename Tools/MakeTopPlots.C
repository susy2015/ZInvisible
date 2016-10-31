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
        fileMap["TTbarSingleLep"] = {ss["TTbarSingleLepT"]};
        //fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        //fileMap["TTbarInc"] = {ss["TTbarInc"]};
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
    Plotter::DatasetSummary ds_ttbar(   "t#bar{t} inc",       fileMap["TTbar"],                             "", "");
    Plotter::DatasetSummary ds_ttbar1l( "t#bar{t} l1",        fileMap["TTbarSingleLep"],                    "", "");
    Plotter::DatasetSummary ds_T2tt(    "T2tt (500, 325)",    fileMap["Signal_T2tt_mStop500_mLSP325"],      "", "");
    Plotter::DatasetSummary ds_T1tttt(  "T1tttt (1200, 800)", fileMap["Signal_T1tttt_mGluino1200_mLSP800"], "", "");
    Plotter::DatasetSummary ds_T1tttt_2("T1tttt (1500, 100)", fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "", "");

    auto PDSCutBind = [](Plotter::DatasetSummary pds, std::string cuts) { pds.setCuts(cuts); pds.label = pds.label + " " + cuts; return pds; };
    auto PDSLabelBind = [](Plotter::DatasetSummary pds, std::string label) { pds.label = pds.label + " " + label; return pds; };

    auto PDCMaker = [&](std::string var) { return Plotter::DataCollection( "single", var, {ds_T1tttt, ds_T1tttt_2, ds_T2tt, ds_Znunu, ds_ttbar, ds_ttbar1l}); };

    auto PDCMaker_nTops   = [&](Plotter::DatasetSummary& pds) { return Plotter::DataCollection( "single", {{"nTops", pds}, {"nTopsNew", pds}, {"nTopsNewMVA", pds}, {"vTopsMatchNewMVA(size)", pds}}); };
    auto vPDCMaker_disc    = [&](Plotter::DatasetSummary& pds) 
    {
        return std::vector<Plotter::DataCollection>( { Plotter::DataCollection( "single", {{"discriminators", pds}}),
                    Plotter::DataCollection( "stack",  {{"discriminatorsParMatch", PDSLabelBind(pds, "2 jet match")}, {"discriminatorsMatch", PDSLabelBind(pds, "3 jet match")} }),
                    } ); 
    };
    auto PDCMaker_topProp = [&](Plotter::DatasetSummary& pds, int iTop = 0, std::string var = "pt") 
    {
        std::string afterVar = ((iTop<0)?(""):("[" + std::to_string(iTop) + "]")) + "(" + var + ")";
        return Plotter::DataCollection( "single", { {"vTops"       + afterVar, PDSLabelBind(pds, "Orig")}, 
                                                    {"vTopsNew"    + afterVar, PDSLabelBind(pds, "New trijet only")},
                                                    {"vTopsNewMVA" + afterVar, PDSLabelBind(pds, "MVA trijet only")} });
    };
    auto vPDCMaker_topProp = [&](Plotter::DatasetSummary& pds, int iTop = 0, std::string var = "pt") 
    {
        std::string afterVar = ((iTop<0)?(""):("[" + std::to_string(iTop) + "]")) + "(" + var + ")";
        return std::vector<Plotter::DataCollection>({
                Plotter::DataCollection( "fill",   { {"vTopsNew"    + afterVar, PDSLabelBind(pds, "New trijet only")}, {"vTopsMatchNew"    + afterVar, PDSLabelBind(pds, "New trijet Match")}, {"vTopsParMatchNew"    + afterVar, PDSLabelBind(pds, "New trijet Match")} }),
                Plotter::DataCollection( "single", { {"vTopsNewMVA" + afterVar, PDSLabelBind(pds, "MVA trijet only")}, {"vTopsMatchNewMVA" + afterVar, PDSLabelBind(pds, "MVA trijet Match")}, {"vTopsParMatchNewMVA" + afterVar, PDSLabelBind(pds, "MVA trijet Match")} }),
                    });
    };
    auto PDCMaker_topMultComp = [&](Plotter::DatasetSummary& pds, std::string var = "pt")
    {
        std::string varStr = "vTopsNewMVA(" + var + ")";
        return Plotter::DataCollection( "single", {{varStr, PDSCutBind(pds, "nTopsNewMVA=0")}, {varStr, PDSCutBind(pds, "nTopsNewMVA=1")}, {varStr, PDSCutBind(pds, "nTopsNewMVA=2")}, {varStr, PDSCutBind(pds, "nTopsNewMVA>2")} });
    };
    auto vPDCMaker_Eff = [&](Plotter::DatasetSummary & pds, std::string var = "pt")
    {
        std::string varStr = "(" + var + ")";
        return std::vector<Plotter::DataCollection>({Plotter::DataCollection( "ratio", {{"vTopsGenMatchNew"    + varStr, PDSLabelBind(pds, "New trijet only")}, {"genTops" + varStr, PDSLabelBind(pds, "New trijet only")}}),
                                                     Plotter::DataCollection( "ratio", {{"vTopsGenMatchNewMVA" + varStr, PDSLabelBind(pds, "MVA trijet only")}, {"genTops" + varStr, PDSLabelBind(pds, "MVA trijet only")}})
                    });
    };
    auto vPDCMaker_Purity = [&](Plotter::DatasetSummary & pds, std::string var = "pt")
    {
        std::string varStr = "(" + var + ")";
        return std::vector<Plotter::DataCollection>({Plotter::DataCollection( "ratio", {{"vTopsMatchNew"    + varStr, PDSLabelBind(pds, "New trijet only")}, {"vTopsNew"    + varStr, PDSLabelBind(pds, "New trijet only")}}),
                                                     Plotter::DataCollection( "ratio", {{"vTopsMatchNewMVA" + varStr, PDSLabelBind(pds, "MVA trijet only")}, {"vTopsNewMva" + varStr, PDSLabelBind(pds, "MVA trijet only")}})
                    });
    };
    auto vPDCMaker_fake = [&](std::string var = "met")
    {
        return std::vector<Plotter::DataCollection>({Plotter::DataCollection( "ratio", {{var, PDSLabelBind(PDSCutBind(ds_Znunu, "nTopsNew>0.5"),    "New 3j")}, {var, PDSLabelBind(ds_Znunu, "New 3j")}}),
                                                     Plotter::DataCollection( "ratio", {{var, PDSLabelBind(PDSCutBind(ds_Znunu, "nTopsNewMVA>0.5"), "MVA 3j")}, {var, PDSLabelBind(ds_Znunu, "MVA 3j")}})
                    });
    };

    std::string cuts = "passLeptVeto;passBJets;passnJets;passMET;nTaggerJets>4";                                          
    vh.push_back(PHS("nTops",    {PDCMaker("nTops")},       {1, 2}, cuts,  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));
    vh.push_back(PHS("nTopsNew", {PDCMaker("nTopsNewMVA")}, {1, 2}, cuts,  10, 0, 10,  true,  true,  "N_{T}", "Norm Events"));

    std::vector<std::pair<PDS, std::string>> pdsVec = {{ds_T1tttt, "T1tttt_1"}, {ds_T1tttt_2, "T1tttt_2"}, {ds_T2tt, "T2tt"}, {ds_Znunu, "Znunu"}, {ds_ttbar, "ttbar"}, {ds_ttbar1l, "ttbar1l"}};
    for(auto& pds : pdsVec)
    {
        vh.push_back(PHS("nTops_"          + pds.second,  {PDCMaker_topProp(pds.first, -1, "size")},  {1, 2}, cuts + "",  10,   0, 10,    false,  false,  "N_{T}",             "Events"));
        vh.push_back(PHS("disc_"           + pds.second,  vPDCMaker_disc(pds.first),                  {1, 2}, cuts + "",   50,  0, 1,     false,  false,  "disc",              "Events"));
        vh.push_back(PHS("tAll_pt_"        + pds.second,  {PDCMaker_topProp(pds.first, -1, "pt")},    {1, 2}, cuts + "",   50,  0, 2000,  false,  false,  "p_{T} top [GeV]",   "Events"));
        vh.push_back(PHS("t1_pt_"          + pds.second,  {PDCMaker_topProp(pds.first,  0, "pt")},    {1, 2}, cuts + "",   50,  0, 2000,  false,  false,  "p_{T} top 1 [GeV]", "Events"));
        vh.push_back(PHS("t2_pt_"          + pds.second,  {PDCMaker_topProp(pds.first,  1, "pt")},    {1, 2}, cuts + "",   50,  0, 2000,  false,  false,  "p_{T} top 2 [GeV]", "Events"));
        vh.push_back(PHS("t3_pt_"          + pds.second,  {PDCMaker_topProp(pds.first,  2, "pt")},    {1, 2}, cuts + "",   50,  0, 2000,  false,  false,  "p_{T} top 3 [GeV]", "Events"));
        vh.push_back(PHS("t4_pt_"          + pds.second,  {PDCMaker_topProp(pds.first,  3, "pt")},    {1, 2}, cuts + "",   50,  0, 2000,  false,  false,  "p_{T} top 4 [GeV]", "Events"));
        vh.push_back(PHS("tAll_eta_"       + pds.second,  {PDCMaker_topProp(pds.first, -1, "eta")},   {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top",           "Events"));
        vh.push_back(PHS("t1_eta_"         + pds.second,  {PDCMaker_topProp(pds.first,  0, "eta")},   {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top 1",         "Events"));
        vh.push_back(PHS("t2_eta_"         + pds.second,  {PDCMaker_topProp(pds.first,  1, "eta")},   {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top 2",         "Events"));
        vh.push_back(PHS("t3_eta_"         + pds.second,  {PDCMaker_topProp(pds.first,  2, "eta")},   {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top 3",         "Events"));
        vh.push_back(PHS("t4_eta_"         + pds.second,  {PDCMaker_topProp(pds.first,  3, "eta")},   {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top 4",         "Events"));
        vh.push_back(PHS("match_nTops_"    + pds.second,  vPDCMaker_topProp(pds.first, -1, "size"),   {1, 1}, cuts + "",  10,   0, 10,    false,  false,  "N_{T}",             "Events"));
        vh.push_back(PHS("match_tAll_pt_"  + pds.second,  vPDCMaker_topProp(pds.first, -1, "pt"),     {1, 1}, cuts + "",  100,  0, 2000,  false,  false,  "p_{T} top [GeV]",   "Events"));
        vh.push_back(PHS("match_tAll_eta_" + pds.second,  vPDCMaker_topProp(pds.first, -1, "eta"),    {1, 1}, cuts + "",  100, -5, 5,     false,  false,  "eta top",           "Events"));
        vh.push_back(PHS("tops_pt_"        + pds.second,  {PDCMaker_topMultComp(pds.first,  "pt")},   {1, 2}, cuts + "",  100, -5, 2000,  false,  false,  "p_{T} top [GeV]",   "Events"));
        vh.push_back(PHS("tops_m_"         + pds.second,  {PDCMaker_topMultComp(pds.first,  "M")},    {1, 2}, cuts + "",  100, -5, 500,   false,  false,  "mass top [GeV]",    "Events"));
        vh.push_back(PHS("tops_eta_"       + pds.second,  {PDCMaker_topMultComp(pds.first,  "eta")},  {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top",           "Events"));
        vh.push_back(PHS("tops_eta_"       + pds.second,  {PDCMaker_topMultComp(pds.first,  "eta")},  {1, 2}, cuts + "",  100, -5, 5,     false,  false,  "eta top",           "Events"));
        
        vh.push_back(PHS("eff_eta_" + pds.second,  vPDCMaker_Eff(pds.first, "eta"),  {1, 1}, cuts + "",  100, -5,    5,     false,  false,  "top eta",          "Efficiency"));
        vh.push_back(PHS("eff_pt_"  + pds.second,  vPDCMaker_Eff(pds.first, "pt"),   {1, 1}, cuts + "",  100,  0, 2000,     false,  false,  "top pt [GeV]",     "Efficiency"));
        vh.push_back(PHS("eff_m_"   + pds.second,  vPDCMaker_Eff(pds.first, "M"),    {1, 1}, cuts + "",  100,  0,  500,     false,  false,  "top mass [GeV]",   "Efficiency"));

        vh.push_back(PHS("purity_eta_" + pds.second,  vPDCMaker_Purity(pds.first, "eta"),  {1, 1}, cuts + "",  100, -5,    5,     false,  false,  "top eta",          "Purity"));
        vh.push_back(PHS("purity_pt_"  + pds.second,  vPDCMaker_Purity(pds.first, "pt"),   {1, 1}, cuts + "",  100,  0, 2000,     false,  false,  "top pt [GeV]",     "Purity"));
        vh.push_back(PHS("purity_m_"   + pds.second,  vPDCMaker_Purity(pds.first, "M"),    {1, 1}, cuts + "",  100,  0,  500,     false,  false,  "top mass [GeV]",   "Purity"));
    }

    vh.push_back(PHS("fakerate_met",   vPDCMaker_fake("met"),                {1, 1}, cuts + "",  100,  0,  2000,     false,  false,  "met [GeV]",   "Fakerate"));
    vh.push_back(PHS("fakerate_njet",  vPDCMaker_fake("cntNJetsPt50Eta24"),  {1, 1}, cuts + "",   20,  0,    20,     false,  false,  "N Jets",      "Fakerate"));
    vh.push_back(PHS("fakerate_nbjet", vPDCMaker_fake("nbjets"),             {1, 1}, cuts + "",   10,  0,    10,     false,  false,  "N b-jets",    "Fakerate"));

    vh.push_back(PHS("met", {PDCMaker("met")}, {1, 2}, "",  100, 0, 1000,  true,  true,  "MET [GeV]", "Norm Events"));

    vh.push_back(PHS("topMass1", {PDCMaker("vTops[0](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass2", {PDCMaker("vTops[1](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass3", {PDCMaker("vTops[2](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass4", {PDCMaker("vTops[3](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));

    vh.push_back(PHS("topMassAllNew", {PDCMaker("vTopsNewMVA(M)")},    {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass1New",   {PDCMaker("vTopsNewMVA[0](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass2New",   {PDCMaker("vTopsNewMVA[1](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass3New",   {PDCMaker("vTopsNewMVA[2](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));
    vh.push_back(PHS("topMass4New",   {PDCMaker("vTopsNewMVA[3](M)")}, {1, 2}, cuts,  100, 0, 500,  true,  true,  "Top Mass [GeV]", "Norm Events"));

    vh.push_back(PHS("topPtAllNew", {PDCMaker("vTopsNewMVA(pt)")},    {1, 2}, cuts,  100, 0, 2000,  true,  true,  "Top p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("topPt1New",   {PDCMaker("vTopsNewMVA[0](pt)")}, {1, 2}, cuts,  100, 0, 2000,  true,  true,  "Top p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("topPt2New",   {PDCMaker("vTopsNewMVA[1](pt)")}, {1, 2}, cuts,  100, 0, 2000,  true,  true,  "Top p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("topPt3New",   {PDCMaker("vTopsNewMVA[2](pt)")}, {1, 2}, cuts,  100, 0, 2000,  true,  true,  "Top p_{T} [GeV]", "Norm Events"));
    vh.push_back(PHS("topPt4New",   {PDCMaker("vTopsNewMVA[3](pt)")}, {1, 2}, cuts,  100, 0, 2000,  true,  true,  "Top p_{T} [GeV]", "Norm Events"));

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
