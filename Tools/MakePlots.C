#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/samples.h"

#include <getopt.h>
#include <iostream>

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

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
        {"filename",   required_argument, 0, 'I'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
        {"luminosity", required_argument, 0, 'L'},
        {"sbEra",      required_argument, 0, 'S'}
    };

    bool runOnCondor = false;
    bool doWeights = false;
    bool doPhotons = true;
    bool doSearchBins = false;
    bool doPlots = true;
    bool doSave = true;
    bool doTuple = true;
    bool fromTuple = true;
    string filename = "histoutput.root", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi = AnaSamples::luminosity;
    std::string sbEra = "SB_v1_2017";//"SB_v1_2017";

    while((opt = getopt_long(argc, argv, "pstfcI:D:N:M:E:P:L:S:", long_options, &option_index)) != -1)
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

        case 'I':
            filename = optarg;
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
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
        //filename = thistFile;
        //doSave = true;
        //doPlots = false;
        std::cout << "Filename modified for use with condor: " << filename << std::endl;
        fromTuple = true;
        sampleloc = "condor";
    }

    std::cout << "input filename: " << filename << std::endl;
    std::cout << "Sample location: " << sampleloc << std::endl;

    // follow this syntax; order matters for your arguments
    
    //SampleSet::SampleSet(std::string file, bool isCondor, double lumi)
    AnaSamples::SampleSet        ss("sampleSets.txt", runOnCondor, AnaSamples::luminosity);
    
    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples) : ss_(samples)
    AnaSamples::SampleCollection sc("sampleCollections.txt", ss);

    const double zAcc = 1.0;
    //    const double zAcc = 0.5954;
    //    const double zAcc = 0.855;
    const double znunu_mumu_ratio = 5.942;
    const double znunu_ee_ratio   = 5.942;

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_1200to2500"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["DYJetsToLL_HT_600to800"] = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["ZJetsToNuNu_HT_2500toInf"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["DYJetsToLL_HT_600to800"] = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["IncDY"] = {ss["DYJetsToLL_Inc"]}; 
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon"]};
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

    // Number of searchbins
    SearchBins sb(sbEra);
    int NSB = sb.nSearchBins();//37; // 45

    // Shortcuts for axis labels
    std::string label_met = "p_{T}^{miss} [GeV]";
    std::string label_ht  = "H_{T} [GeV]";
    std::string label_mht = "MH_{T} [GeV]";
    std::string label_nj  = "N_{jets}";
    std::string label_nb  = "N_{b}";
    std::string label_nt  = "N_{t}";
    std::string label_mt2 = "M_{T2} [GeV]";
    std::string label_eta = "#eta";
    std::string label_mupt = "#mu p_{T} [GeV]";
    std::string label_mueta = "#mu #eta";
    std::string label_genmupt  = "gen #mu p_{T} [GeV]";
    std::string label_genmueta = "gen #mu #eta";
    std::string label_mu1pt = "#mu_{1} p_{T} [GeV]";
    std::string label_mu2pt = "#mu_{2} p_{T} [GeV]";
    std::string label_mu1eta = "#mu_{1} #eta";
    std::string label_mu2eta = "#mu_{2} #eta";
    std::string label_elpt = "e p_{T} [GeV]";
    std::string label_eleta = "e #eta";
    std::string label_el1pt = "e_{1} p_{T} [GeV]";
    std::string label_el2pt = "e_{2} p_{T} [GeV]";
    std::string label_photonpt = "#gamma p_{T} [GeV]";
    std::string label_photoneta = "#gamma #eta";
    std::string label_jpt  = "j p_{T} [GeV]";
    std::string label_j1pt = "j_{1} p_{T} [GeV]";
    std::string label_j2pt = "j_{2} p_{T} [GeV]";
    std::string label_j3pt = "j_{3} p_{T} [GeV]";
    std::string label_mll  = "m_{ll} [GeV]";
    std::string label_topPt = "top p_{T} [GeV]";
    std::string label_genTopPt = "gen top p_{T} [GeV]";
    std::string label_phopt = "p_{T}^{#gamma} [GeV]";
    std::string label_metg = "p_{T}^{#gamma (miss)} [GeV]";
    std::string label_acc = "acc/gen";
    std::string label_reco = "reco/acc";
    std::string label_iso = "iso/reco";

    vector<Plotter::HistSummary> vh;

    Plotter::DatasetSummary dsDY_nunu(            "Z#rightarrow#nu#nu",                fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    // Datasetsummaries we are using                                                                                                        
    // no weight (genWeight deals with negative weights); also add btag weights here                                                        
    Plotter::DatasetSummary dsData_SingleMuon("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "");
    Plotter::DatasetSummary dsDY_mu(          "DY #mu",         fileMap["DYJetsToLL"],      "",        "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor");
    Plotter::DatasetSummary dsDYInc_mu(       "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor");
    Plotter::DatasetSummary dsDY_elec(        "DY e",         fileMap["DYJetsToLL"],      "",          "bTagSF_EventWeightSimple_Central;_PUweightFactor"); // do not use muTrigWgt for electrons (it is 0.0)
    Plotter::DatasetSummary dsDYInc_elec(     "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "bTagSF_EventWeightSimple_Central;_PUweightFactor"); // do not use muTrigWgt for electrons (it is 0.0)
    Plotter::DatasetSummary dsPhoton(         "#gamma+ jets", fileMap["GJets"],         "",            "bTagSF_EventWeightSimple_Central;_PUweightFactor");
    Plotter::DatasetSummary dstt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;isr_Unc_Cent;_PUweightFactor");
    Plotter::DatasetSummary dstW(             "Single t",   fileMap["SingleTopZinv"],   "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    Plotter::DatasetSummary dsttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    Plotter::DatasetSummary dsT1tttt_gluino1200_lsp800("T1tttt_gluino1200_lsp800",     fileMap["Signal_T1tttt_mGluino1200_mLSP800"], "",  "");
    Plotter::DatasetSummary dsT1tttt_gluino1500_lsp100("T1tttt_gluino1500_lsp100",     fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "",  "");
    Plotter::DatasetSummary dsT1tttt_gluino2000_lsp100("T1tttt_gluino2000_lsp100",     fileMap["Signal_T1tttt_mGluino2000_mLSP100"], "",  "");
    Plotter::DatasetSummary dsVV(             "Diboson",    fileMap["Diboson"],        "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    Plotter::DatasetSummary dsRare(           "Rare ",      fileMap["Rare"],           "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    std::vector<std::vector<Plotter::DatasetSummary>> stack_MC = {{dsDY_mu, dsDYInc_mu}, {dstt2l}, {dstW}, {dsRare, dsVV, dsttZ}};

    // Apply data/mc njet weight for DY and ttbar                                                                                                                                    
    Plotter::DatasetSummary dswDY(             "DY",         fileMap["DYJetsToLL"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;njWGJets;_PUweightFactor");
    Plotter::DatasetSummary dswDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "muTrigWgt;bTagSF_EventWeightSimple_Central;njWGJets;_PUweightFactor");
    Plotter::DatasetSummary dswtt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;nJetWgtTTbar;isr_Unc_Cent;_PUweightFactor");
    Plotter::DatasetSummary dswtW(             "Single t",   fileMap["SingleTopZinv"],   "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    Plotter::DatasetSummary dswttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    Plotter::DatasetSummary dswVV(             "Diboson",    fileMap["Diboson"],         "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    Plotter::DatasetSummary dswRare(           "Rare ",      fileMap["Rare"],            "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    std::vector<std::vector<Plotter::DatasetSummary>> stackw_MC = {{dswDY, dswDYInc}, {dswtt2l}, {dswtW}, {dswRare, dswVV, dswttZ}};

    Plotter::DatasetSummary dswwDY(             "DY",         fileMap["DYJetsToLL"],      "",            "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    Plotter::DatasetSummary dswwDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    std::vector<std::vector<Plotter::DatasetSummary>> stackww_MC = {{dswwDY, dswwDYInc}, {dswtt2l}, {dswtW}, {dswttZ}, {dswVV}, {dswRare, dswVV, dswttZ}};
 
    //lambda is your friend
    //for electrons do not use muTrigWgt (it is 0.0 for electrons)
    auto makePDSMu     = [&](const std::string& label) {return Plotter::DatasetSummary("DY "+label, fileMap["DYJetsToLL"], "", "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor"); };
    auto makePDSElec   = [&](const std::string& label) {return Plotter::DatasetSummary("DY "+label, fileMap["DYJetsToLL"], "", "bTagSF_EventWeightSimple_Central;_PUweightFactor"); };
    auto makePDSPhoton = [&](const std::string& label) {return Plotter::DatasetSummary("#gamma + jets "+label, fileMap["GJets"], "", "bTagSF_EventWeightSimple_Central;_PUweightFactor"); };
    
    // acceptance
    // muons
    Plotter::DataCollection dcMC_ngenMu(            "single", "ngenMu",                     {dsDY_mu});
    Plotter::DataCollection dcMC_ngenMatchMu(       "single", "ngenMatchMu",                {dsDY_mu});
    Plotter::DataCollection dcMC_ngenMuInAcc(       "single", "ngenMuInAcc",                {dsDY_mu});
    Plotter::DataCollection dcMC_ngenMatchMuInAcc(  "single", "ngenMatchMuInAcc",           {dsDY_mu});
    Plotter::DataCollection dcMC_genMuPt(           "single", "genMu(pt)",                  {dsDY_mu});
    Plotter::DataCollection dcMC_genMuInAccPt(      "single", "genMuInAcc(pt)",             {dsDY_mu});
    Plotter::DataCollection dcMC_genMatchMuInAccPt( "single", "genMatchMuInAcc(pt)",        {dsDY_mu});
    Plotter::DataCollection dcMC_genMuEta(          "single", "genMu(eta)",                 {dsDY_mu});
    Plotter::DataCollection dcMC_genMuInAccEta(     "single", "genMuInAcc(eta)",            {dsDY_mu});
    Plotter::DataCollection dcMC_genMatchMuInAccEta("single", "genMatchMuInAcc(eta)",       {dsDY_mu});
    // electrons
    Plotter::DataCollection dcMC_ngenElec(            "single", "ngenElec",                 {dsDY_elec});
    Plotter::DataCollection dcMC_ngenMatchElec(       "single", "ngenMatchElec",            {dsDY_elec});
    Plotter::DataCollection dcMC_ngenElecInAcc(       "single", "ngenElecInAcc",            {dsDY_elec});
    Plotter::DataCollection dcMC_ngenMatchElecInAcc(  "single", "ngenMatchElecInAcc",       {dsDY_elec});
    Plotter::DataCollection dcMC_genElecPt(           "single", "genElec(pt)",              {dsDY_elec});
    Plotter::DataCollection dcMC_genElecInAccPt(      "single", "genElecInAcc(pt)",         {dsDY_elec});
    Plotter::DataCollection dcMC_genMatchElecInAccPt( "single", "genMatchElecInAcc(pt)",    {dsDY_elec});
    Plotter::DataCollection dcMC_genElecEta(          "single", "genElec(eta)",             {dsDY_elec});
    Plotter::DataCollection dcMC_genElecInAccEta(     "single", "genElecInAcc(eta)",        {dsDY_elec});
    Plotter::DataCollection dcMC_genMatchElecInAccEta("single", "genMatchElecInAcc(eta)",   {dsDY_elec});
    // magic lambda functions... give it pt, eta, etc
    // leptons
    // acceptance = genInAcc / gen (acceptance / MC)
    auto makePDCElecAcc     = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genElecInAcc("+var+")", makePDSElec("e acc")}, {"genElec("+var+")", makePDSElec("e gen")}}); };
    auto makePDCMuAcc       = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genMuInAcc("+var+")",   makePDSMu("#mu acc")}, {"genMu("+var+")",   makePDSMu("#mu gen")}}); };
    // efficiency = genMatchInAcc / genInAcc (reco / acceptance)
    auto makePDCElecRecoEff = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genMatchElecInAcc("+var+")", makePDSElec("e reco")}, {"genElecInAcc("+var+")", makePDSElec("e acc")}}); };
    auto makePDCMuRecoEff   = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genMatchMuInAcc("+var+")",   makePDSMu("#mu reco")}, {"genMuInAcc("+var+")",   makePDSMu("#mu acc")}}); };
    // iso efficiency = genMatchIsoInAcc / genMatchInAcc (iso / reco)
    auto makePDCElecIsoEff  = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genMatchIsoElecInAcc("+var+")", makePDSElec("e iso")}, {"genMatchElecInAcc("+var+")", makePDSElec("e reco")}}); };
    auto makePDCMuIsoEff    = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"genMatchIsoMuInAcc("+var+")",   makePDSMu("#mu iso")}, {"genMatchMuInAcc("+var+")",   makePDSMu("#mu reco")}}); };
    // photons
    // acceptance = gammaLVecGenAcc / gen (acceptance / MC)
    auto makePDCPhotonAcc     = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"gammaLVecGenAcc("+var+")", makePDSPhoton("acc")}, {"gammaLVecGen("+var+")", makePDSPhoton("gen")}}); };
    // efficiency = promptPhotons / gammaLVecGenAcc (reco / acceptance)
    auto makePDCPhotonRecoEff = [&](const std::string& var, const std::string& style) {return Plotter::DataCollection(style, {{"promptPhotons("+var+")", makePDSPhoton("reco")}, {"gammaLVecGenAcc("+var+")", makePDSPhoton("acc")}}); };

    // photons
    Plotter::DataCollection dcMC_genPhotonPt(     "single", "gammaLVecGenAcc(pt)",      {dsPhoton});
    Plotter::DataCollection dcMC_genPhotonEta(    "single", "gammaLVecGenAcc(eta)",     {dsPhoton});
    Plotter::DataCollection dcMC_matchedPhotonPt( "single", "promptPhotons(pt)",        {dsPhoton});
    Plotter::DataCollection dcMC_matchedPhotonEta("single", "promptPhotons(eta)",       {dsPhoton});

    // tops
    Plotter::DataCollection dcMC_T1tttt("single",  "genTops(pt)", {dsT1tttt_gluino1200_lsp800, dsT1tttt_gluino1500_lsp100, dsT1tttt_gluino2000_lsp100});

    // nj                                                                                                                                                                            
    Plotter::DataCollection dcData_SingleMuon_nj("data",   "cntNJetsPt30Eta24Zinv", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_nj(             "stack",  "cntNJetsPt30Eta24Zinv", stack_MC);
    Plotter::DataCollection dcwMC_nj(            "stack",  "cntNJetsPt30Eta24Zinv", stackw_MC);
    Plotter::DataCollection dcwwMC_nj(           "stack",  "cntNJetsPt30Eta24Zinv", stackww_MC);
 
    // gen Z pt: genZPt
    Plotter::DataCollection dcData_DY_gen_pt("data",  "genZPt", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_DY_gen_pt("stack",   "genZPt", stack_MC);
  
    // reco Z pt: bestRecoZPt
    Plotter::DataCollection dcData_DY_reco_pt("data",  "bestRecoZPt", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_DY_reco_pt("stack",   "bestRecoZPt", stack_MC);

    // gen Z eta
    Plotter::DataCollection dcData_DY_eta("data",  "genZEta", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_DY_eta("stack",   "genZEta", stack_MC);

    // gen Z phi
    Plotter::DataCollection dcData_DY_phi("data",  "genZPhi", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_DY_phi("stack",   "genZPhi", stack_MC);

    // invariant mass
    Plotter::DataCollection dcData_DY_mass("data",  "genZmass", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_DY_mass("stack",   "genZmass", stack_MC);

    // met                                                                                                                                                                           
    Plotter::DataCollection dcData_SingleMuon_met("data",   "cleanMetPt", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_met(             "stack",  "cleanMetPt", stack_MC);
    Plotter::DataCollection dcwMC_met(            "stack",  "cleanMetPt", stackw_MC);
    Plotter::DataCollection dcwwMC_met(           "stack",  "cleanMetPt", stackww_MC);
    // ntops                                                                                                                                                                         
    Plotter::DataCollection dcData_SingleMuon_nt("data",   "nTopCandSortedCntZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_nt(             "stack",  "nTopCandSortedCntZinv", stack_MC);
    Plotter::DataCollection dcwMC_nt(            "stack",  "nTopCandSortedCntZinv", stackw_MC);
    Plotter::DataCollection dcwwMC_nt(           "stack",  "nTopCandSortedCntZinv", stackww_MC);
    // MT2                                                                                                                                                                           
    Plotter::DataCollection dcData_SingleMuon_mt2("data",   "best_had_brJet_MT2Zinv", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mt2(             "stack",  "best_had_brJet_MT2Zinv", stack_MC);
    Plotter::DataCollection dcwMC_mt2(            "stack",  "best_had_brJet_MT2Zinv", stackw_MC);
    Plotter::DataCollection dcwwMC_mt2(           "stack",  "best_had_brJet_MT2Zinv", stackww_MC);
    // nb                                                                                                                                  
    Plotter::DataCollection dcData_SingleMuon_nb("data",   "cntCSVSZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_nb(             "stack",  "cntCSVSZinv", stack_MC);
    Plotter::DataCollection dcwMC_nb(            "stack",  "cntCSVSZinv", stackw_MC);
    Plotter::DataCollection dcwwMC_nb(           "stack",  "cntCSVSZinv", stackww_MC);
    // ht                                                                                                                                                                            
    Plotter::DataCollection dcData_SingleMuon_ht("data",   "HTZinv", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_ht(             "stack",  "HTZinv", stack_MC);
    Plotter::DataCollection dcwMC_ht(            "stack",  "HTZinv", stackw_MC);
    Plotter::DataCollection dcwwMC_ht(           "stack",  "HTZinv", stackww_MC);
    // gen mu pt
    //Plotter::DataCollection dcData_SingleMuon_genMuPt("data",  "genMuPt", {dsData_SingleMuon});
    //Plotter::DataCollection dcMC_genMuPt("stack",   "genMuPt", stack_MC);
    // gen mu eta
    //Plotter::DataCollection dcData_SingleMuon_genMuEta("data",  "genMuEta", {dsData_SingleMuon});
    //Plotter::DataCollection dcMC_genMuEta("stack",   "genMuEta", stack_MC);
    // mu1pt                                                                                                                                                                         
    Plotter::DataCollection dcData_SingleMuon_mu1pt("data",   "cutMuPt1", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mu1pt(             "stack",  "cutMuPt1", stack_MC);
    Plotter::DataCollection dcwMC_mu1pt(            "stack",  "cutMuPt1", stackw_MC);
    // mu1eta                                                                                                                                                                         
    Plotter::DataCollection dcData_SingleMuon_mu1eta("data",   "cutMuEta1", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mu1eta(             "stack",  "cutMuEta1", stack_MC);
    Plotter::DataCollection dcwMC_mu1eta(            "stack",  "cutMuEta1", stackw_MC);
    // mu2pt                                                                                                                                                                         
    Plotter::DataCollection dcData_SingleMuon_mu2pt("data",   "cutMuPt2", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mu2pt(             "stack",  "cutMuPt2", stack_MC);
    Plotter::DataCollection dcwMC_mu2pt(            "stack",  "cutMuPt2", stackw_MC);
    // mu2eta                                                                                                                                                                         
    Plotter::DataCollection dcData_SingleMuon_mu2eta("data",   "cutMuEta2", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mu2eta(             "stack",  "cutMuEta2", stack_MC);
    Plotter::DataCollection dcwMC_mu2eta(            "stack",  "cutMuEta2", stackw_MC);
    // mll                                                                                                                                                                           
    Plotter::DataCollection dcData_SingleMuon_mll("data",   "bestRecoZM", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_mll(             "stack",  "bestRecoZM", stack_MC);
    Plotter::DataCollection dcwMC_mll(            "stack",  "bestRecoZM", stackw_MC);
    // nsearchbins                                                                                                                                                                   
    Plotter::DataCollection dcData_SingleMuon_nSearchBin("data",   "nSearchBin", {dsData_SingleMuon});
    Plotter::DataCollection dcMC_nSearchBin(             "stack",  "nSearchBin", stack_MC);
    Plotter::DataCollection dcwMC_nSearchBin(            "stack",  "nSearchBin", stackw_MC);
    Plotter::DataCollection dcwwMC_nSearchBin(           "stack",  "nSearchBin", stackww_MC);

    Plotter::DataCollection dcData_SingleMuon_mht("data",   "cleanMHt", {dsData_SingleMuon});
    Plotter::DataCollection dcwwMC_mht(           "stack",  "cleanMHt", stackww_MC);
    Plotter::DataCollection dcwMC_mht(            "stack",  "cleanMHt", stackw_MC);

    // electron cut levels
    std::vector<std::pair<std::string,std::string>> cutlevels_electrons = {
        //{"nothing",                   ""},
        {"nosel",                     "passNoiseEventFilterZinv"},
        //{"elecZinv",                  "passNoiseEventFilterZinv;passElecZinvSel"},
    };

    // muon cut levels; commented some cut leves to make less plots
    std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
        //{"nothing",                   ""},
        {"nosel",                     "passNoiseEventFilterZinv"},
        //{"muZinv",                    "passNoiseEventFilterZinv;passMuZinvSel"},
        //{"muZinv_blnotag",            "passMuZinvSel;passBaselineNoTagZinv"},
        //{"muZinv_bl",                 "passMuZinvSel;passBaselineZinv"},
        //{"muZinv_0b_blnotag",         "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
        //{"muZinv_g1b_blnotag",        "passMuZinvSel;passBJetsZinv;passBaselineNoTagZinv"},
        //{"muZinv_g1b",                "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv"},
        //{"muZinv_0b_loose0",          "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>300;passnJetsZinv;passdPhisZinv"},
        //{"muZinv_g1b_loose0",         "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"},
        //{"muZinv_loose0",             "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    };


    // loop over electron cut levels
    for(std::pair<std::string,std::string>& cut : cutlevels_electrons)
    {
        //individual plots
        //vh.push_back(PHS("MC_ngenElec_"                        +cut.first,  {dcMC_ngenElec},                                   {1, 1}, cut.second, 20, 0, 20,   true, false, "number of electrons",  "Events"));
        //vh.push_back(PHS("MC_ngenElecInAcc_"                   +cut.first,  {dcMC_ngenElecInAcc},                              {1, 1}, cut.second, 20, 0, 20,   true, false, "number of electrons",  "Events"));
        //vh.push_back(PHS("MC_ngenMatchElecInAcc_"              +cut.first,  {dcMC_ngenMatchElecInAcc},                         {1, 1}, cut.second, 20, 0, 20,   true, false, "number of electrons",  "Events"));
        //vh.push_back(PHS("MC_genElecPt_"                       +cut.first,  {dcMC_genElecPt},                                  {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt,  "Events"));
        //vh.push_back(PHS("MC_genElecInAccPt_"                  +cut.first,  {dcMC_genElecInAccPt},                             {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt,  "Events"));
        //vh.push_back(PHS("MC_genMatchElecInAccPt_"             +cut.first,  {dcMC_genMatchElecInAccPt},                        {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt,  "Events"));
        //vh.push_back(PHS("MC_genElecEta_"                      +cut.first,  {dcMC_genElecEta},                                 {1, 1}, cut.second, 40, -5, 5, true, false, label_eleta, "Events"));
        //vh.push_back(PHS("MC_genElecInAccEta_"                 +cut.first,  {dcMC_genElecInAccEta},                            {1, 1}, cut.second, 40, -5, 5, true, false, label_eleta, "Events"));
        //vh.push_back(PHS("MC_genMatchElecInAccEta_"            +cut.first,  {dcMC_genMatchElecInAccEta},                       {1, 1}, cut.second, 40, -5, 5, true, false, label_eleta, "Events"));

        //ratios
        
        // Acceptance
        vh.push_back(PHS("MC_ElecAccPt_ratio_"            +cut.first,  {makePDCElecAcc("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_elpt, label_acc));
        vh.push_back(PHS("MC_ElecAccEta_ratio_"           +cut.first,  {makePDCElecAcc("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, label_acc));
        vh.push_back(PHS("MC_ElecAccPt_single_"           +cut.first,  {makePDCElecAcc("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt, "Events"));
        vh.push_back(PHS("MC_ElecAccEta_single_"          +cut.first,  {makePDCElecAcc("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, "Events"));
        // Reco Efficiency
        vh.push_back(PHS("MC_ngenElecEff_original_"       +cut.first,  {dcMC_ngenMatchElecInAcc, dcMC_ngenElecInAcc},     {1, 2}, cut.second, 20, 0, 20, true, false, "number of electrons", "Events"));
        vh.push_back(PHS("MC_genElecPtEff_original_"      +cut.first,  {dcMC_genMatchElecInAccPt, dcMC_genElecInAccPt},   {1, 2}, cut.second, 60, 0, 1500, true, false, label_elpt, "Events"));
        vh.push_back(PHS("MC_genElecEtaEff_original_"     +cut.first,  {dcMC_genMatchElecInAccEta, dcMC_genElecInAccEta}, {1, 2}, cut.second, 40, -5, 5, true, false, label_eleta, "Events"));
        vh.push_back(PHS("MC_ElecRecoEffPt_ratio_"        +cut.first,  {makePDCElecRecoEff("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_elpt, label_reco));
        vh.push_back(PHS("MC_ElecRecoEffEta_ratio_"       +cut.first,  {makePDCElecRecoEff("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, label_reco));
        vh.push_back(PHS("MC_ElecRecoEffPt_single_"       +cut.first,  {makePDCElecRecoEff("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt, "Events"));
        vh.push_back(PHS("MC_ElecRecoEffEta_single_"      +cut.first,  {makePDCElecRecoEff("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, "Events"));
        // Iso Efficiency
        vh.push_back(PHS("MC_ElecIsoEffPt_ratio_"         +cut.first,  {makePDCElecIsoEff("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_elpt, label_iso));
        vh.push_back(PHS("MC_ElecIsoEffEta_ratio_"        +cut.first,  {makePDCElecIsoEff("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, label_iso));
        vh.push_back(PHS("MC_ElecIsoEffPt_single_"        +cut.first,  {makePDCElecIsoEff("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_elpt, "Events"));
        vh.push_back(PHS("MC_ElecIsoEffEta_single_"       +cut.first,  {makePDCElecIsoEff("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_eleta, "Events"));
    }

    // loop over muon cut levels
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
        //individual plots 
        //vh.push_back(PHS("MC_ngenMu_"                        +cut.first,  {dcMC_ngenMu},                                   {1, 1}, cut.second, 20, 0, 20,   true, false, "number of muons",  "Events"));
        //vh.push_back(PHS("MC_ngenMuInAcc_"                   +cut.first,  {dcMC_ngenMuInAcc},                              {1, 1}, cut.second, 20, 0, 20,   true, false, "number of muons",  "Events"));
        //vh.push_back(PHS("MC_ngenMatchMuInAcc_"              +cut.first,  {dcMC_ngenMatchMuInAcc},                         {1, 1}, cut.second, 20, 0, 20,   true, false, "number of muons",  "Events"));
        //vh.push_back(PHS("MC_genMuPt_"                       +cut.first,  {dcMC_genMuPt},                                  {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMuInAccPt_"                  +cut.first,  {dcMC_genMuInAccPt},                             {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMatchMuInAccPt_"             +cut.first,  {dcMC_genMatchMuInAccPt},                        {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMuEta_"                      +cut.first,  {dcMC_genMuEta},                                 {1, 1}, cut.second, 40, -5, 5, true, false, label_mueta, "Events"));
        //vh.push_back(PHS("MC_genMuInAccEta_"                 +cut.first,  {dcMC_genMuInAccEta},                            {1, 1}, cut.second, 40, -5, 5, true, false, label_mueta, "Events"));
        //vh.push_back(PHS("MC_genMatchMuInAccEta_"            +cut.first,  {dcMC_genMatchMuInAccEta},                       {1, 1}, cut.second, 40, -5, 5, true, false, label_mueta, "Events"));
        
        //ratios
        
        //vh.push_back(PHS("MC_ngenMuEff_original_"                   +cut.first,  {dcMC_ngenMatchMuInAcc, dcMC_ngenMuInAcc},       {1, 2}, cut.second, 20, 0, 20, true, false,  "number of muons",  "Events"));
        //vh.push_back(PHS("MC_genMuPtEff_origianl_"                  +cut.first,  {dcMC_genMatchMuInAccPt, dcMC_genMuInAccPt},     {1, 2}, cut.second, 60, 0, 1500, true, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMuEtaEff_original_"                 +cut.first,  {dcMC_genMatchMuInAccEta, dcMC_genMuInAccEta},   {1, 2}, cut.second, 40, -5, 5, true, false, label_mueta, "Events"));
        //vh.push_back(PHS("MC_genMuPtEff_ratio_"                     +cut.first,  {makePDCMuRecoEff("pt","ratio")},    {1, 1}, cut.second, 60, 0, 1500, false, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMuEtaEff_ratio_"                    +cut.first,  {makePDCMuRecoEff("eta","ratio")},   {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, "Events"));
        //vh.push_back(PHS("MC_genMuPtEff_single_"                    +cut.first,  {makePDCMuRecoEff("pt","single")},    {1, 1}, cut.second, 60, 0, 1500, false, false, label_mupt,  "Events"));
        //vh.push_back(PHS("MC_genMuEtaEff_single_"                   +cut.first,  {makePDCMuRecoEff("eta","single")},   {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, "Events"));
        
        // Acceptance
        vh.push_back(PHS("MC_MuAccPt_ratio_"            +cut.first,  {makePDCMuAcc("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_mupt, label_acc));
        vh.push_back(PHS("MC_MuAccEta_ratio_"           +cut.first,  {makePDCMuAcc("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, label_acc));
        vh.push_back(PHS("MC_MuAccPt_single_"           +cut.first,  {makePDCMuAcc("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt, "Events"));
        vh.push_back(PHS("MC_MuAccEta_single_"          +cut.first,  {makePDCMuAcc("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, "Events"));
        // Reco Efficiency
        vh.push_back(PHS("MC_ngenMuEff_original_"       +cut.first,  {dcMC_ngenMatchMuInAcc, dcMC_ngenMuInAcc},     {1, 2}, cut.second, 20, 0, 20, true, false,  "number of muons", "Events"));
        vh.push_back(PHS("MC_genMuPtEff_original_"      +cut.first,  {dcMC_genMatchMuInAccPt, dcMC_genMuInAccPt},   {1, 2}, cut.second, 60, 0, 1500, true, false, label_mupt, "Events"));
        vh.push_back(PHS("MC_genMuEtaEff_original_"     +cut.first,  {dcMC_genMatchMuInAccEta, dcMC_genMuInAccEta}, {1, 2}, cut.second, 40, -5, 5, true, false, label_mueta, "Events"));
        vh.push_back(PHS("MC_MuRecoEffPt_ratio_"        +cut.first,  {makePDCMuRecoEff("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_mupt, label_reco));
        vh.push_back(PHS("MC_MuRecoEffEta_ratio_"       +cut.first,  {makePDCMuRecoEff("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, label_reco));
        vh.push_back(PHS("MC_MuRecoEffPt_single_"       +cut.first,  {makePDCMuRecoEff("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt, "Events"));
        vh.push_back(PHS("MC_MuRecoEffEta_single_"      +cut.first,  {makePDCMuRecoEff("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, "Events"));
        // Iso Efficiency
        vh.push_back(PHS("MC_MuIsoEffPt_ratio_"         +cut.first,  {makePDCMuIsoEff("pt","ratio")},                {1, 1}, cut.second, 60, 0, 1500, false, false, label_mupt, label_iso));
        vh.push_back(PHS("MC_MuIsoEffEta_ratio_"        +cut.first,  {makePDCMuIsoEff("eta","ratio")},               {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, label_iso));
        vh.push_back(PHS("MC_MuIsoEffPt_single_"        +cut.first,  {makePDCMuIsoEff("pt","single")},               {1, 1}, cut.second, 60, 0, 1500, true, false, label_mupt, "Events"));
        vh.push_back(PHS("MC_MuIsoEffEta_single_"       +cut.first,  {makePDCMuIsoEff("eta","single")},              {1, 1}, cut.second, 40, -5, 5, false, false, label_mueta, "Events"));
        
        //Z things (DY)
        vh.push_back(PHS("DataMC_DY_gen_pt_"                 +cut.first,  {dcData_DY_gen_pt, dcMC_DY_gen_pt},              {1, 2}, cut.second, 60, 0, 1500, true, false,   "gen Z Pt",   "Events"));
        vh.push_back(PHS("DataMC_DY_reco_pt_"                +cut.first,  {dcData_DY_reco_pt, dcMC_DY_reco_pt},            {1, 2}, cut.second, 60, 0, 1500, true, false,   "reco Z Pt",  "Events"));
        vh.push_back(PHS("DataMC_DY_eta_"                    +cut.first,  {dcData_DY_eta, dcMC_DY_eta},                    {1, 2}, cut.second, 60, -10, 10, true, false,   "gen Z Eta",  "Events"));
        vh.push_back(PHS("DataMC_DY_phi_"                    +cut.first,  {dcData_DY_phi, dcMC_DY_phi},                    {1, 2}, cut.second, 60, -10, 10, true, false,   "gen Z Phi",  "Events"));
        vh.push_back(PHS("DataMC_DY_mass_"                   +cut.first,  {dcData_DY_mass, dcMC_DY_mass},                  {1, 2}, cut.second, 60, 0, 200, true, false,    "gen Z Mass", "Events"));
        if (doWeights)
        {
            //no weights
            vh.push_back(PHS("DataMC_SingleMuon_met_"            +cut.first,  {dcData_SingleMuon_met,   dcMC_met},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_met,      "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_ht_"             +cut.first,  {dcData_SingleMuon_ht,    dcMC_ht},              {1, 2}, cut.second, 60, 0, 1500, true, false,  label_ht,       "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nt_"             +cut.first,  {dcData_SingleMuon_nt,    dcMC_nt},              {1, 2}, cut.second, 5,  0, 5,    false, false,  label_nt,      "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mt2_"            +cut.first,  {dcData_SingleMuon_mt2,   dcMC_mt2},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mt2,      "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nb_"             +cut.first,  {dcData_SingleMuon_nb,    dcMC_nb},              {1, 2}, cut.second, 10, 0, 10,   true, false,  label_nb,       "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nj_"             +cut.first,  {dcData_SingleMuon_nj,    dcMC_nj},              {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,       "Events"));
            //vh.push_back(PHS("DataMC_SingleMuon_genMuPt_"        +cut.first,  {dcData_SingleMuon_genMuPt, dcMC_genMuPt},       {1, 2}, cut.second, 50, 0, 1000, true, false,  label_genmupt,  "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu1pt_"          +cut.first,  {dcData_SingleMuon_mu1pt, dcMC_mu1pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,    "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu2pt_"          +cut.first,  {dcData_SingleMuon_mu2pt, dcMC_mu2pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,    "Events"));
            ///vh.push_back(PHS("DataMC_SingleMuon_genMuEta_"       +cut.first,  {dcData_SingleMuon_genMuEta, dcMC_genMuEta},     {1, 2}, cut.second, 40, -10, 10, true, false,  label_genmueta, "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu1eta_"         +cut.first,  {dcData_SingleMuon_mu1eta, dcMC_mu1eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_mu1eta,   "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu2eta_"         +cut.first,  {dcData_SingleMuon_mu2eta, dcMC_mu2eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_mu2eta,   "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mll_"            +cut.first,  {dcData_SingleMuon_mll,   dcMC_mll},             {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll,      "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nSearchBin_"     +cut.first,  {dcData_SingleMuon_nSearchBin, dcMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin",   "Events"));
            //Shape correction weights
            vh.push_back(PHS("DataMC_SingleMuon_met_Wgt_"        +cut.first,  {dcData_SingleMuon_met,   dcwMC_met},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_met,  "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_ht_Wgt_"         +cut.first,  {dcData_SingleMuon_ht,    dcwMC_ht},              {1, 2}, cut.second, 60, 0, 1500, true, false,  label_ht,   "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nt_Wgt_"         +cut.first,  {dcData_SingleMuon_nt,    dcwMC_nt},              {1, 2}, cut.second, 5,  0, 5,    false, false,  label_nt,  "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mt2_Wgt_"        +cut.first,  {dcData_SingleMuon_mt2,   dcwMC_mt2},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mt2,  "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nb_Wgt_"         +cut.first,  {dcData_SingleMuon_nb,    dcwMC_nb},              {1, 2}, cut.second, 10, 0, 10,   true, false,  label_nb,   "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nj_Wgt_"         +cut.first,  {dcData_SingleMuon_nj,    dcwMC_nj},              {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,   "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu1pt_Wgt_"      +cut.first,  {dcData_SingleMuon_mu1pt, dcwMC_mu1pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt, "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu2pt_Wgt_"      +cut.first,  {dcData_SingleMuon_mu2pt, dcwMC_mu2pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt, "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu1eta_Wgt_"     +cut.first,  {dcData_SingleMuon_mu1eta, dcwMC_mu1eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_mu1eta, "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mu2eta_Wgt_"     +cut.first,  {dcData_SingleMuon_mu2eta, dcwMC_mu2eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_mu2eta, "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_mll_Wgt_"        +cut.first,  {dcData_SingleMuon_mll,   dcwMC_mll},             {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll,  "Events"));
            vh.push_back(PHS("DataMC_SingleMuon_nSearchBin_Wgt_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", "Events"));
        }
    }
    if (doWeights)
    {
        for(std::pair<std::string,std::string>& cut : cutlevels_muon)
        {
            vector<double> metBins = {0, 50, 100, 150, 200, 275, 300, 350, 400, 450, 2000};
            vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
            vh.push_back(PHS("DataMCw_SingleMuon_met_"        +cut.first,  {dcData_SingleMuon_met,        dcwMC_met},        {1, 2}, cut.second, 25, 0, 1500, true, false,  label_met,             "Events / 60 GeV"));
            vh.push_back(PHS("DataMCw_SingleMuon_rebin_met_"  +cut.first,  {dcData_SingleMuon_met,        dcwMC_met},        {1, 2}, cut.second, metBins,     true, false,  label_met,             "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_ht_"         +cut.first,  {dcData_SingleMuon_ht,         dcwMC_ht},         {1, 2}, cut.second, 60, 0, 1500, true, false,  label_ht,              "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_mht_"        +cut.first,  {dcData_SingleMuon_mht,        dcwMC_mht},        {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mht,             "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nt_"         +cut.first,  {dcData_SingleMuon_nt,         dcwMC_nt},         {1, 2}, cut.second, 5,  0, 5,    true, false,  label_nt,              "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_mt2_"        +cut.first,  {dcData_SingleMuon_mt2,        dcwMC_mt2},        {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mt2,             "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt2_"  +cut.first,  {dcData_SingleMuon_mt2,        dcwMC_mt2},        {1, 2}, cut.second, mt2Bins,     true, false,  label_mt2,             "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nb_"         +cut.first,  {dcData_SingleMuon_nb,         dcwMC_nb},         {1, 2}, cut.second, 7, 0, 7,   true, false,  label_nb,                "Events / bin"));
            vh.push_back(PHS("DataMCw_SingleMuon_nj_"         +cut.first,  {dcData_SingleMuon_nj,         dcwMC_nj},         {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,              "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin",          "Events"));

            //// Normalization weight applied, only for blnotag selections
            if(cut.first.rfind("blnotag") == (cut.first.size()-7))
            {
                // DataMC weights applied
                vh.push_back(PHS("DataMCww_SingleMuon_met_"   +cut.first,  {dcData_SingleMuon_met,   dcwwMC_met},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met,                                  "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_ht_"    +cut.first,  {dcData_SingleMuon_ht,    dcwwMC_ht},    {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,                                   "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_mht_"   +cut.first,  {dcData_SingleMuon_mht,   dcwwMC_mht},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht,                                  "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nt_"    +cut.first,  {dcData_SingleMuon_nt,    dcwwMC_nt},    {1, 2}, cut.second, 5,  0, 5,    true, false,  label_nt,                                   "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_mt2_"   +cut.first,  {dcData_SingleMuon_mt2,   dcwwMC_mt2},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2,                                  "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nb_"    +cut.first,  {dcData_SingleMuon_nb,    dcwwMC_nb},    {1, 2}, cut.second, 10, 0, 10,   true, false,  label_nb,                                   "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nj_"    +cut.first,  {dcData_SingleMuon_nj,    dcwwMC_nj},    {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,                                   "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin",               "Events"));
            }
        }
    }

    // photons
    if (doPhotons)
    {
        // photons: promptPhotons
        //vh.push_back(PHS("MC_ngenPhotonEff_"                     +cut.first,  {},       {1, 2}, cut.second, 20, 0, 20, true, false,  "number of photon",  "Events"));
        //vh.push_back(PHS("MC_genPhotonPtEff",   {dcMC_matchedPhotonPt, dcMC_genPhotonPt},     {1, 2}, "", 60, 0, 1500, true, false, label_photonpt,  "Events"));
        //vh.push_back(PHS("MC_genPhotonEtaEff",  {dcMC_matchedPhotonEta, dcMC_genPhotonEta},   {1, 2}, "", 40, -5, 5,   true, false, label_photoneta, "Events"));
        // Acceptance: use makePDCPhotonAcc
        vh.push_back(PHS("MC_genPhotonPtAcc_ratio",    {makePDCPhotonAcc("pt","ratio")},     {1, 1}, "", 60, 0, 1500, false, false, label_photonpt,  label_reco));
        vh.push_back(PHS("MC_genPhotonEtaAcc_ratio",   {makePDCPhotonAcc("eta","ratio")},    {1, 1}, "", 40, -5, 5,   false, false, label_photoneta, label_reco));
        vh.push_back(PHS("MC_genPhotonPtAcc_single",   {makePDCPhotonAcc("pt","single")},    {1, 1}, "", 60, 0, 1500, true,  false, label_photonpt,  "Events"));
        vh.push_back(PHS("MC_genPhotonEtaAcc_single",  {makePDCPhotonAcc("eta","single")},   {1, 1}, "", 40, -5, 5,   false, false, label_photoneta, "Events"));
        // Reco Efficiency
        vh.push_back(PHS("MC_genPhotonPtRecoEff_ratio",    {makePDCPhotonRecoEff("pt","ratio")},     {1, 1}, "", 60, 0, 1500, false, false, label_photonpt,  label_reco));
        vh.push_back(PHS("MC_genPhotonEtaRecoEff_ratio",   {makePDCPhotonRecoEff("eta","ratio")},    {1, 1}, "", 40, -5, 5,   false, false, label_photoneta, label_reco));
        vh.push_back(PHS("MC_genPhotonPtRecoEff_single",   {makePDCPhotonRecoEff("pt","single")},    {1, 1}, "", 60, 0, 1500, true,  false, label_photonpt,  "Events"));
        vh.push_back(PHS("MC_genPhotonEtaRecoEff_single",  {makePDCPhotonRecoEff("eta","single")},   {1, 1}, "", 40, -5, 5,   false, false, label_photoneta, "Events"));
    }

    //tops
    vh.push_back(PHS("DataMC_T1tttt", {dcMC_T1tttt}, {1, 1}, "passNoiseEventFilterZinv", 60, 0, 1500, true, false, label_genTopPt, "Events"));

    // nothing
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerCentral(          "Z#rightarrow#nu#nu Trigger Central ",        fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerUp(               "Z#rightarrow#nu#nu Trigger Up ",             fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerDown(             "Z#rightarrow#nu#nu Trigger Down ",           fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    Plotter::DatasetSummary dsDY_nunu_njetnorm(                         "Z#rightarrow#nu#nu Njet+norm ",              fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    // only b jet scale factor
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerCentral_scaled(   "Z#rightarrow#nu#nu Trigger scale Central ", fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerUp_scaled(        "Z#rightarrow#nu#nu Trigger scale Up ",      fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerDown_scaled(      "Z#rightarrow#nu#nu Trigger scale Down ",    fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_scaled(                  "Z#rightarrow#nu#nu Njet+norm scale ",       fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    // bjet scale factor, shape factor and normalization
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerCentral_weighted( "Z#rightarrow#nu#nu Trigger weight Central ", fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffMC");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerUp_weighted(      "Z#rightarrow#nu#nu Trigger weight Up ",      fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffUpMC");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerDown_weighted(    "Z#rightarrow#nu#nu Trigger weight Down ",    fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffDownMC");
    Plotter::DatasetSummary dsDY_nunu_njetnorm_weighted(                "Z#rightarrow#nu#nu Njet+norm weight ",       fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    Plotter::DataCollection trigger_nSearchBin( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown}, {"nSearchBin",    dsDY_nunu_njetnorm}  });
    Plotter::DataCollection trigger_nSearchBin_scaled( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_scaled}  });
    Plotter::DataCollection trigger_nSearchBin_weighted( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_weighted}  });


    if (doSearchBins)
    {
        vh.push_back(PHS("Trigger_",         {trigger_nSearchBin},           {2, 1}, "passBaseline",     NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
        vh.push_back(PHS("TriggerScl_",      {trigger_nSearchBin_scaled},    {2, 1}, "passBaseline",     NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
        vh.push_back(PHS("TriggerWgt_",      {trigger_nSearchBin_weighted},  {2, 1}, "passBaseline",     NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
        vh.push_back(PHS("Trigger_Zinv_",    {trigger_nSearchBin},           {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
        vh.push_back(PHS("TriggerScl_Zinv_", {trigger_nSearchBin_scaled},    {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
        vh.push_back(PHS("TriggerWgt_Zinv_", {trigger_nSearchBin_weighted},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
    }
    //Generate cutflows 
    vector<string> cfsZ = {"",
                           "passNoiseEventFilterZinv",
                           "passNoiseEventFilterZinv;passLeptVeto",
                           "passNoiseEventFilterZinv;passLeptVeto",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passBJetsZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passBJetsZinv;passTaggerZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passBJetsZinv;passTaggerZinv;passMETZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passTaggerZinv;passMETZinv;passBJetsZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
                           "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
                           "passLeptVeto;passBaselineZinv"};
    vector<Plotter::CutFlowSummary> cutFlowSummaries;

    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("ZtoNuNu",           PDC("", "", {dsDY_nunu}),           cfsZ));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, sbEra);

    std::cout << "Creating Plotter to make... well... you know... plots." << std::endl;
    Plotter plotter(vh, vvf, fromTuple, filename, nFiles, startFile, nEvts);
    plotter.setCutFlows(cutFlowSummaries);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
