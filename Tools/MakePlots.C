#include "Plotter.h"
#include "samples.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SusyUtility.h"
#include "TMath.h"

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
        {"dophotons",        no_argument, 0, 'g'},
        {"doleptons",        no_argument, 0, 'l'},
        {"verbose",          no_argument, 0, 'v'},
        {"filename",   required_argument, 0, 'I'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
        {"luminosity", required_argument, 0, 'L'},
        {"refLumi",    required_argument, 0, 'R'},
        {"sbEra",      required_argument, 0, 'S'},
        {"era",        required_argument, 0, 'Y'}
    };
    bool runOnCondor        = false;
    bool doDataMCElectron   = true;
    bool doDataMCMuon       = true;
    bool doDataMCPhoton     = true;
    bool doWeights = false;
    bool doLeptons = false;
    bool doPhotons = false;
    bool doGJetsAndZnunu = false;
    bool doDYAndZnunu = false;
    bool doSearchBins = true;
    bool doPlots = true;
    bool doSave = true;
    bool doTuple = true;
    bool fromTuple = true;
    bool verbose = false;
    int startFile   = 0;
    int nFiles      = -1;
    int nEvts       = -1;
    double lumi     = -1.0;
    std::string refLumi = "";
    std::string filename = "histoutput.root";
    std::string dataSets = "";
    std::string sampleloc = AnaSamples::fileDir;
    std::string plotDir = "plots";
    std::string sbEra = "SB_v1_2017";
    std::string era  = "";
    std::string year = "";
    while((opt = getopt_long(argc, argv, "pstfcglvI:D:N:M:E:P:L:R:S:Y:", long_options, &option_index)) != -1)
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
        
        case 'g':
            doPhotons = true;
            break;

        case 'l':
            doLeptons = true;
            break;

        case 'v':
            verbose = true;
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
        
        case 'R':
            refLumi = optarg;
            break;

        case 'S':
            sbEra = optarg;
            break;
        
        case 'Y':
            era = optarg;
            break;
        }
    }

    // get year from era
    // if era is 2018_PreHEM, year is 2018
    if (era.length() >= 4)
    {
        year = era.substr(0, 4);
    }
    
    // datasets
    std::string ElectronDataset = "Data_SingleElectron";
    std::string MuonDataset     = "Data_SingleMuon";
    std::string PhotonDataset   = "Data_SinglePhoton";
    // year and periods
    std::string eraTag    = "_" + era; 
    std::string yearTag   = "_" + year; 
    std::string periodTag = ""; 
    // HEM veto for 2018 PostHEM
    bool doHEMVeto = false;
    std::string HEMVeto                 = "";
    std::string HEMVeto_drLeptonCleaned = "";
    std::string HEMVeto_drPhotonCleaned = "";
    std::string semicolon_HEMVeto                 = "";
    std::string semicolon_HEMVeto_drLeptonCleaned = "";
    std::string semicolon_HEMVeto_drPhotonCleaned = "";
    // Note: Don't apply Flag_ecalBadCalibFilter to 2016, but apply it to 2017 and 2018
    std::string Flag_ecalBadCalibFilter = "";
    std::string PrefireWeight           = "";
    std::string puWeight                = ";puWeight";
    std::string SAT_Pass_lowDM          = "";
    std::string SAT_Pass_highDM         = "";
    std::string SAT_Pass_lowDM_Loose    = "";
    std::string SAT_Pass_highDM_Loose   = "";
    std::string SAT_Pass_lowDM_Mid      = "";
    std::string SAT_Pass_highDM_Mid     = "";

    // lumi for Plotter
    if (era.compare("2016") == 0)
    {
        PrefireWeight   = ";PrefireWeight";
    }
    else if (era.compare("2017") == 0)
    {
        PrefireWeight           = ";PrefireWeight";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2017_BE") == 0)
    {
        puWeight                = ";17BtoEpuWeight";
        PrefireWeight           = ";PrefireWeight";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2017_F") == 0)
    {
        puWeight                = ";17FpuWeight";
        PrefireWeight           = ";PrefireWeight";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2018") == 0)
    {
        ElectronDataset         = "Data_EGamma";
        PhotonDataset           = "Data_EGamma";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2018_PreHEM") == 0)
    {
        periodTag               = "_PreHEM";
        ElectronDataset         = "Data_EGamma";
        PhotonDataset           = "Data_EGamma";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2018_PostHEM") == 0)
    {
        // HEM vetos: use ";veto_name" so that it can be appended to cuts
        periodTag                         = "_PostHEM";
        ElectronDataset                   = "Data_EGamma";
        PhotonDataset                     = "Data_EGamma";
        Flag_ecalBadCalibFilter           = ";Flag_ecalBadCalibFilter";
        HEMVeto                           = "SAT_Pass_HEMVeto";
        HEMVeto_drLeptonCleaned           = "SAT_Pass_HEMVeto_drLeptonCleaned";
        HEMVeto_drPhotonCleaned           = "SAT_Pass_HEMVeto_drPhotonCleaned";
        semicolon_HEMVeto                 = ";" + HEMVeto;
        semicolon_HEMVeto_drLeptonCleaned = ";" + HEMVeto_drLeptonCleaned;
        semicolon_HEMVeto_drPhotonCleaned = ";" + HEMVeto_drPhotonCleaned;
        doHEMVeto                         = true;
    }
    else
    {
        std::cout << "Please enter 2016, 2017, 2017_BE, 2017_F, 2018, 2018_PreHEM or 2018_PostHEM for the era using the -Y flag." << std::endl;
        exit(1);
    }

    // if luminosity not set with -L flag
    if (lumi < 0.0)
    {
        if (refLumi.empty())
        {
            std::cout << "Please enter either a luminosity (-L flag) or a data sample collection name (-R flag) for use as a reference luminosity." << std::endl;
            exit(1);
        }
        else
        {
            // Hack to get luminosity from reference
            AnaSamples::SampleSet        SS_temp("sampleSets_PostProcessed_"+year+".cfg", runOnCondor, lumi);
            AnaSamples::SampleCollection SC_temp("sampleCollections_"+year+".cfg", SS_temp);
            lumi = SC_temp.getSampleLumi(refLumi);
        }
    }
    
    std::cout << "Using luminosity " << lumi << std::endl;
    
    // add year and period tags
    ElectronDataset = ElectronDataset + yearTag + periodTag;
    MuonDataset     = MuonDataset     + yearTag + periodTag;
    PhotonDataset   = PhotonDataset   + yearTag + periodTag;
    // testing
    //printf("ElectronDataset: %s\n", ElectronDataset.c_str());
    //printf("MuonDataset: %s\n",     MuonDataset.c_str());
    //printf("PhotonDataset: %s\n",   PhotonDataset.c_str());

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
        std::cout << "Filename modified for use with condor: " << filename << std::endl;
        //doSave = true;
        //doPlots = false;
        fromTuple = true;
        doPhotons = true;
        doLeptons = false;
        sampleloc = "condor";
    }

    struct sampleStruct
    {
        AnaSamples::SampleSet           sample_set;
        AnaSamples::SampleCollection    sample_collection;
        std::string                     sample_year;
    };

    // --- follow the syntax; order matters for your arguments --- //
    
    //SampleSet::SampleSet(std::string file, bool isCondor, double lumi)
    AnaSamples::SampleSet        SS_2016("sampleSets_PostProcessed_2016.cfg", runOnCondor, lumi);
    AnaSamples::SampleSet        SS_2017("sampleSets_PostProcessed_2017.cfg", runOnCondor, lumi);
    AnaSamples::SampleSet        SS_2018("sampleSets_PostProcessed_2018.cfg", runOnCondor, lumi);
    
    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples)
    AnaSamples::SampleCollection SC_2016("sampleCollections_2016.cfg", SS_2016);
    AnaSamples::SampleCollection SC_2017("sampleCollections_2017.cfg", SS_2017);
    AnaSamples::SampleCollection SC_2018("sampleCollections_2018.cfg", SS_2018);

    // Warning: keep years together when you add them to sampleList:
    std::vector<sampleStruct> sampleList;
    sampleList.push_back({SS_2016, SC_2016, "2016"});
    sampleList.push_back({SS_2017, SC_2017, "2017"});
    sampleList.push_back({SS_2018, SC_2018, "2018"});
    
    const double zAcc = 1.0;
    // const double zAcc = 0.5954;
    // const double zAcc = 0.855;
    const double znunu_mumu_ratio = 5.942;
    const double znunu_ee_ratio   = 5.942;

    map<string, vector<AFS>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL" + yearTag]                = {SS_2016["DYJetsToLL_HT_1200to2500" + yearTag]};
        fileMap["ZJetsToNuNu" + yearTag]               = {SS_2016["ZJetsToNuNu_HT_2500toInf" + yearTag]};
        fileMap["DYJetsToLL_HT_600to800" + yearTag]    = {SS_2016["DYJetsToLL_HT_600to800" + yearTag]};
        fileMap["ZJetsToNuNu_HT_2500toInf" + yearTag]  = {SS_2016["ZJetsToNuNu_HT_2500toInf" + yearTag]};
        fileMap["TTbarDiLep" + yearTag]                = {SS_2016["TTbarDiLep" + yearTag]};
        fileMap["TTbarNoHad" + yearTag]                = {SS_2016["TTbarDiLep" + yearTag]};
        fileMap[MuonDataset]                           = {SS_2016[MuonDataset]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL" + yearTag]              = {SS_2016["DYJetsToLL_HT_600to800" + yearTag]};
        fileMap["DYJetsToLL_HT_600to800" + yearTag]  = {SS_2016["DYJetsToLL_HT_600to800" + yearTag]};
        fileMap["IncDY" + yearTag]                   = {SS_2016["DYJetsToLL_Inc" + yearTag]}; 
        fileMap["TTbarDiLep" + yearTag]              = {SS_2016["TTbarDiLep" + yearTag]};
        fileMap["TTbarNoHad" + yearTag]              = {SS_2016["TTbarDiLep" + yearTag]};
        fileMap[MuonDataset]                         = {SS_2016[MuonDataset]};
    }
    else
    {
        for (const auto& sample : sampleList)
        {
            AnaSamples::SampleSet           ss = sample.sample_set;
            AnaSamples::SampleCollection    sc = sample.sample_collection;
            std::string                     sy = sample.sample_year; 
            // --- calculate total luminosity for data --- //
            //printf("year: %s\n", sy.c_str());
            //printf("%s: lumi = %f\n", (ElectronDataset).c_str(),    sc.getSampleLumi(ElectronDataset));
            //printf("%s: lumi = %f\n", (MuonDataset).c_str(),        sc.getSampleLumi(MuonDataset));
            //printf("%s: lumi = %f\n", (PhotonDataset).c_str(),      sc.getSampleLumi(PhotonDataset));
            // ------------------------------------------- // 
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
    }


    //SearchBins sb(sbEra);
    // Number of searchbins
    int NSB = 204;
    // search bins for low and high dm
    // Low DM, 53 bins: 0 - 52
    // High DM, 151 bins: 53 - 203
    // Total 204 bins: 0 - 203
    int min_sb_low_dm = 0;
    int max_sb_low_dm = 53;
    int min_sb_high_dm = 53;
    int max_sb_high_dm = 204;
    //Validation Bins
    // Low DM, 15 bins: 0 - 14
    // Low DM High MET, 4 bins: 15 - 18
    // High DM, 24 bins: 22 - 45
    // Total 43 bins: 0 - 18 and 22 - 45
    int min_vb_low_dm = 0;
    int max_vb_low_dm = 15;
    int min_vb_low_dm_high_met = 15;
    int max_vb_low_dm_high_met = 19;
    int min_vb_high_dm = 22;
    int max_vb_high_dm = 46;

    // min and max values for histos
    int nBins = 40;
    // p_t in GeV
    double minPt = 0.0;
    double maxPt = 1000.0;
    // Energy in GeV
    double minEnergy = 0.0;
    double maxEnergy = 2000.0;
    int minJets = 0;
    int maxJets = 20;
    // mass in GeV
    // mass of electron: 0.511 MeV = 5.11 * 10^-4 GeV
    // mass of muon: 106 MeV = 0.106 GeV
    // mass of photon: 0.00 eV
    double minMassElec = 0.0;
    double maxMassElec = TMath::Power(10.0, -3);
    double minMassMu = 0.0;
    double maxMassMu = 0.2;
    double minMassPhoton = -2.0 * TMath::Power(10.0, -4);
    double maxMassPhoton =  2.0 * TMath::Power(10.0, -4);
    double minEta = -5.0;
    double maxEta = 5.0;
    double minPhi = -1.0 * TMath::Pi();
    double maxPhi = TMath::Pi();

    // met bin edges
    std::vector<double> metBinEdges      = {250.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0};
    std::vector<double> photonPtBinEdges = {220.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0};

    // Shortcuts for axis labels
    std::string label_Events = "Events";
    //std::string label_met = "p_{T}^{miss} [GeV]";
    std::string label_met = "#slash{E}_{T} [GeV]";
    std::string label_metWithLL = "#slash{E}_{T}^{LL} [GeV]";
    std::string label_metWithPhoton = "#slash{E}_{T}^{#gamma} [GeV]";
    std::string label_metphi = "#phi_{MET}";
    std::string label_metphiWithLL = "#phi_{MET}^{LL}";
    std::string label_metphiWithPhoton = "#phi_{MET}^{#gamma}";
    std::string label_bestRecoZM  = "m_{LL} [GeV]";
    std::string label_bestRecoZPt = "p_{T}(LL) [GeV]";
    std::string label_ht  = "H_{T} [GeV]";
    std::string label_mtb = "M_{T}(b_{1,2}, #slash{E}_{T}) [GeV]";
    std::string label_ptb = "p_{T}(b) [GeV]";
    std::string label_ISRJetPt = "ISR Jet p_{T} [GeV]"; 
    std::string label_mht = "MH_{T} [GeV]";
    std::string label_nj  = "N_{jets}";
    std::string label_nb  = "N_{bottoms}";
    std::string label_nt  = "N_{tops}";
    std::string label_dr  = "#DeltaR";
    // start with dphi1 for leading jet 1
    std::string label_dphi1  = "#Delta#phi_{1}";
    std::string label_dphi2  = "#Delta#phi_{2}";
    std::string label_dphi3  = "#Delta#phi_{3}";
    std::string label_dphi4  = "#Delta#phi_{4}";
    std::string label_dphi5  = "#Delta#phi_{5}";
    std::vector<std::string> vec_label_dphi = {label_dphi1, label_dphi2, label_dphi3, label_dphi4, label_dphi5}; 
    std::string label_mt2 = "M_{T2} [GeV]";
    std::string label_eta = "#eta";
    std::string label_MuPt = "p_{T}^{#mu} [GeV]";
    std::string label_MuEnergy = "E^{#mu} [GeV]";
    std::string label_MuMass = "m^{#mu} [GeV]";
    std::string label_MuEta = "#eta^{#mu}";
    std::string label_MuPhi = "#phi^{#mu}";
    std::string label_genmupt  = "gen #mu p_{T} [GeV]";
    std::string label_genmueta = "gen #mu #eta";
    std::string label_MuPt1 = "#mu_{1} p_{T} [GeV]";
    std::string label_MuPt2 = "#mu_{2} p_{T} [GeV]";
    std::string label_MuEta1 = "#mu_{1} #eta";
    std::string label_MuEta2 = "#mu_{2} #eta";
    std::string label_ElecPt = "p_{T}^{e} [GeV]";
    std::string label_ElecEnergy = "E^{e} [GeV]";
    std::string label_ElecMass = "m^{e} [GeV]";
    std::string label_ElecEta = "#eta^{e}";
    std::string label_ElecPhi = "#phi^{e}";
    std::string label_ElecPt1 = "e_{1} p_{T} [GeV]";
    std::string label_ElecPt2 = "e_{2} p_{T} [GeV]";
    std::string label_ElecEta1 = "e_{1} #eta";
    std::string label_ElecEta2 = "e_{2} #eta";
    std::string label_PhotonPt = "p_{T}^{#gamma} [GeV]";
    std::string label_PhotonEnergy = "E^{#gamma} [GeV]";
    std::string label_PhotonMass = "m^{#gamma} [GeV]";
    std::string label_PhotonEta = "#eta^{#gamma}";
    std::string label_PhotonPhi = "#phi^{#gamma}";
    std::string label_jetpt  = "jet p_{T} [GeV]";
    std::string label_jeteta = "jet #eta [GeV]";
    std::string label_jetphi = "jet #phi [GeV]";
    std::string label_jetE   = "jet E [GeV]";
    std::string label_j1pt = "j_{1} p_{T} [GeV]";
    std::string label_j2pt = "j_{2} p_{T} [GeV]";
    std::string label_j3pt = "j_{3} p_{T} [GeV]";
    std::string label_mll  = "m_{ll} [GeV]";
    std::string label_topPt = "top p_{T} [GeV]";
    std::string label_genTopPt = "gen top p_{T} [GeV]";
    std::string label_phopt = "p_{T}^{#gamma} [GeV]";
    std::string label_metg = "p_{T}^{#gamma (miss)} [GeV]";
    std::string label_ptcut_single = "GenEtaPt & GenPt";
    std::string label_ptcut_ratio  = "GenEtaPt / GenPt";
    std::string label_acc_single = "RecoEta & Reco";
    std::string label_acc_ratio  = "RecoEta / Reco";
    std::string label_matched_single = "RecoEtaPtMatched & RecoEtaPt";
    std::string label_matched_ratio  = "RecoEtaPtMatched / RecoEtaPt";
    std::string label_iso_single = "RecoIso & RecoEtaPtMatched";
    std::string label_iso_ratio  = "RecoIso / RecoEtaPtMatched";
    // make a map of labels
    // there is no Gen Iso
    //std::vector<std::pair<std::string,std::string>> cutlevels_electrons
    std::map<std::string, std::string> label_map = {
        {"GenAcc_single",     "GenEta & Gen"},
        {"GenAcc_ratio",      "GenEta / Gen"},
        {"GenMatch_single",   "GenEtaPtMatched & GenEtaPt"},
        {"GenMatch_ratio",    "GenEtaPtMatched / GenEtaPt"},
        {"RecoAcc_single",    "RecoEta & Reco"},
        {"RecoAcc_ratio",     "RecoEta / Reco"},
        {"RecoIso_single",    "RecoIso & RecoEtaPt"},
        {"RecoIso_ratio",     "RecoIso / RecoEtaPt"},
        {"RecoMatch_single",  "RecoEtaPtMatched & RecoIso"},
        {"RecoMatch_ratio",   "RecoEtaPtMatched / RecoIso"},
    };

    //vector<Plotter::HistSummary> vh;
    vector<PHS> vh;
    
    ////////////////////////////////////////////////
    // --- version with weights for reference --- //
    ////////////////////////////////////////////////
    
    /*
    // Datasetsummaries we are using                                                                                                        
    // no weight (genWeight deals with negative weights); also add btag weights here                                                        
    PDS dsData_SingleMuon("Data",         fileMap["Data_SingleMuon"], "passMuonTrigger",   "");
    PDS dsDY_mu(          "DY #mu",       fileMap["DYJetsToLL"],      "",        "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor");
    PDS dsDYInc_mu(       "DY HT<100",    fileMap["IncDY"],           "genHT<100",   "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor");
    PDS dsDY_elec(        "DY e",         fileMap["DYJetsToLL"],      "",          "bTagSF_EventWeightSimple_Central;_PUweightFactor"); // do not use muTrigWgt for electrons (it is 0.0)
    PDS dsDYInc_elec(     "DY HT<100",    fileMap["IncDY"],           "genHT<100",   "bTagSF_EventWeightSimple_Central;_PUweightFactor"); // do not use muTrigWgt for electrons (it is 0.0)
    PDS dsPhoton(         "#gamma+ jets", fileMap["GJets"],         "",            "bTagSF_EventWeightSimple_Central;_PUweightFactor");
    PDS dstt2l(           "t#bar{t}",     fileMap["TTbarNoHad"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;isr_Unc_Cent;_PUweightFactor");
    PDS dstW(             "Single t",     fileMap["SingleTopZinv"],   "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    PDS dsttZ(            "t#bar{t}Z",    fileMap["TTZ"],             "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    PDS dsT1tttt_gluino1200_lsp800("T1tttt_gluino1200_lsp800",     fileMap["Signal_T1tttt_mGluino1200_mLSP800"], "",  "");
    PDS dsT1tttt_gluino1500_lsp100("T1tttt_gluino1500_lsp100",     fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "",  "");
    PDS dsT1tttt_gluino2000_lsp100("T1tttt_gluino2000_lsp100",     fileMap["Signal_T1tttt_mGluino2000_mLSP100"], "",  "");
    PDS dsVV(             "Diboson",      fileMap["Diboson"],        "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    PDS dsRare(           "Rare",         fileMap["Rare"],           "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    std::vector<std::vector<PDS>> stack_MC = {{dsDY_mu, dsDYInc_mu}, {dstt2l}, {dstW}, {dsRare, dsVV, dsttZ}};

    // Apply data/mc njet weight for DY and ttbar                                                                                                                                    
    PDS dswDY(             "DY",         fileMap["DYJetsToLL"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;njWGJets;_PUweightFactor");
    PDS dswDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "muTrigWgt;bTagSF_EventWeightSimple_Central;njWGJets;_PUweightFactor");
    PDS dswtt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;nJetWgtTTbar;isr_Unc_Cent;_PUweightFactor");
    PDS dswtW(             "Single t",   fileMap["SingleTopZinv"],   "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    PDS dswttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    PDS dswVV(             "Diboson",    fileMap["Diboson"],         "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor;genWeight");
    PDS dswRare(           "Rare",       fileMap["Rare"],            "",            "muTrigWgt;bTagSF_EventWeightSimple_Central;genWeight;_PUweightFactor");
    std::vector<std::vector<PDS>> stackw_MC = {{dswDY, dswDYInc}, {dswtt2l}, {dswtW}, {dswRare, dswVV, dswttZ}};

    PDS dswwDY(             "DY",         fileMap["DYJetsToLL"],      "",            "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    PDS dswwDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",   "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    std::vector<std::vector<PDS>> stackww_MC = {{dswwDY, dswwDYInc}, {dswtt2l}, {dswtW}, {dswttZ}, {dswVV}, {dswRare, dswVV, dswttZ}};
    */

    /////////////////////////////////////
    // --- version without weights --- //
    /////////////////////////////////////
    
    // Datasetsummaries we are using                                                                                                        
    // no weight (genWeight deals with negative weights); also add btag weights here                                                        
    PDS dsData_SingleMuon("Data",         fileMap[MuonDataset],  "passMuonTrigger",   "");
    PDS dsDY_mu(          "DY #mu",       fileMap["DYJetsToLL" + yearTag],      "",   "");
    PDS dsDYInc_mu(       "DY HT<100",    fileMap["IncDY" + yearTag],           "",   "");
    PDS dsDY_elec(        "DY e",         fileMap["DYJetsToLL" + yearTag],      "",   ""); 
    PDS dsDYInc_elec(     "DY HT<100",    fileMap["IncDY" + yearTag],           "",   ""); 
    PDS dsPhoton(         "#gamma+ jets", fileMap["GJets" + yearTag],           "",   "");
    PDS dstt2l(           "t#bar{t}",     fileMap["TTbarNoHad" + yearTag],      "",   "");
    PDS dstW(             "Single t",     fileMap["SingleTopZinv" + yearTag],   "",   "");
    PDS dsttZ(            "t#bar{t}Z",    fileMap["TTZ" + yearTag],             "",   "");
    PDS dsVV(             "Diboson",      fileMap["Diboson" + yearTag],        "",    "");
    PDS dsRare(           "Rare",         fileMap["Rare" + yearTag],           "",    "");
    PDS dsT1tttt_gluino1200_lsp800("T1tttt_gluino1200_lsp800",     fileMap["Signal_T1tttt_mGluino1200_mLSP800"], "",  "");
    PDS dsT1tttt_gluino1500_lsp100("T1tttt_gluino1500_lsp100",     fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "",  "");
    PDS dsT1tttt_gluino2000_lsp100("T1tttt_gluino2000_lsp100",     fileMap["Signal_T1tttt_mGluino2000_mLSP100"], "",  "");
    std::vector<std::vector<PDS>> stack_MC = {{dsDY_mu, dsDYInc_mu}, {dstt2l}, {dstW}, {dsRare, dsVV, dsttZ}};

    // Apply data/mc njet weight for DY and ttbar                                                                                                                                    
    PDS dswDY(             "DY",         fileMap["DYJetsToLL" + yearTag],      "",            "");
    PDS dswDYInc(          "DY HT<100",  fileMap["IncDY" + yearTag],           "",            "");
    PDS dswtt2l(           "t#bar{t}",   fileMap["TTbarNoHad" + yearTag],      "",            "");
    PDS dswtW(             "Single t",   fileMap["SingleTopZinv" + yearTag],   "",            "");
    PDS dswttZ(            "t#bar{t}Z",  fileMap["TTZ" + yearTag],             "",            "");
    PDS dswVV(             "Diboson",    fileMap["Diboson" + yearTag],         "",            "");
    PDS dswRare(           "Rare",       fileMap["Rare" + yearTag],            "",            "");
    std::vector<std::vector<PDS>> stackw_MC = {{dswDY, dswDYInc}, {dswtt2l}, {dswtW}, {dswRare, dswVV, dswttZ}};

    PDS dswwDY(             "DY",         fileMap["DYJetsToLL" + yearTag],      "",            "");
    PDS dswwDYInc(          "DY HT<100",  fileMap["IncDY" + yearTag],           "",            "");
    std::vector<std::vector<PDS>> stackww_MC = {{dswwDY, dswwDYInc}, {dswtt2l}, {dswtW}, {dswttZ}, {dswVV}, {dswRare, dswVV, dswttZ}};


    // use DYInc if true, otherwise use DYJetsToLL (HT binned DY)
    bool useDYInc = false;
    // apply ISRWeight to ttbar only... 
    std::string ISRWeight = ";ISRWeight";
    
    auto makeStackMC_DiLepton = [&](const std::string& cuts, const std::string& weights)
    {
        PDS dsDYInc(         "DY Inc",       fileMap["IncDY" + yearTag],           cuts,   weights);
        PDS dsDY(            "DY",           fileMap["DYJetsToLL" + yearTag],      cuts,   weights);
        PDS dsTTbar(         "t#bar{t}",     fileMap["TTbarNoHad" + yearTag],      cuts,   weights + ISRWeight);
        PDS dsSingleTopZinv( "Single t",     fileMap["SingleTopZinv" + yearTag],   cuts,   weights);
        PDS dsRare(          "Rare",         fileMap["Rare" + yearTag],            cuts,   weights);
        PDS dsDiboson(       "Diboson",      fileMap["Diboson" + yearTag],         cuts,   weights);
        PDS dsTTZ(           "t#bar{t}Z",    fileMap["TTZ" + yearTag],             cuts,   weights);
        if (useDYInc)
        {
            std::vector<std::vector<PDS>> StackMC = {{dsDYInc}, {dsTTbar}, {dsSingleTopZinv}, {dsRare, dsDiboson, dsTTZ}};
            return StackMC;
        }
        else
        {
            std::vector<std::vector<PDS>> StackMC = {{dsDY}, {dsTTbar}, {dsSingleTopZinv}, {dsRare, dsDiboson, dsTTZ}};
            return StackMC;
        }
    };
    
    auto makeStackMC_DiLepton_Normalization = [&](const std::string& cuts, const std::string& weights)
    {
        PDS dsDYInc(            "IncZToLL",         fileMap["IncDY" + yearTag],           cuts,   weights);
        PDS dsDY(               "ZToLL",            fileMap["DYJetsToLL" + yearTag],      cuts,   weights);
        PDS dsTTbar(            "NoZToLL",          fileMap["TTbarNoHad" + yearTag],      cuts,   weights + ISRWeight);
        PDS dsSingleTopZinv(    "Single t",         fileMap["SingleTopZinv" + yearTag],   cuts,   weights);
        PDS dsRareZ(            "RareZ",            fileMap["RareZ" + yearTag],           cuts,   weights);
        PDS dsRareNoZ(          "RareNoZ",          fileMap["RareNoZ" + yearTag],         cuts,   weights);
        PDS dsDibosonZToLL(     "DibosonZToLL",     fileMap["DibosonZToLL" + yearTag],    cuts,   weights);
        PDS dsDibosonNoZToLL(   "DibosonNoZToLL",   fileMap["DibosonNoZToLL" + yearTag],  cuts,   weights);
        PDS dsTTZ(              "t#bar{t}Z",        fileMap["TTZ" + yearTag],             cuts,   weights);
        if (useDYInc)
        {
            // "Z to LL" MC and "No Z to LL" MC
            std::vector<std::vector<PDS>> StackMC = {{dsDYInc, dsRareZ, dsDibosonZToLL}, {dsTTbar, dsSingleTopZinv, dsRareNoZ, dsDibosonNoZToLL, dsTTZ}};
            return StackMC;
        }
        else
        {
            // "Z to LL" MC and "No Z to LL" MC
            std::vector<std::vector<PDS>> StackMC = {{dsDY, dsRareZ, dsDibosonZToLL}, {dsTTbar, dsSingleTopZinv, dsRareNoZ, dsDibosonNoZToLL, dsTTZ}};
            return StackMC;
        }
    };
    
    // standard MC
    auto makeStackMC_Photon = [&](const std::string& cuts, const std::string& weights)
    {
        PDS dsGJets(      "#gamma+jets",            fileMap["GJets" + yearTag],         cuts,   weights);
        PDS dsQCD(        "QCD",                    fileMap["QCD_Photon" + yearTag],    cuts,   weights);
        PDS dsWJetsToLNu( "W(l#nu)+jets",           fileMap["WJetsToLNu" + yearTag],    cuts,   weights);
        PDS dsTTG(        "t#bar{t}#gamma+jets",    fileMap["TTG" + yearTag],           cuts,   weights);
        PDS dsTTbar(      "t#bar{t}",               fileMap["TTbar" + yearTag],         cuts,   weights + ISRWeight);
        PDS dstW(         "tW",                     fileMap["tW" + yearTag],            cuts,   weights);
        PDS dsRare(       "Rare",                   fileMap["Rare_Photon" + yearTag],   cuts,   weights);
        PDS dsDiboson(    "Diboson",                fileMap["Diboson" + yearTag],       cuts,   weights);
        PDS dsTTZ(        "t#bar{t}Z",              fileMap["TTZ" + yearTag],           cuts,   weights);
        //std::vector<std::vector<PDS>> stack_gammaMC = {{dsGJets},{dsQCD},{dsWJets},{dsTTG},{dstt2l},{dstW},{dsVV},{dsRare,dsttZ}}; // from MakePhotonPlots.C for reference
        std::vector<std::vector<PDS>> StackMC = {{dsGJets}, {dsQCD}, {dsWJetsToLNu}, {dsTTG}, {dsTTbar}, {dstW}, {dsRare, dsDiboson, dsTTZ}};
        return StackMC;
    };
    // using fake and fragmented photons
    //auto makeStackMC_Photon = [&](const std::string& cuts, const std::string& weights)
    //{
    //    PDS dsGJets(            "#gamma+jets",            fileMap["GJets" + yearTag],         cuts + ";passPhotonSelectionDirect",      weights);
    //    PDS dsQCDFragmented(    "QCD Fragmented",         fileMap["QCD_Photon" + yearTag],    cuts + ";passPhotonSelectionFragmented",  weights);
    //    PDS dsQCDFake(          "QCD Fake",               fileMap["QCD_Photon" + yearTag],    cuts + ";passPhotonSelectionFake",        weights);
    //    PDS dsTTG(              "t#bar{t}#gamma+jets",    fileMap["TTG" + yearTag],           cuts + ";passPhotonSelectionDirect",      weights);
    //    PDS dsWJetsToLNu(       "W(l#nu)+jets",           fileMap["WJetsToLNu" + yearTag],    cuts,   weights);
    //    PDS dsTTbar(            "t#bar{t}",               fileMap["TTbar" + yearTag],         cuts,   weights + ISRWeight);
    //    PDS dstW(               "tW",                     fileMap["tW" + yearTag],            cuts,   weights);
    //    PDS dsRare(             "Rare",                   fileMap["Rare_Photon" + yearTag],   cuts,   weights);
    //    PDS dsDiboson(          "Diboson",                fileMap["Diboson" + yearTag],       cuts,   weights);
    //    PDS dsTTZ(              "t#bar{t}Z",              fileMap["TTZ" + yearTag],           cuts,   weights);
    //    //std::vector<std::vector<PDS>> stack_gammaMC = {{dsGJets},{dsQCD},{dsWJets},{dsTTG},{dstt2l},{dstW},{dsVV},{dsRare,dsttZ}}; // from MakePhotonPlots.C for reference
    //    std::vector<std::vector<PDS>> StackMC = {{dsGJets}, {dsQCDFragmented}, {dsQCDFake}, {dsTTG}, {dsWJetsToLNu},  {dsTTbar}, {dstW}, {dsRare, dsDiboson, dsTTZ}};
    //    return StackMC;
    //};



    std::vector<Plotter::CutFlowSummary> cutFlowSummaries;
    
    // Z to NuNu
    std::vector<std::string> Cuts_MC_ZNuNu = {
                              "",
                              "Pass_LeptonVeto",
                              HEMVeto,
                              "SAT_Pass_JetID",
                              "SAT_Pass_EventFilter" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET",
                              "SAT_Pass_HT",
                              "SAT_Pass_NJets20",
                             };

    // Electron
    std::vector<std::string> Cuts_Data_Electron = {
                              "",
                              "passElectronTrigger",
                              "Pass_MuonVeto",
                              "passElecPt",
                              "passDiElecSel",
                              "passElecZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Electron_Loose
    std::vector<std::string> Cuts_Data_Electron_Loose = {
                              "",
                              "passElectronTrigger",
                              "Pass_MuonVeto",
                              "passElecPt",
                              "passDiElecSel",
                              "passElecZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_Loose_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Electron_Mid
    std::vector<std::string> Cuts_Data_Electron_Mid = {
                              "",
                              "passElectronTrigger",
                              "Pass_MuonVeto",
                              "passElecPt",
                              "passDiElecSel",
                              "passElecZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_Mid_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Muon
    std::vector<std::string> Cuts_Data_Muon = {
                              "",
                              "passMuonTrigger",
                              "Pass_ElecVeto",
                              "passMuPt",
                              "passDiMuSel",
                              "passMuZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Muon_Loose
    std::vector<std::string> Cuts_Data_Muon_Loose = {
                              "",
                              "passMuonTrigger",
                              "Pass_ElecVeto",
                              "passMuPt",
                              "passDiMuSel",
                              "passMuZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_Loose_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Muon_Mid
    std::vector<std::string> Cuts_Data_Muon_Mid = {
                              "",
                              "passMuonTrigger",
                              "Pass_ElecVeto",
                              "passMuPt",
                              "passDiMuSel",
                              "passMuZinvSel",
                              HEMVeto_drLeptonCleaned,
                              "SAT_Pass_JetID_drLeptonCleaned",
                              "SAT_Pass_EventFilter_drLeptonCleaned" + Flag_ecalBadCalibFilter,
                              "SAT_Pass_MET_Mid_drLeptonCleaned",
                              "SAT_Pass_HT_drLeptonCleaned",
                              "SAT_Pass_NJets20_drLeptonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drLeptonCleaned",
                             };
    // Photon
    // pre-baseline cuts for Data and MC
    // apply passPhotonTrigger and Flag_eeBadScFilter to Data but not to MC
    std::string preBaselineCutsPhotonData = "passPhotonTrigger;Pass_LeptonVeto;passPhotonSelection;SAT_Pass_JetID_drPhotonCleaned;SAT_Pass_EventFilter_drPhotonCleaned;Flag_eeBadScFilter;MET_pt<250" + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned;
    std::string preBaselineCutsPhotonMC   = "Pass_LeptonVeto;passPhotonSelection;SAT_Pass_JetID_drPhotonCleaned;SAT_Pass_EventFilter_drPhotonCleaned;MET_pt<250" + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned;
    //standard baseline
    std::vector<std::string> Cuts_Data_Photon = {
                              "",
                              "passPhotonTrigger",
                              "Pass_LeptonVeto",
                              "passPhotonSelection",
                              HEMVeto_drPhotonCleaned,
                              "SAT_Pass_JetID_drPhotonCleaned",
                              "SAT_Pass_EventFilter_drPhotonCleaned" + Flag_ecalBadCalibFilter,
                              "MET_pt<250",
                              "SAT_Pass_MET_drPhotonCleaned",
                              "SAT_Pass_HT_drPhotonCleaned",
                              "SAT_Pass_NJets20_drPhotonCleaned",
                              "SAT_Pass_dPhiMETLowDM_drPhotonCleaned",
                             };
    // all low dm and high dm cuts
    std::vector<std::string> Cuts_Data_Photon_LowDM_All = {
                              preBaselineCutsPhotonData,
                              "SAT_Pass_Baseline_drPhotonCleaned",
                              "nMergedTops_drPhotonCleaned=0",
                              "nResolvedTops_drPhotonCleaned=0",
                              "nWs_drPhotonCleaned=0",
                              "SAT_Pass_ISR_drPhotonCleaned",
                              "SAT_Pass_S_MET_drPhotonCleaned",
                              "SAT_Pass_MTB_LowDM_drPhotonCleaned",
                             };
    std::vector<std::string> Cuts_Data_Photon_HighDM_All = {
                              preBaselineCutsPhotonData,
                              "SAT_Pass_Baseline_drPhotonCleaned",
                              "SAT_Pass_dPhiMETHighDM_drPhotonCleaned",
                              "nBottoms_drPhotonCleaned>=1",
                              "nJets_drPhotonCleaned>=5",
                             };
    
    
    
    // Z to NuNu in validation bins
    std::vector<std::string> Cuts_MC_ZNuNu_LowDM          = Cuts_MC_ZNuNu;
    std::vector<std::string> Cuts_MC_ZNuNu_LowDM_HighMET  = Cuts_MC_ZNuNu;
    std::vector<std::string> Cuts_MC_ZNuNu_HighDM         = Cuts_MC_ZNuNu;
    Cuts_MC_ZNuNu_LowDM.push_back("SAT_Pass_dPhiMETLowDM");
    Cuts_MC_ZNuNu_LowDM.push_back("SAT_Pass_lowDM");
    Cuts_MC_ZNuNu_LowDM_HighMET.push_back("SAT_Pass_mid_dPhiMETLowDM");
    Cuts_MC_ZNuNu_LowDM_HighMET.push_back("SAT_Pass_lowDM_mid_dPhi");
    Cuts_MC_ZNuNu_HighDM.push_back("SAT_Pass_mid_dPhiMETHighDM");
    Cuts_MC_ZNuNu_HighDM.push_back("SAT_Pass_highDM_mid_dPhi");
    std::vector<std::string> CutLevels_MC_ZNuNu_LowDM         = SusyUtility::getCutLevels(Cuts_MC_ZNuNu_LowDM);
    std::vector<std::string> CutLevels_MC_ZNuNu_LowDM_HighMET = SusyUtility::getCutLevels(Cuts_MC_ZNuNu_LowDM_HighMET);
    std::vector<std::string> CutLevels_MC_ZNuNu_HighDM        = SusyUtility::getCutLevels(Cuts_MC_ZNuNu_HighDM);
    
    // Electron
    std::vector<std::string> Cuts_Data_Electron_LowDM  = Cuts_Data_Electron;
    std::vector<std::string> Cuts_Data_Electron_HighDM = Cuts_Data_Electron;
    Cuts_Data_Electron_LowDM.push_back("SAT_Pass_lowDM_drLeptonCleaned");
    Cuts_Data_Electron_HighDM.push_back("SAT_Pass_highDM_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Electron_LowDM   = Cuts_Data_Electron_LowDM;
    std::vector<std::string> Cuts_MC_Electron_HighDM  = Cuts_Data_Electron_HighDM;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Electron_LowDM[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Electron_HighDM[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data
    Cuts_MC_Electron_LowDM[1]  = "";
    Cuts_MC_Electron_HighDM[1] = "";
    std::vector<std::string> CutLevels_Data_Electron_LowDM  = SusyUtility::getCutLevels(Cuts_Data_Electron_LowDM);
    std::vector<std::string> CutLevels_Data_Electron_HighDM = SusyUtility::getCutLevels(Cuts_Data_Electron_HighDM);
    std::vector<std::string> CutLevels_MC_Electron_LowDM    = SusyUtility::getCutLevels(Cuts_MC_Electron_LowDM);
    std::vector<std::string> CutLevels_MC_Electron_HighDM   = SusyUtility::getCutLevels(Cuts_MC_Electron_HighDM);

    // Electron_Loose
    std::vector<std::string> Cuts_Data_Electron_LowDM_Loose  = Cuts_Data_Electron_Loose;
    std::vector<std::string> Cuts_Data_Electron_HighDM_Loose = Cuts_Data_Electron_Loose;
    Cuts_Data_Electron_LowDM_Loose.push_back("SAT_Pass_lowDM_Loose_drLeptonCleaned");
    Cuts_Data_Electron_HighDM_Loose.push_back("SAT_Pass_highDM_Loose_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Electron_LowDM_Loose   = Cuts_Data_Electron_LowDM_Loose;
    std::vector<std::string> Cuts_MC_Electron_HighDM_Loose  = Cuts_Data_Electron_HighDM_Loose;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Electron_LowDM_Loose[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Electron_HighDM_Loose[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data (Loose)
    Cuts_MC_Electron_LowDM_Loose[1]  = "";
    Cuts_MC_Electron_HighDM_Loose[1] = "";
    std::vector<std::string> CutLevels_Data_Electron_LowDM_Loose  = SusyUtility::getCutLevels(Cuts_Data_Electron_LowDM_Loose);
    std::vector<std::string> CutLevels_Data_Electron_HighDM_Loose = SusyUtility::getCutLevels(Cuts_Data_Electron_HighDM_Loose);
    std::vector<std::string> CutLevels_MC_Electron_LowDM_Loose    = SusyUtility::getCutLevels(Cuts_MC_Electron_LowDM_Loose);
    std::vector<std::string> CutLevels_MC_Electron_HighDM_Loose   = SusyUtility::getCutLevels(Cuts_MC_Electron_HighDM_Loose);

    // Electron_Mid
    std::vector<std::string> Cuts_Data_Electron_LowDM_Mid  = Cuts_Data_Electron_Mid;
    std::vector<std::string> Cuts_Data_Electron_HighDM_Mid = Cuts_Data_Electron_Mid;
    Cuts_Data_Electron_LowDM_Mid.push_back("SAT_Pass_lowDM_Mid_drLeptonCleaned");
    Cuts_Data_Electron_HighDM_Mid.push_back("SAT_Pass_highDM_Mid_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Electron_LowDM_Mid   = Cuts_Data_Electron_LowDM_Mid;
    std::vector<std::string> Cuts_MC_Electron_HighDM_Mid  = Cuts_Data_Electron_HighDM_Mid;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Electron_LowDM_Mid[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Electron_HighDM_Mid[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data (Mid)
    Cuts_MC_Electron_LowDM_Mid[1]  = "";
    Cuts_MC_Electron_HighDM_Mid[1] = "";
    std::vector<std::string> CutLevels_Data_Electron_LowDM_Mid  = SusyUtility::getCutLevels(Cuts_Data_Electron_LowDM_Mid);
    std::vector<std::string> CutLevels_Data_Electron_HighDM_Mid = SusyUtility::getCutLevels(Cuts_Data_Electron_HighDM_Mid);
    std::vector<std::string> CutLevels_MC_Electron_LowDM_Mid    = SusyUtility::getCutLevels(Cuts_MC_Electron_LowDM_Mid);
    std::vector<std::string> CutLevels_MC_Electron_HighDM_Mid   = SusyUtility::getCutLevels(Cuts_MC_Electron_HighDM_Mid);

    // Muon
    std::vector<std::string> Cuts_Data_Muon_LowDM  = Cuts_Data_Muon;
    std::vector<std::string> Cuts_Data_Muon_HighDM = Cuts_Data_Muon;
    Cuts_Data_Muon_LowDM.push_back("SAT_Pass_lowDM_drLeptonCleaned");
    Cuts_Data_Muon_HighDM.push_back("SAT_Pass_highDM_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Muon_LowDM   = Cuts_Data_Muon_LowDM;
    std::vector<std::string> Cuts_MC_Muon_HighDM  = Cuts_Data_Muon_HighDM;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Muon_LowDM[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Muon_HighDM[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data
    Cuts_MC_Muon_LowDM[1]  = "";
    Cuts_MC_Muon_HighDM[1] = "";
    std::vector<std::string> CutLevels_Data_Muon_LowDM  = SusyUtility::getCutLevels(Cuts_Data_Muon_LowDM);
    std::vector<std::string> CutLevels_Data_Muon_HighDM = SusyUtility::getCutLevels(Cuts_Data_Muon_HighDM);
    std::vector<std::string> CutLevels_MC_Muon_LowDM    = SusyUtility::getCutLevels(Cuts_MC_Muon_LowDM);
    std::vector<std::string> CutLevels_MC_Muon_HighDM   = SusyUtility::getCutLevels(Cuts_MC_Muon_HighDM);

    // Muon_Loose
    std::vector<std::string> Cuts_Data_Muon_LowDM_Loose  = Cuts_Data_Muon_Loose;
    std::vector<std::string> Cuts_Data_Muon_HighDM_Loose = Cuts_Data_Muon_Loose;
    Cuts_Data_Muon_LowDM_Loose.push_back("SAT_Pass_lowDM_Loose_drLeptonCleaned");
    Cuts_Data_Muon_HighDM_Loose.push_back("SAT_Pass_highDM_Loose_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Muon_LowDM_Loose   = Cuts_Data_Muon_LowDM_Loose;
    std::vector<std::string> Cuts_MC_Muon_HighDM_Loose  = Cuts_Data_Muon_HighDM_Loose;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Muon_LowDM_Loose[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Muon_HighDM_Loose[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data (Loose)
    Cuts_MC_Muon_LowDM_Loose[1]  = "";
    Cuts_MC_Muon_HighDM_Loose[1] = "";
    std::vector<std::string> CutLevels_Data_Muon_LowDM_Loose  = SusyUtility::getCutLevels(Cuts_Data_Muon_LowDM_Loose);
    std::vector<std::string> CutLevels_Data_Muon_HighDM_Loose = SusyUtility::getCutLevels(Cuts_Data_Muon_HighDM_Loose);
    std::vector<std::string> CutLevels_MC_Muon_LowDM_Loose    = SusyUtility::getCutLevels(Cuts_MC_Muon_LowDM_Loose);
    std::vector<std::string> CutLevels_MC_Muon_HighDM_Loose   = SusyUtility::getCutLevels(Cuts_MC_Muon_HighDM_Loose);

    // Muon_Mid
    std::vector<std::string> Cuts_Data_Muon_LowDM_Mid  = Cuts_Data_Muon_Mid;
    std::vector<std::string> Cuts_Data_Muon_HighDM_Mid = Cuts_Data_Muon_Mid;
    Cuts_Data_Muon_LowDM_Mid.push_back("SAT_Pass_lowDM_Mid_drLeptonCleaned");
    Cuts_Data_Muon_HighDM_Mid.push_back("SAT_Pass_highDM_Mid_drLeptonCleaned");
    std::vector<std::string> Cuts_MC_Muon_LowDM_Mid   = Cuts_Data_Muon_LowDM_Mid;
    std::vector<std::string> Cuts_MC_Muon_HighDM_Mid  = Cuts_Data_Muon_HighDM_Mid;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Muon_LowDM_Mid[8]  += ";Flag_eeBadScFilter";
    Cuts_Data_Muon_HighDM_Mid[8] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data (Mid)
    Cuts_MC_Muon_LowDM_Mid[1]  = "";
    Cuts_MC_Muon_HighDM_Mid[1] = "";
    std::vector<std::string> CutLevels_Data_Muon_LowDM_Mid  = SusyUtility::getCutLevels(Cuts_Data_Muon_LowDM_Mid);
    std::vector<std::string> CutLevels_Data_Muon_HighDM_Mid = SusyUtility::getCutLevels(Cuts_Data_Muon_HighDM_Mid);
    std::vector<std::string> CutLevels_MC_Muon_LowDM_Mid    = SusyUtility::getCutLevels(Cuts_MC_Muon_LowDM_Mid);
    std::vector<std::string> CutLevels_MC_Muon_HighDM_Mid   = SusyUtility::getCutLevels(Cuts_MC_Muon_HighDM_Mid);

    // Photon
    std::vector<std::string> Cuts_Data_Photon_LowDM  = Cuts_Data_Photon;
    std::vector<std::string> Cuts_Data_Photon_HighDM = Cuts_Data_Photon;
    Cuts_Data_Photon_LowDM.push_back("SAT_Pass_lowDM_drPhotonCleaned");
    Cuts_Data_Photon_HighDM.push_back("SAT_Pass_highDM_drPhotonCleaned");
    std::vector<std::string> Cuts_MC_Photon_LowDM   = Cuts_Data_Photon_LowDM;
    std::vector<std::string> Cuts_MC_Photon_HighDM  = Cuts_Data_Photon_HighDM;
    // Note: Apply Flag_eeBadScFilter to Data but not MC
    Cuts_Data_Photon_LowDM[6]  += ";Flag_eeBadScFilter";
    Cuts_Data_Photon_HighDM[6] += ";Flag_eeBadScFilter";
    // no trigger cut for MC; placeholder so that bins match data
    Cuts_MC_Photon_LowDM[1]  = "";
    Cuts_MC_Photon_HighDM[1] = "";
    std::vector<std::string> CutLevels_Data_Photon_LowDM  = SusyUtility::getCutLevels(Cuts_Data_Photon_LowDM);
    std::vector<std::string> CutLevels_Data_Photon_HighDM = SusyUtility::getCutLevels(Cuts_Data_Photon_HighDM);
    std::vector<std::string> CutLevels_MC_Photon_LowDM    = SusyUtility::getCutLevels(Cuts_MC_Photon_LowDM);
    std::vector<std::string> CutLevels_MC_Photon_HighDM   = SusyUtility::getCutLevels(Cuts_MC_Photon_HighDM);
    // All low dm and high dm cuts
    std::vector<std::string> Cuts_MC_Photon_LowDM_All     = Cuts_Data_Photon_LowDM_All;
    std::vector<std::string> Cuts_MC_Photon_HighDM_All    = Cuts_Data_Photon_HighDM_All;
    Cuts_MC_Photon_LowDM_All[0]  = preBaselineCutsPhotonMC;
    Cuts_MC_Photon_HighDM_All[0] = preBaselineCutsPhotonMC;
    std::vector<std::string> CutLevels_Data_Photon_LowDM_All   = SusyUtility::getCutLevels(Cuts_Data_Photon_LowDM_All);
    std::vector<std::string> CutLevels_Data_Photon_HighDM_All  = SusyUtility::getCutLevels(Cuts_Data_Photon_HighDM_All);
    std::vector<std::string> CutLevels_MC_Photon_LowDM_All     = SusyUtility::getCutLevels(Cuts_MC_Photon_LowDM_All);
    std::vector<std::string> CutLevels_MC_Photon_HighDM_All    = SusyUtility::getCutLevels(Cuts_MC_Photon_HighDM_All);
            
    // MC DY comparison of electron and muon selection
            
    // use IncDY or DYJetsToLL
    std::string DY_sample = "";
    if (useDYInc) DY_sample = "IncDY";
    else          DY_sample = "DYJetsToLL";

    // only di-lepton selection
    PDS dsDY_LowDM("DY",                    fileMap[DY_sample + yearTag],  "",               "");
    PDS dsDY_LowDM_Electron("DY Electron",  fileMap[DY_sample + yearTag],  "passDiElecSel",  "DiElecSF");
    PDS dsDY_LowDM_Muon("DY Muon",          fileMap[DY_sample + yearTag],  "passDiMuSel",    "DiMuSF");
    PDS dsDY_HighDM("DY",                   fileMap[DY_sample + yearTag],  "",               "");
    PDS dsDY_HighDM_Electron("DY Electron", fileMap[DY_sample + yearTag],  "passDiElecSel",  "DiElecSF");
    PDS dsDY_HighDM_Muon("DY Muon",         fileMap[DY_sample + yearTag],  "passDiMuSel",    "DiMuSF");

    // number of Electrons
    PDC dcMC_LowDM_nElectrons("single",           "nElectrons_Stop0l", {dsDY_LowDM});
    PDC dcMC_LowDM_Electron_nElectrons("single",  "nElectrons_Stop0l", {dsDY_LowDM_Electron});
    PDC dcMC_LowDM_Muon_nElectrons("single",      "nElectrons_Stop0l", {dsDY_LowDM_Muon});
    PDC dcMC_HighDM_nElectrons("single",          "nElectrons_Stop0l", {dsDY_HighDM});
    PDC dcMC_HighDM_Electron_nElectrons("single", "nElectrons_Stop0l", {dsDY_HighDM_Electron});
    PDC dcMC_HighDM_Muon_nElectrons("single",     "nElectrons_Stop0l", {dsDY_HighDM_Muon});
    // number of Muons
    PDC dcMC_LowDM_nMuons("single",           "nMuons_Stop0l", {dsDY_LowDM});
    PDC dcMC_LowDM_Electron_nMuons("single",  "nMuons_Stop0l", {dsDY_LowDM_Electron});
    PDC dcMC_LowDM_Muon_nMuons("single",      "nMuons_Stop0l", {dsDY_LowDM_Muon});
    PDC dcMC_HighDM_nMuons("single",          "nMuons_Stop0l", {dsDY_HighDM});
    PDC dcMC_HighDM_Electron_nMuons("single", "nMuons_Stop0l", {dsDY_HighDM_Electron});
    PDC dcMC_HighDM_Muon_nMuons("single",     "nMuons_Stop0l", {dsDY_HighDM_Muon});
    // bestRecoZM
    PDC dcMC_LowDM_bestRecoZM("single",           "bestRecoZM", {dsDY_LowDM});
    PDC dcMC_LowDM_Electron_bestRecoZM("single",  "bestRecoZM", {dsDY_LowDM_Electron});
    PDC dcMC_LowDM_Muon_bestRecoZM("single",      "bestRecoZM", {dsDY_LowDM_Muon});
    PDC dcMC_HighDM_bestRecoZM("single",          "bestRecoZM", {dsDY_HighDM});
    PDC dcMC_HighDM_Electron_bestRecoZM("single", "bestRecoZM", {dsDY_HighDM_Electron});
    PDC dcMC_HighDM_Muon_bestRecoZM("single",     "bestRecoZM", {dsDY_HighDM_Muon});
    
    vh.push_back(PHS("MC_LowDM_nElectrons" + eraTag,  {dcMC_LowDM_nElectrons, dcMC_LowDM_Electron_nElectrons, dcMC_LowDM_Muon_nElectrons},    {2, 3}, "", 10,  0,  10, true, false, "nElectrons", "Events"));
    vh.push_back(PHS("MC_HighDM_nElectrons" + eraTag, {dcMC_HighDM_nElectrons, dcMC_HighDM_Electron_nElectrons, dcMC_HighDM_Muon_nElectrons}, {2, 3}, "", 10,  0,  10, true, false, "nElectrons", "Events"));
    vh.push_back(PHS("MC_LowDM_nMuons" + eraTag,      {dcMC_LowDM_nMuons, dcMC_LowDM_Electron_nMuons, dcMC_LowDM_Muon_nMuons},                {2, 3}, "", 10,  0,  10, true, false, "nMuons", "Events"));
    vh.push_back(PHS("MC_HighDM_nMuons" + eraTag,     {dcMC_HighDM_nMuons, dcMC_HighDM_Electron_nMuons, dcMC_HighDM_Muon_nMuons},             {2, 3}, "", 10,  0,  10, true, false, "nMuons", "Events"));
    vh.push_back(PHS("MC_LowDM_bestRecoZM" + eraTag,  {dcMC_LowDM_bestRecoZM, dcMC_LowDM_Electron_bestRecoZM, dcMC_LowDM_Muon_bestRecoZM},    {2, 3}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
    vh.push_back(PHS("MC_HighDM_bestRecoZM" + eraTag, {dcMC_HighDM_bestRecoZM, dcMC_HighDM_Electron_bestRecoZM, dcMC_HighDM_Muon_bestRecoZM}, {2, 3}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
    
    // Jet pt cuts
    std::vector<std::string> JetPtCuts = {"_jetpt20", "_jetpt30", "_jetpt40"};
    // Photon ID Selections
    std::vector<std::string> PhotonIDSelections = {"passPhotonSelectionLoose", "passPhotonSelectionMedium", "passPhotonSelectionTight"};
    
    // Validation Bins Selection
    std::map<std::string, std::string> validation_cuts_low_dm = {
        {"NBeq0_NSVeq0",     ";nBottoms_drLeptonCleaned=0;nSoftBottoms_drLeptonCleaned=0"},
        {"NBeq0_NSVge1",     ";nBottoms_drLeptonCleaned=0;nSoftBottoms_drLeptonCleaned>=1"},
        {"NBeq1_NSVeq0",     ";nBottoms_drLeptonCleaned=1;nSoftBottoms_drLeptonCleaned=0"},
        {"NBeq1_NSVge1",     ";nBottoms_drLeptonCleaned=1;nSoftBottoms_drLeptonCleaned>=1"},
        {"NBge1_NSVeq0",     ";nBottoms_drLeptonCleaned>=1;nSoftBottoms_drLeptonCleaned=0"},
        {"NBge1_NSVge1",     ";nBottoms_drLeptonCleaned>=1;nSoftBottoms_drLeptonCleaned>=1"},
        {"NBge2",            ";nBottoms_drLeptonCleaned>=2"},
    };
    std::map<std::string, std::string> validation_cuts_high = {
        {"NBeq1",     ";nBottoms_drLeptonCleaned=1"},
        {"NBge2",     ";nBottoms_drLeptonCleaned>=2"},
    };


    // begin loop over jet pt cuts 
    for (const auto& JetPtCut : JetPtCuts) 
    {
        //suffix
        std::string suffix = JetPtCut + eraTag;
        // baseline selections for lepton control regions
        SAT_Pass_lowDM        = ";SAT_Pass_lowDM_drLeptonCleaned"  + JetPtCut;
        SAT_Pass_highDM       = ";SAT_Pass_highDM_drLeptonCleaned" + JetPtCut;
        SAT_Pass_lowDM_Loose  = ";SAT_Pass_lowDM_Loose_drLeptonCleaned"  + JetPtCut;
        SAT_Pass_highDM_Loose = ";SAT_Pass_highDM_Loose_drLeptonCleaned" + JetPtCut;
        SAT_Pass_lowDM_Mid    = ";SAT_Pass_lowDM_Mid_drLeptonCleaned"  + JetPtCut;
        SAT_Pass_highDM_Mid   = ";SAT_Pass_highDM_Mid_drLeptonCleaned" + JetPtCut;
        if (doHEMVeto)
        {
            semicolon_HEMVeto_drLeptonCleaned = ";SAT_Pass_HEMVeto_drLeptonCleaned" + JetPtCut;
        }
        // di-electron
        if (doDataMCElectron)
        {
            // use DiElecTriggerEffPt instead of Stop0l_trigger_eff_Zee_pt because it applies a Z mass cut
            std::string ElectronWeights = "DiElecTriggerEffPt;DiElecSF;BTagWeight" + PrefireWeight + puWeight;
            // bestRecoZM used to calculate normalization
            // Validation Bins Selection
            for (const auto& cut : validation_cuts_low_dm)
            {
                PDS dsData_Electron_LowDM_noZMassCut("Data",  fileMap[ElectronDataset],        "Flag_eeBadScFilter;passElecZinvSel;passElectronTrigger" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, "");
                std::vector<std::vector<PDS>> StackMC_Electron_LowDM_noZMassCut                = makeStackMC_DiLepton(                "passElecZinvSel" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, ElectronWeights);
                std::vector<std::vector<PDS>> StackMC_Electron_LowDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(  "passElecZinvSel" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, ElectronWeights);
                // bestRecoZM
                PDC dcData_Electron_LowDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Electron_LowDM_noZMassCut});
                PDC dcMC_Electron_LowDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_LowDM_noZMassCut);
                PDC dcMC_Electron_LowDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_LowDM_Normalization_noZMassCut);
                vh.push_back(PHS("DataMC_Electron_LowDM_bestRecoZM_50to250_" + cut.first + suffix,                  {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_LowDM_bestRecoZM_0to400_" + cut.first + suffix,                   {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_LowDM_Normalization_bestRecoZM_50to250_" + cut.first + suffix,    {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400_" + cut.first + suffix,     {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            }
            for (const auto& cut : validation_cuts_high)
            {
                PDS dsData_Electron_HighDM_noZMassCut("Data",  fileMap[ElectronDataset],        "Flag_eeBadScFilter;passElecZinvSel;passElectronTrigger" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, "");
                std::vector<std::vector<PDS>> StackMC_Electron_HighDM_noZMassCut                = makeStackMC_DiLepton(                "passElecZinvSel" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, ElectronWeights);
                std::vector<std::vector<PDS>> StackMC_Electron_HighDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(  "passElecZinvSel" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, ElectronWeights);
                // bestRecoZM
                PDC dcData_Electron_HighDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Electron_HighDM_noZMassCut});
                PDC dcMC_Electron_HighDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_HighDM_noZMassCut);
                PDC dcMC_Electron_HighDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_HighDM_Normalization_noZMassCut);
                vh.push_back(PHS("DataMC_Electron_HighDM_bestRecoZM_50to250_" + cut.first + suffix,                  {dcData_Electron_HighDM_bestRecoZM,   dcMC_Electron_HighDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_HighDM_bestRecoZM_0to400_" + cut.first + suffix,                   {dcData_Electron_HighDM_bestRecoZM,   dcMC_Electron_HighDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_HighDM_Normalization_bestRecoZM_50to250_" + cut.first + suffix,    {dcData_Electron_HighDM_bestRecoZM,   dcMC_Electron_HighDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Electron_HighDM_Normalization_bestRecoZM_0to400_" + cut.first + suffix,     {dcData_Electron_HighDM_bestRecoZM,   dcMC_Electron_HighDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            }

            // Data
            // Note: Apply Flag_eeBadScFilter to Data but not MC
            // without Z mass cut
            PDS dsData_Electron_LowDM_noZMassCut("Data",  fileMap[ElectronDataset],     "Flag_eeBadScFilter;passElecZinvSel;passElectronTrigger" + SAT_Pass_lowDM  + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_HighDM_noZMassCut("Data", fileMap[ElectronDataset],     "Flag_eeBadScFilter;passElecZinvSel;passElectronTrigger" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            // with on Z mass peak cut
            PDS dsData_Electron_LowDM("Data",  fileMap[ElectronDataset],                "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_HighDM("Data", fileMap[ElectronDataset],                "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_LowDM_Loose("Data",  fileMap[ElectronDataset],          "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_lowDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_HighDM_Loose("Data", fileMap[ElectronDataset],          "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_highDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_LowDM_Mid("Data",  fileMap[ElectronDataset],            "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_lowDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Electron_HighDM_Mid("Data", fileMap[ElectronDataset],            "Flag_eeBadScFilter;passElecZinvSelOnZMassPeak;passElectronTrigger"  + SAT_Pass_highDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            
            // MC
            // without Z mass cut
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM_noZMassCut                = makeStackMC_DiLepton(                "passElecZinvSel"  + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM_noZMassCut               = makeStackMC_DiLepton(                "passElecZinvSel"  + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(  "passElecZinvSel"  + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM_Normalization_noZMassCut = makeStackMC_DiLepton_Normalization(  "passElecZinvSel"  + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            // with on Z mass peak cut
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM                = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM               = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM_Loose          = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_lowDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM_Loose         = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_highDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM_Mid            = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_lowDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM_Mid           = makeStackMC_DiLepton(                 "passElecZinvSelOnZMassPeak" + SAT_Pass_highDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_LowDM_Normalization  = makeStackMC_DiLepton_Normalization(   "passElecZinvSelOnZMassPeak" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            std::vector<std::vector<PDS>> StackMC_Electron_HighDM_Normalization = makeStackMC_DiLepton_Normalization(   "passElecZinvSelOnZMassPeak" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, ElectronWeights);
            
            // n_jets
            PDC dcData_Electron_LowDM_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_nj( "data",   "nJets_drLeptonCleaned", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Electron_HighDM);
            
            // n_jets_Loose
            PDC dcData_Electron_LowDM_Loose_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_nj( "data",   "nJets_drLeptonCleaned", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Electron_HighDM_Loose);
            
            // n_jets_Mid
            PDC dcData_Electron_LowDM_Mid_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_nj( "data",   "nJets_drLeptonCleaned", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Electron_HighDM_Mid);
            
            // HT
            PDC dcData_Electron_LowDM_ht(  "data",   "HT_drLeptonCleaned", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ht( "data",   "HT_drLeptonCleaned", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Electron_HighDM);

            // HT_Loose
            PDC dcData_Electron_LowDM_Loose_ht(  "data",   "HT_drLeptonCleaned", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_ht( "data",   "HT_drLeptonCleaned", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Electron_HighDM_Loose);

            // HT_Mid
            PDC dcData_Electron_LowDM_Mid_ht(  "data",   "HT_drLeptonCleaned", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_ht( "data",   "HT_drLeptonCleaned", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Electron_HighDM_Mid);

            // met
            PDC dcData_Electron_LowDM_met(  "data",   "metWithLL", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_met( "data",   "metWithLL", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_met(    "stack",  "metWithLL", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_met(   "stack",  "metWithLL", StackMC_Electron_HighDM);
            PDC dcMC_Electron_LowDM_Normalization_met(    "stack",  "metWithLL", StackMC_Electron_LowDM_Normalization);
            PDC dcMC_Electron_HighDM_Normalization_met(   "stack",  "metWithLL", StackMC_Electron_HighDM_Normalization);
            
            // met_Loose
            PDC dcData_Electron_LowDM_Loose_met(  "data",   "metWithLL", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_met( "data",   "metWithLL", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_met(    "stack",  "metWithLL", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_met(   "stack",  "metWithLL", StackMC_Electron_HighDM_Loose);
            
            // met_Mid
            PDC dcData_Electron_LowDM_Mid_met(  "data",   "metWithLL", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_met( "data",   "metWithLL", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_met(    "stack",  "metWithLL", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_met(   "stack",  "metWithLL", StackMC_Electron_HighDM_Mid);
            
            // metphi
            PDC dcData_Electron_LowDM_metphi(  "data",   "metphiWithLL", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_metphi( "data",   "metphiWithLL", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_metphi(    "stack",  "metphiWithLL", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_metphi(   "stack",  "metphiWithLL", StackMC_Electron_HighDM);
            
            // jetpt
            PDC dcData_Electron_LowDM_jetpt(  "data",   "JetTLV(pt)", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_jetpt( "data",   "JetTLV(pt)", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_jetpt(    "stack",  "JetTLV(pt)", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_jetpt(   "stack",  "JetTLV(pt)", StackMC_Electron_HighDM);

            // jetpt_Loose
            PDC dcData_Electron_LowDM_Loose_jetpt(  "data",   "JetTLV(pt)", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_jetpt( "data",   "JetTLV(pt)", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_jetpt(    "stack",  "JetTLV(pt)", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_jetpt(   "stack",  "JetTLV(pt)", StackMC_Electron_HighDM_Loose);

            // jetpt_Mid
            PDC dcData_Electron_LowDM_Mid_jetpt(  "data",   "JetTLV(pt)", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_jetpt( "data",   "JetTLV(pt)", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_jetpt(    "stack",  "JetTLV(pt)", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_jetpt(   "stack",  "JetTLV(pt)", StackMC_Electron_HighDM_Mid);

            // jetpt_drLeptonCleaned
            PDC dcData_Electron_LowDM_jetpt_drLeptonCleaned(  "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_jetpt_drLeptonCleaned( "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_jetpt_drLeptonCleaned(    "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_jetpt_drLeptonCleaned(   "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_HighDM);

            // jetpt_drLeptonCleaned_Loose
            PDC dcData_Electron_LowDM_Loose_jetpt_drLeptonCleaned(  "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_jetpt_drLeptonCleaned( "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_jetpt_drLeptonCleaned(    "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_jetpt_drLeptonCleaned(   "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_HighDM_Loose);

            // jetpt_drLeptonCleaned_Mid
            PDC dcData_Electron_LowDM_Mid_jetpt_drLeptonCleaned(  "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_jetpt_drLeptonCleaned( "data",   "JetTLV_drLeptonCleaned(pt)", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_jetpt_drLeptonCleaned(    "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_jetpt_drLeptonCleaned(   "stack",  "JetTLV_drLeptonCleaned(pt)", StackMC_Electron_HighDM_Mid);

            // jeteta
            PDC dcData_Electron_LowDM_jeteta(  "data",   "JetTLV(eta)", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_jeteta( "data",   "JetTLV(eta)", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_jeteta(    "stack",  "JetTLV(eta)", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_jeteta(   "stack",  "JetTLV(eta)", StackMC_Electron_HighDM);

            // jeteta_drLeptonCleaned
            PDC dcData_Electron_LowDM_jeteta_drLeptonCleaned(  "data",   "JetTLV_drLeptonCleaned(eta)", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_jeteta_drLeptonCleaned( "data",   "JetTLV_drLeptonCleaned(eta)", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_jeteta_drLeptonCleaned(    "stack",  "JetTLV_drLeptonCleaned(eta)", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_jeteta_drLeptonCleaned(   "stack",  "JetTLV_drLeptonCleaned(eta)", StackMC_Electron_HighDM);

            // electron pt
            PDC dcData_Electron_LowDM_ElecPt1(  "data",   "cutElecPt1", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ElecPt1( "data",   "cutElecPt1", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ElecPt1(    "stack",  "cutElecPt1", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ElecPt1(   "stack",  "cutElecPt1", StackMC_Electron_HighDM);
            PDC dcData_Electron_LowDM_ElecPt2(  "data",   "cutElecPt2", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ElecPt2( "data",   "cutElecPt2", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ElecPt2(    "stack",  "cutElecPt2", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ElecPt2(   "stack",  "cutElecPt2", StackMC_Electron_HighDM);

            // electron eta
            PDC dcData_Electron_LowDM_ElecEta1(  "data",   "cutElecEta1", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ElecEta1( "data",   "cutElecEta1", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ElecEta1(    "stack",  "cutElecEta1", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ElecEta1(   "stack",  "cutElecEta1", StackMC_Electron_HighDM);
            PDC dcData_Electron_LowDM_ElecEta2(  "data",   "cutElecEta2", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ElecEta2( "data",   "cutElecEta2", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ElecEta2(    "stack",  "cutElecEta2", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ElecEta2(   "stack",  "cutElecEta2", StackMC_Electron_HighDM);

            // bestRecoZPt
            PDC dcData_Electron_LowDM_bestRecoZPt(  "data",   "bestRecoZPt", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_bestRecoZPt( "data",   "bestRecoZPt", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_bestRecoZPt(    "stack",  "bestRecoZPt", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_bestRecoZPt(   "stack",  "bestRecoZPt", StackMC_Electron_HighDM);
            
            // bestRecoZM
            PDC dcData_Electron_LowDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Electron_LowDM_noZMassCut});
            PDC dcData_Electron_HighDM_bestRecoZM( "data",   "bestRecoZM", {dsData_Electron_HighDM_noZMassCut});
            PDC dcMC_Electron_LowDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_LowDM_noZMassCut);
            PDC dcMC_Electron_HighDM_bestRecoZM(   "stack",  "bestRecoZM", StackMC_Electron_HighDM_noZMassCut);
            PDC dcMC_Electron_LowDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Electron_LowDM_Normalization_noZMassCut);
            PDC dcMC_Electron_HighDM_Normalization_bestRecoZM(   "stack",  "bestRecoZM", StackMC_Electron_HighDM_Normalization_noZMassCut);
            
            // mtb
            PDC dcData_Electron_LowDM_mtb(  "data",   "mtb", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_mtb( "data",   "mtb", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_mtb(    "stack",  "mtb", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_mtb(   "stack",  "mtb", StackMC_Electron_HighDM);
            
            // mtb_drLeptonCleaned
            PDC dcData_Electron_LowDM_mtb_drLeptonCleaned(  "data",   "mtb_drLeptonCleaned", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_mtb_drLeptonCleaned( "data",   "mtb_drLeptonCleaned", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_mtb_drLeptonCleaned(    "stack",  "mtb_drLeptonCleaned", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_mtb_drLeptonCleaned(   "stack",  "mtb_drLeptonCleaned", StackMC_Electron_HighDM);
            
            // mtb_drLeptonCleaned_Loose
            PDC dcData_Electron_LowDM_Loose_mtb_drLeptonCleaned(  "data",   "mtb_drLeptonCleaned", {dsData_Electron_LowDM_Loose});
            PDC dcData_Electron_HighDM_Loose_mtb_drLeptonCleaned( "data",   "mtb_drLeptonCleaned", {dsData_Electron_HighDM_Loose});
            PDC dcMC_Electron_LowDM_Loose_mtb_drLeptonCleaned(    "stack",  "mtb_drLeptonCleaned", StackMC_Electron_LowDM_Loose);
            PDC dcMC_Electron_HighDM_Loose_mtb_drLeptonCleaned(   "stack",  "mtb_drLeptonCleaned", StackMC_Electron_HighDM_Loose);
            
            // mtb_drLeptonCleaned_Mid
            PDC dcData_Electron_LowDM_Mid_mtb_drLeptonCleaned(  "data",   "mtb_drLeptonCleaned", {dsData_Electron_LowDM_Mid});
            PDC dcData_Electron_HighDM_Mid_mtb_drLeptonCleaned( "data",   "mtb_drLeptonCleaned", {dsData_Electron_HighDM_Mid});
            PDC dcMC_Electron_LowDM_Mid_mtb_drLeptonCleaned(    "stack",  "mtb_drLeptonCleaned", StackMC_Electron_LowDM_Mid);
            PDC dcMC_Electron_HighDM_Mid_mtb_drLeptonCleaned(   "stack",  "mtb_drLeptonCleaned", StackMC_Electron_HighDM_Mid);
            
            // Stop0l_Mtb
            PDC dcData_Electron_LowDM_Stop0l_Mtb(  "data",   "Stop0l_Mtb", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_Stop0l_Mtb( "data",   "Stop0l_Mtb", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_Stop0l_Mtb(    "stack",  "Stop0l_Mtb", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_Stop0l_Mtb(   "stack",  "Stop0l_Mtb", StackMC_Electron_HighDM);
            
            // ptb
            PDC dcData_Electron_LowDM_ptb(  "data",   "ptb", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ptb( "data",   "ptb", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ptb(    "stack",  "ptb", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ptb(   "stack",  "ptb", StackMC_Electron_HighDM);
            
            // ptb_drLeptonCleaned
            PDC dcData_Electron_LowDM_ptb_drLeptonCleaned(  "data",   "ptb_drLeptonCleaned", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ptb_drLeptonCleaned( "data",   "ptb_drLeptonCleaned", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ptb_drLeptonCleaned(    "stack",  "ptb_drLeptonCleaned", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ptb_drLeptonCleaned(   "stack",  "ptb_drLeptonCleaned", StackMC_Electron_HighDM);
            
            // Stop0l_Ptb
            PDC dcData_Electron_LowDM_Stop0l_Ptb(  "data",   "Stop0l_Ptb", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_Stop0l_Ptb( "data",   "Stop0l_Ptb", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_Stop0l_Ptb(    "stack",  "Stop0l_Ptb", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_Stop0l_Ptb(   "stack",  "Stop0l_Ptb", StackMC_Electron_HighDM);

            // ISRJetPt
            PDC dcData_Electron_LowDM_ISRJetPt(  "data",   "ISRJetPt", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ISRJetPt( "data",   "ISRJetPt", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ISRJetPt(    "stack",  "ISRJetPt", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ISRJetPt(   "stack",  "ISRJetPt", StackMC_Electron_HighDM);
            
            // ISRJetPt_drLeptonCleaned
            PDC dcData_Electron_LowDM_ISRJetPt_drLeptonCleaned(  "data",   "ISRJetPt_drLeptonCleaned", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_ISRJetPt_drLeptonCleaned( "data",   "ISRJetPt_drLeptonCleaned", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_ISRJetPt_drLeptonCleaned(    "stack",  "ISRJetPt_drLeptonCleaned", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_ISRJetPt_drLeptonCleaned(   "stack",  "ISRJetPt_drLeptonCleaned", StackMC_Electron_HighDM);
            
            // Stop0l_ISRJetPt
            PDC dcData_Electron_LowDM_Stop0l_ISRJetPt(  "data",   "Stop0l_ISRJetPt", {dsData_Electron_LowDM});
            PDC dcData_Electron_HighDM_Stop0l_ISRJetPt( "data",   "Stop0l_ISRJetPt", {dsData_Electron_HighDM});
            PDC dcMC_Electron_LowDM_Stop0l_ISRJetPt(    "stack",  "Stop0l_ISRJetPt", StackMC_Electron_LowDM);
            PDC dcMC_Electron_HighDM_Stop0l_ISRJetPt(   "stack",  "Stop0l_ISRJetPt", StackMC_Electron_HighDM);
            
            // dphi
            std::vector<PDC> dcVecData_Electron_LowDM_dPhi;
            std::vector<PDC> dcVecData_Electron_HighDM_dPhi;
            std::vector<PDC> dcVecMC_Electron_LowDM_dPhi;
            std::vector<PDC> dcVecMC_Electron_HighDM_dPhi;
            for (int i = 0; i < 4; i++)
            {
                std::string var = "dPhiVec_drLeptonCleaned[" + std::to_string(i) + "]";
                dcVecData_Electron_LowDM_dPhi.push_back(    PDC("data", var, {dsData_Electron_LowDM}));
                dcVecData_Electron_HighDM_dPhi.push_back(   PDC("data", var, {dsData_Electron_HighDM}));
                dcVecMC_Electron_LowDM_dPhi.push_back(      PDC("stack", var, StackMC_Electron_LowDM));
                dcVecMC_Electron_HighDM_dPhi.push_back(     PDC("stack", var, StackMC_Electron_HighDM));
            }
            
            // Standard selection
            vh.push_back(PHS("DataMC_Electron_LowDM_nj" + suffix,                                  {dcData_Electron_LowDM_nj,   dcMC_Electron_LowDM_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_nj" + suffix,                                 {dcData_Electron_HighDM_nj,  dcMC_Electron_HighDM_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ht" + suffix,                                  {dcData_Electron_LowDM_ht,   dcMC_Electron_LowDM_ht},   {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ht" + suffix,                                 {dcData_Electron_HighDM_ht,  dcMC_Electron_HighDM_ht},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_met" + suffix,                                 {dcData_Electron_LowDM_met,  dcMC_Electron_LowDM_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_met" + suffix,                                {dcData_Electron_HighDM_met, dcMC_Electron_HighDM_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_metphi" + suffix,                              {dcData_Electron_LowDM_metphi,  dcMC_Electron_LowDM_metphi},  {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_metphi" + suffix,                             {dcData_Electron_HighDM_metphi, dcMC_Electron_HighDM_metphi}, {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_jetpt" + suffix,                               {dcData_Electron_LowDM_jetpt,   dcMC_Electron_LowDM_jetpt},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_jetpt" + suffix,                              {dcData_Electron_HighDM_jetpt,  dcMC_Electron_HighDM_jetpt},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_jetpt_drLeptonCleaned" + suffix,               {dcData_Electron_LowDM_jetpt_drLeptonCleaned,   dcMC_Electron_LowDM_jetpt_drLeptonCleaned},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_jetpt_drLeptonCleaned" + suffix,              {dcData_Electron_HighDM_jetpt_drLeptonCleaned,  dcMC_Electron_HighDM_jetpt_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_jeteta" + suffix,                              {dcData_Electron_LowDM_jeteta,  dcMC_Electron_LowDM_jeteta},  {1, 2}, "",   nBins,  minEta, maxEta, true, false, label_jeteta, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_jeteta" + suffix,                             {dcData_Electron_HighDM_jeteta, dcMC_Electron_HighDM_jeteta}, {1, 2}, "",   nBins,  minEta, maxEta, true, false, label_jeteta, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_jeteta_drLeptonCleaned" + suffix,              {dcData_Electron_LowDM_jeteta_drLeptonCleaned,  dcMC_Electron_LowDM_jeteta_drLeptonCleaned},  {1, 2}, "",   nBins,  minEta, maxEta, true, false, "lepton cleaned "+label_jeteta, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_jeteta_drLeptonCleaned" + suffix,             {dcData_Electron_HighDM_jeteta_drLeptonCleaned, dcMC_Electron_HighDM_jeteta_drLeptonCleaned}, {1, 2}, "",   nBins,  minEta, maxEta, true, false, "lepton cleaned "+label_jeteta, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ElecPt1" + suffix,                             {dcData_Electron_LowDM_ElecPt1,   dcMC_Electron_LowDM_ElecPt1},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ElecPt1, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ElecPt1" + suffix,                            {dcData_Electron_HighDM_ElecPt1,  dcMC_Electron_HighDM_ElecPt1}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ElecPt1, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ElecPt2" + suffix,                             {dcData_Electron_LowDM_ElecPt2,   dcMC_Electron_LowDM_ElecPt2},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ElecPt2, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ElecPt2" + suffix,                            {dcData_Electron_HighDM_ElecPt2,  dcMC_Electron_HighDM_ElecPt2}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ElecPt2, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ElecEta1" + suffix,                            {dcData_Electron_LowDM_ElecEta1,  dcMC_Electron_LowDM_ElecEta1},  {1, 2}, "", nBins,  minEta, maxEta, true, false, label_ElecEta1, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ElecEta1" + suffix,                           {dcData_Electron_HighDM_ElecEta1, dcMC_Electron_HighDM_ElecEta1}, {1, 2}, "", nBins,  minEta, maxEta, true, false, label_ElecEta1, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ElecEta2" + suffix,                            {dcData_Electron_LowDM_ElecEta2,  dcMC_Electron_LowDM_ElecEta2},  {1, 2}, "", nBins,  minEta, maxEta, true, false, label_ElecEta2, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ElecEta2" + suffix,                           {dcData_Electron_HighDM_ElecEta2, dcMC_Electron_HighDM_ElecEta2}, {1, 2}, "", nBins,  minEta, maxEta, true, false, label_ElecEta2, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_bestRecoZPt" + suffix,                         {dcData_Electron_LowDM_bestRecoZPt,  dcMC_Electron_LowDM_bestRecoZPt},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_bestRecoZPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_bestRecoZPt" + suffix,                        {dcData_Electron_HighDM_bestRecoZPt, dcMC_Electron_HighDM_bestRecoZPt}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_bestRecoZPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_bestRecoZM_50to250" + suffix,                  {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_bestRecoZM_50to250" + suffix,                 {dcData_Electron_HighDM_bestRecoZM,  dcMC_Electron_HighDM_bestRecoZM}, {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_bestRecoZM_0to400" + suffix,                   {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_bestRecoZM_0to400" + suffix,                  {dcData_Electron_HighDM_bestRecoZM,  dcMC_Electron_HighDM_bestRecoZM}, {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Normalization_bestRecoZM_50to250" + suffix,    {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Normalization_bestRecoZM_50to250" + suffix,   {dcData_Electron_HighDM_bestRecoZM,  dcMC_Electron_HighDM_Normalization_bestRecoZM}, {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Normalization_bestRecoZM_0to400" + suffix,     {dcData_Electron_LowDM_bestRecoZM,   dcMC_Electron_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Normalization_bestRecoZM_0to400" + suffix,    {dcData_Electron_HighDM_bestRecoZM,  dcMC_Electron_HighDM_Normalization_bestRecoZM}, {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_mtb" + suffix,                                 {dcData_Electron_LowDM_mtb,   dcMC_Electron_LowDM_mtb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_mtb" + suffix,                                {dcData_Electron_HighDM_mtb,  dcMC_Electron_HighDM_mtb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_mtb_drLeptonCleaned" + suffix,                 {dcData_Electron_LowDM_mtb_drLeptonCleaned,   dcMC_Electron_LowDM_mtb_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_mtb_drLeptonCleaned" + suffix,                {dcData_Electron_HighDM_mtb_drLeptonCleaned,  dcMC_Electron_HighDM_mtb_drLeptonCleaned}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Stop0l_Mtb" + suffix,                          {dcData_Electron_LowDM_Stop0l_Mtb,   dcMC_Electron_LowDM_Stop0l_Mtb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Stop0l_Mtb" + suffix,                         {dcData_Electron_HighDM_Stop0l_Mtb,  dcMC_Electron_HighDM_Stop0l_Mtb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ptb" + suffix,                                 {dcData_Electron_LowDM_ptb,   dcMC_Electron_LowDM_ptb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ptb" + suffix,                                {dcData_Electron_HighDM_ptb,  dcMC_Electron_HighDM_ptb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ptb_drLeptonCleaned" + suffix,                 {dcData_Electron_LowDM_ptb_drLeptonCleaned,       dcMC_Electron_LowDM_ptb_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ptb_drLeptonCleaned" + suffix,                {dcData_Electron_HighDM_ptb_drLeptonCleaned,      dcMC_Electron_HighDM_ptb_drLeptonCleaned}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Stop0l_Ptb" + suffix,                          {dcData_Electron_LowDM_Stop0l_Ptb,                dcMC_Electron_LowDM_Stop0l_Ptb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Stop0l_Ptb" + suffix,                         {dcData_Electron_HighDM_Stop0l_Ptb,               dcMC_Electron_HighDM_Stop0l_Ptb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ISRJetPt" + suffix,                            {dcData_Electron_LowDM_ISRJetPt,                  dcMC_Electron_LowDM_ISRJetPt},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ISRJetPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ISRJetPt" + suffix,                           {dcData_Electron_HighDM_ISRJetPt,                 dcMC_Electron_HighDM_ISRJetPt}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ISRJetPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_ISRJetPt_drLeptonCleaned" + suffix,            {dcData_Electron_LowDM_ISRJetPt_drLeptonCleaned,  dcMC_Electron_LowDM_ISRJetPt_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_ISRJetPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_ISRJetPt_drLeptonCleaned" + suffix,           {dcData_Electron_HighDM_ISRJetPt_drLeptonCleaned, dcMC_Electron_HighDM_ISRJetPt_drLeptonCleaned}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_ISRJetPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Stop0l_ISRJetPt" + suffix,                     {dcData_Electron_LowDM_Stop0l_ISRJetPt,           dcMC_Electron_LowDM_Stop0l_ISRJetPt},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_ISRJetPt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Stop0l_ISRJetPt" + suffix,                    {dcData_Electron_HighDM_Stop0l_ISRJetPt,          dcMC_Electron_HighDM_Stop0l_ISRJetPt}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "Stop0l "+label_ISRJetPt, "Events"));
            
            // Loose selection
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_nj" + suffix,                            {dcData_Electron_LowDM_Loose_nj,   dcMC_Electron_LowDM_Loose_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_nj" + suffix,                           {dcData_Electron_HighDM_Loose_nj,  dcMC_Electron_HighDM_Loose_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_ht" + suffix,                            {dcData_Electron_LowDM_Loose_ht,   dcMC_Electron_LowDM_Loose_ht},   {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_ht" + suffix,                           {dcData_Electron_HighDM_Loose_ht,  dcMC_Electron_HighDM_Loose_ht},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_met" + suffix,                           {dcData_Electron_LowDM_Loose_met,  dcMC_Electron_LowDM_Loose_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_met" + suffix,                          {dcData_Electron_HighDM_Loose_met, dcMC_Electron_HighDM_Loose_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_jetpt" + suffix,                         {dcData_Electron_LowDM_Loose_jetpt,   dcMC_Electron_LowDM_Loose_jetpt},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_jetpt" + suffix,                        {dcData_Electron_HighDM_Loose_jetpt,  dcMC_Electron_HighDM_Loose_jetpt},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_jetpt_drLeptonCleaned" + suffix,         {dcData_Electron_LowDM_Loose_jetpt_drLeptonCleaned,   dcMC_Electron_LowDM_Loose_jetpt_drLeptonCleaned},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_jetpt_drLeptonCleaned" + suffix,        {dcData_Electron_HighDM_Loose_jetpt_drLeptonCleaned,  dcMC_Electron_HighDM_Loose_jetpt_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Loose_mtb_drLeptonCleaned" + suffix,           {dcData_Electron_LowDM_Loose_mtb_drLeptonCleaned,   dcMC_Electron_LowDM_Loose_mtb_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Loose_mtb_drLeptonCleaned" + suffix,          {dcData_Electron_HighDM_Loose_mtb_drLeptonCleaned,  dcMC_Electron_HighDM_Loose_mtb_drLeptonCleaned}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));

            // Mid selection
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_nj" + suffix,                              {dcData_Electron_LowDM_Mid_nj,   dcMC_Electron_LowDM_Mid_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_nj" + suffix,                             {dcData_Electron_HighDM_Mid_nj,  dcMC_Electron_HighDM_Mid_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_ht" + suffix,                              {dcData_Electron_LowDM_Mid_ht,   dcMC_Electron_LowDM_Mid_ht},   {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_ht" + suffix,                             {dcData_Electron_HighDM_Mid_ht,  dcMC_Electron_HighDM_Mid_ht},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_met" + suffix,                             {dcData_Electron_LowDM_Mid_met,  dcMC_Electron_LowDM_Mid_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_met" + suffix,                            {dcData_Electron_HighDM_Mid_met, dcMC_Electron_HighDM_Mid_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_jetpt" + suffix,                           {dcData_Electron_LowDM_Mid_jetpt,   dcMC_Electron_LowDM_Mid_jetpt},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_jetpt" + suffix,                          {dcData_Electron_HighDM_Mid_jetpt,  dcMC_Electron_HighDM_Mid_jetpt},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_jetpt_drLeptonCleaned" + suffix,           {dcData_Electron_LowDM_Mid_jetpt_drLeptonCleaned,   dcMC_Electron_LowDM_Mid_jetpt_drLeptonCleaned},   {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_jetpt_drLeptonCleaned" + suffix,          {dcData_Electron_HighDM_Mid_jetpt_drLeptonCleaned,  dcMC_Electron_HighDM_Mid_jetpt_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_jetpt, "Events"));
            vh.push_back(PHS("DataMC_Electron_LowDM_Mid_mtb_drLeptonCleaned" + suffix,             {dcData_Electron_LowDM_Mid_mtb_drLeptonCleaned,   dcMC_Electron_LowDM_Mid_mtb_drLeptonCleaned},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Electron_HighDM_Mid_mtb_drLeptonCleaned" + suffix,            {dcData_Electron_HighDM_Mid_mtb_drLeptonCleaned,  dcMC_Electron_HighDM_Mid_mtb_drLeptonCleaned}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, "lepton cleaned "+label_mtb, "Events"));

            // dphi
            for (int i = 0; i < 4; i++)
            {
                std::string nameLowDM  = "DataMC_Electron_LowDM_dPhi"  + std::to_string(i+1) + suffix;
                std::string nameHighDM = "DataMC_Electron_HighDM_dPhi" + std::to_string(i+1) + suffix;
                vh.push_back(PHS(nameLowDM,  {dcVecData_Electron_LowDM_dPhi[i], dcVecMC_Electron_LowDM_dPhi[i]},    {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
                vh.push_back(PHS(nameHighDM, {dcVecData_Electron_HighDM_dPhi[i], dcVecMC_Electron_HighDM_dPhi[i]},  {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
            }
        
            // cut flow: name, DataCollection, cutLevels
            // Plotter::CutFlowSummary::CutFlowSummary(std::string n, Plotter::DataCollection ns, std::vector<std::string> cutLevels)
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Electron_LowDM_met"  + suffix,   dcData_Electron_LowDM_met,  CutLevels_Data_Electron_LowDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Electron_HighDM_met" + suffix,   dcData_Electron_HighDM_met, CutLevels_Data_Electron_HighDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Electron_LowDM_met"    + suffix,   dcMC_Electron_LowDM_met,    CutLevels_MC_Electron_LowDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Electron_HighDM_met"   + suffix,   dcMC_Electron_HighDM_met,   CutLevels_MC_Electron_HighDM));
        }
        
        // di-muon
        if (doDataMCMuon)
        {
            // Use DiMuTriggerEffPt intead of Stop0l_trigger_eff_Zmumu_pt because it applies a Z mass cut
            std::string MuonWeights = "DiMuTriggerEffPt;DiMuSF;BTagWeight" + PrefireWeight + puWeight;
            // bestRecoZM used to calculate normalization
            // Validation Bins Selection
            for (const auto& cut : validation_cuts_low_dm)
            {
                PDS dsData_Muon_LowDM_noZMassCut("Data",  fileMap[MuonDataset],                "Flag_eeBadScFilter;passMuZinvSel;passMuonTrigger"      + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, "");
                std::vector<std::vector<PDS>> StackMC_Muon_LowDM_noZMassCut                = makeStackMC_DiLepton(                "passMuZinvSel"      + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, MuonWeights);
                std::vector<std::vector<PDS>> StackMC_Muon_LowDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(  "passMuZinvSel"      + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, MuonWeights);
                // bestRecoZM
                PDC dcData_Muon_LowDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Muon_LowDM_noZMassCut});
                PDC dcMC_Muon_LowDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_LowDM_noZMassCut);
                PDC dcMC_Muon_LowDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_LowDM_Normalization_noZMassCut);
                vh.push_back(PHS("DataMC_Muon_LowDM_bestRecoZM_50to250_" + cut.first + suffix,                  {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_LowDM_bestRecoZM_0to400_" + cut.first + suffix,                   {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_LowDM_Normalization_bestRecoZM_50to250_" + cut.first + suffix,    {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_LowDM_Normalization_bestRecoZM_0to400_" + cut.first + suffix,     {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            }
            for (const auto& cut : validation_cuts_high)
            {
                PDS dsData_Muon_HighDM_noZMassCut("Data",  fileMap[MuonDataset],                "Flag_eeBadScFilter;passMuZinvSel;passMuonTrigger"   + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, "");
                std::vector<std::vector<PDS>> StackMC_Muon_HighDM_noZMassCut                = makeStackMC_DiLepton(                "passMuZinvSel"   + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, MuonWeights);
                std::vector<std::vector<PDS>> StackMC_Muon_HighDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(  "passMuZinvSel"   + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned + cut.second, MuonWeights);
                // bestRecoZM
                PDC dcData_Muon_HighDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Muon_HighDM_noZMassCut});
                PDC dcMC_Muon_HighDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_HighDM_noZMassCut);
                PDC dcMC_Muon_HighDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_HighDM_Normalization_noZMassCut);
                vh.push_back(PHS("DataMC_Muon_HighDM_bestRecoZM_50to250_" + cut.first + suffix,                  {dcData_Muon_HighDM_bestRecoZM,   dcMC_Muon_HighDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_HighDM_bestRecoZM_0to400_" + cut.first + suffix,                   {dcData_Muon_HighDM_bestRecoZM,   dcMC_Muon_HighDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_HighDM_Normalization_bestRecoZM_50to250_" + cut.first + suffix,    {dcData_Muon_HighDM_bestRecoZM,   dcMC_Muon_HighDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
                vh.push_back(PHS("DataMC_Muon_HighDM_Normalization_bestRecoZM_0to400_" + cut.first + suffix,     {dcData_Muon_HighDM_bestRecoZM,   dcMC_Muon_HighDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            }
            
            // Data
            // Note: Apply Flag_eeBadScFilter to Data but not MC
            // without Z mass cut
            PDS dsData_Muon_LowDM_noZMassCut("Data",   fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSel;passMuonTrigger" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_HighDM_noZMassCut("Data",  fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSel;passMuonTrigger" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            // with on Z mass peak cut
            PDS dsData_Muon_LowDM("Data",        fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_HighDM("Data",       fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_LowDM_Loose("Data",  fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_lowDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_HighDM_Loose("Data", fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_highDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_LowDM_Mid("Data",    fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_lowDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            PDS dsData_Muon_HighDM_Mid("Data",   fileMap[MuonDataset], "Flag_eeBadScFilter;passMuZinvSelOnZMassPeak;passMuonTrigger"  + SAT_Pass_highDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, "");
            // MC
            // without Z mass cut
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM_noZMassCut                = makeStackMC_DiLepton(                 "passMuZinvSel" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM_noZMassCut               = makeStackMC_DiLepton(                 "passMuZinvSel" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM_Normalization_noZMassCut  = makeStackMC_DiLepton_Normalization(   "passMuZinvSel" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM_Normalization_noZMassCut = makeStackMC_DiLepton_Normalization(   "passMuZinvSel" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            // with on Z mass peak cut
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM  = makeStackMC_DiLepton(                               "passMuZinvSelOnZMassPeak" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM = makeStackMC_DiLepton(                               "passMuZinvSelOnZMassPeak" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM_Loose  = makeStackMC_DiLepton(                         "passMuZinvSelOnZMassPeak" + SAT_Pass_lowDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM_Loose = makeStackMC_DiLepton(                         "passMuZinvSelOnZMassPeak" + SAT_Pass_highDM_Loose + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM_Mid  = makeStackMC_DiLepton(                           "passMuZinvSelOnZMassPeak" + SAT_Pass_lowDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM_Mid = makeStackMC_DiLepton(                           "passMuZinvSelOnZMassPeak" + SAT_Pass_highDM_Mid + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_LowDM_Normalization  = makeStackMC_DiLepton_Normalization(   "passMuZinvSelOnZMassPeak" + SAT_Pass_lowDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            std::vector<std::vector<PDS>> StackMC_Muon_HighDM_Normalization = makeStackMC_DiLepton_Normalization(   "passMuZinvSelOnZMassPeak" + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drLeptonCleaned, MuonWeights);
            
            // n_jets
            PDC dcData_Muon_LowDM_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_nj( "data",   "nJets_drLeptonCleaned", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Muon_HighDM);
            
            // n_jets_Loose
            PDC dcData_Muon_LowDM_Loose_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Muon_LowDM_Loose});
            PDC dcData_Muon_HighDM_Loose_nj( "data",   "nJets_drLeptonCleaned", {dsData_Muon_HighDM_Loose});
            PDC dcMC_Muon_LowDM_Loose_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Muon_LowDM_Loose);
            PDC dcMC_Muon_HighDM_Loose_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Muon_HighDM_Loose);

            // n_jets_Mid
            PDC dcData_Muon_LowDM_Mid_nj(  "data",   "nJets_drLeptonCleaned", {dsData_Muon_LowDM_Mid});
            PDC dcData_Muon_HighDM_Mid_nj( "data",   "nJets_drLeptonCleaned", {dsData_Muon_HighDM_Mid});
            PDC dcMC_Muon_LowDM_Mid_nj(    "stack",  "nJets_drLeptonCleaned", StackMC_Muon_LowDM_Mid);
            PDC dcMC_Muon_HighDM_Mid_nj(   "stack",  "nJets_drLeptonCleaned", StackMC_Muon_HighDM_Mid);

            // HT
            PDC dcData_Muon_LowDM_ht(  "data",   "HT_drLeptonCleaned", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_ht( "data",   "HT_drLeptonCleaned", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Muon_HighDM);

            // HT_Loose
            PDC dcData_Muon_LowDM_Loose_ht(  "data",   "HT_drLeptonCleaned", {dsData_Muon_LowDM_Loose});
            PDC dcData_Muon_HighDM_Loose_ht( "data",   "HT_drLeptonCleaned", {dsData_Muon_HighDM_Loose});
            PDC dcMC_Muon_LowDM_Loose_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Muon_LowDM_Loose);
            PDC dcMC_Muon_HighDM_Loose_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Muon_HighDM_Loose);

            // HT_Mid
            PDC dcData_Muon_LowDM_Mid_ht(  "data",   "HT_drLeptonCleaned", {dsData_Muon_LowDM_Mid});
            PDC dcData_Muon_HighDM_Mid_ht( "data",   "HT_drLeptonCleaned", {dsData_Muon_HighDM_Mid});
            PDC dcMC_Muon_LowDM_Mid_ht(    "stack",  "HT_drLeptonCleaned", StackMC_Muon_LowDM_Mid);
            PDC dcMC_Muon_HighDM_Mid_ht(   "stack",  "HT_drLeptonCleaned", StackMC_Muon_HighDM_Mid);

            // met
            PDC dcData_Muon_LowDM_met(  "data",   "metWithLL", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_met( "data",   "metWithLL", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_met(    "stack",  "metWithLL", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_met(   "stack",  "metWithLL", StackMC_Muon_HighDM);
            PDC dcMC_Muon_LowDM_Normalization_met(    "stack",  "metWithLL", StackMC_Muon_LowDM_Normalization);
            PDC dcMC_Muon_HighDM_Normalization_met(   "stack",  "metWithLL", StackMC_Muon_HighDM_Normalization);

            // met_Loose
            PDC dcData_Muon_LowDM_Loose_met(  "data",   "metWithLL", {dsData_Muon_LowDM_Loose});
            PDC dcData_Muon_HighDM_Loose_met( "data",   "metWithLL", {dsData_Muon_HighDM_Loose});
            PDC dcMC_Muon_LowDM_Loose_met(    "stack",  "metWithLL", StackMC_Muon_LowDM_Loose);
            PDC dcMC_Muon_HighDM_Loose_met(   "stack",  "metWithLL", StackMC_Muon_HighDM_Loose);
            
            // met_Mid
            PDC dcData_Muon_LowDM_Mid_met(  "data",   "metWithLL", {dsData_Muon_LowDM_Mid});
            PDC dcData_Muon_HighDM_Mid_met( "data",   "metWithLL", {dsData_Muon_HighDM_Mid});
            PDC dcMC_Muon_LowDM_Mid_met(    "stack",  "metWithLL", StackMC_Muon_LowDM_Mid);
            PDC dcMC_Muon_HighDM_Mid_met(   "stack",  "metWithLL", StackMC_Muon_HighDM_Mid);
            
            // metphi
            PDC dcData_Muon_LowDM_metphi(  "data",   "metphiWithLL", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_metphi( "data",   "metphiWithLL", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_metphi(    "stack",  "metphiWithLL", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_metphi(   "stack",  "metphiWithLL", StackMC_Muon_HighDM);

            // muon pt
            PDC dcData_Muon_LowDM_MuPt1(  "data",   "cutMuPt1", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_MuPt1( "data",   "cutMuPt1", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_MuPt1(    "stack",  "cutMuPt1", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_MuPt1(   "stack",  "cutMuPt1", StackMC_Muon_HighDM);
            PDC dcData_Muon_LowDM_MuPt2(  "data",   "cutMuPt2", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_MuPt2( "data",   "cutMuPt2", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_MuPt2(    "stack",  "cutMuPt2", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_MuPt2(   "stack",  "cutMuPt2", StackMC_Muon_HighDM);

            // muon eta
            PDC dcData_Muon_LowDM_MuEta1(  "data",   "cutMuEta1", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_MuEta1( "data",   "cutMuEta1", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_MuEta1(    "stack",  "cutMuEta1", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_MuEta1(   "stack",  "cutMuEta1", StackMC_Muon_HighDM);
            PDC dcData_Muon_LowDM_MuEta2(  "data",   "cutMuEta2", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_MuEta2( "data",   "cutMuEta2", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_MuEta2(    "stack",  "cutMuEta2", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_MuEta2(   "stack",  "cutMuEta2", StackMC_Muon_HighDM);

            // bestRecoZPt
            PDC dcData_Muon_LowDM_bestRecoZPt(  "data",   "bestRecoZPt", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_bestRecoZPt( "data",   "bestRecoZPt", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_bestRecoZPt(    "stack",  "bestRecoZPt", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_bestRecoZPt(   "stack",  "bestRecoZPt", StackMC_Muon_HighDM);
            
            // bestRecoZM
            PDC dcData_Muon_LowDM_bestRecoZM(  "data",   "bestRecoZM", {dsData_Muon_LowDM_noZMassCut});
            PDC dcData_Muon_HighDM_bestRecoZM( "data",   "bestRecoZM", {dsData_Muon_HighDM_noZMassCut});
            PDC dcMC_Muon_LowDM_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_LowDM_noZMassCut);
            PDC dcMC_Muon_HighDM_bestRecoZM(   "stack",  "bestRecoZM", StackMC_Muon_HighDM_noZMassCut);
            PDC dcMC_Muon_LowDM_Normalization_bestRecoZM(    "stack",  "bestRecoZM", StackMC_Muon_LowDM_Normalization_noZMassCut);
            PDC dcMC_Muon_HighDM_Normalization_bestRecoZM(   "stack",  "bestRecoZM", StackMC_Muon_HighDM_Normalization_noZMassCut);
            
            // mtb
            PDC dcData_Muon_LowDM_mtb(  "data",   "mtb_drLeptonCleaned", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_mtb( "data",   "mtb_drLeptonCleaned", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_mtb(    "stack",  "mtb_drLeptonCleaned", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_mtb(   "stack",  "mtb_drLeptonCleaned", StackMC_Muon_HighDM);
            
            // ptb
            PDC dcData_Muon_LowDM_ptb(  "data",   "ptb_drLeptonCleaned", {dsData_Muon_LowDM});
            PDC dcData_Muon_HighDM_ptb( "data",   "ptb_drLeptonCleaned", {dsData_Muon_HighDM});
            PDC dcMC_Muon_LowDM_ptb(    "stack",  "ptb_drLeptonCleaned", StackMC_Muon_LowDM);
            PDC dcMC_Muon_HighDM_ptb(   "stack",  "ptb_drLeptonCleaned", StackMC_Muon_HighDM);
            
            // dphi
            std::vector<PDC> dcVecData_Muon_LowDM_dPhi;
            std::vector<PDC> dcVecData_Muon_HighDM_dPhi;
            std::vector<PDC> dcVecMC_Muon_LowDM_dPhi;
            std::vector<PDC> dcVecMC_Muon_HighDM_dPhi;

            for (int i = 0; i < 4; i++)
            {
                std::string var = "dPhiVec_drLeptonCleaned[" + std::to_string(i) + "]";
                dcVecData_Muon_LowDM_dPhi.push_back(    PDC("data", var, {dsData_Muon_LowDM}));
                dcVecData_Muon_HighDM_dPhi.push_back(   PDC("data", var, {dsData_Muon_HighDM}));
                dcVecMC_Muon_LowDM_dPhi.push_back(      PDC("stack", var, StackMC_Muon_LowDM));
                dcVecMC_Muon_HighDM_dPhi.push_back(     PDC("stack", var, StackMC_Muon_HighDM));
            }
                    
            // Standard selection
            vh.push_back(PHS("DataMC_Muon_LowDM_nj" + suffix,                                  {dcData_Muon_LowDM_nj,   dcMC_Muon_LowDM_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_nj" + suffix,                                 {dcData_Muon_HighDM_nj,  dcMC_Muon_HighDM_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_ht" + suffix,                                  {dcData_Muon_LowDM_ht,   dcMC_Muon_LowDM_ht},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_ht" + suffix,                                 {dcData_Muon_HighDM_ht,  dcMC_Muon_HighDM_ht}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_met" + suffix,                                 {dcData_Muon_LowDM_met,  dcMC_Muon_LowDM_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_met" + suffix,                                {dcData_Muon_HighDM_met, dcMC_Muon_HighDM_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_metphi" + suffix,                              {dcData_Muon_LowDM_metphi,  dcMC_Muon_LowDM_metphi},  {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_metphi" + suffix,                             {dcData_Muon_HighDM_metphi, dcMC_Muon_HighDM_metphi}, {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_MuPt1" + suffix,                               {dcData_Muon_LowDM_MuPt1,   dcMC_Muon_LowDM_MuPt1},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_MuPt1, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_MuPt1" + suffix,                              {dcData_Muon_HighDM_MuPt1,  dcMC_Muon_HighDM_MuPt1}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_MuPt1, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_MuPt2" + suffix,                               {dcData_Muon_LowDM_MuPt2,   dcMC_Muon_LowDM_MuPt2},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_MuPt2, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_MuPt2" + suffix,                              {dcData_Muon_HighDM_MuPt2,  dcMC_Muon_HighDM_MuPt2}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_MuPt2, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_MuEta1" + suffix,                              {dcData_Muon_LowDM_MuEta1,  dcMC_Muon_LowDM_MuEta1},  {1, 2}, "", nBins,  minEta, maxEta, true, false, label_MuEta1, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_MuEta1" + suffix,                             {dcData_Muon_HighDM_MuEta1, dcMC_Muon_HighDM_MuEta1}, {1, 2}, "", nBins,  minEta, maxEta, true, false, label_MuEta1, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_MuEta2" + suffix,                              {dcData_Muon_LowDM_MuEta2,  dcMC_Muon_LowDM_MuEta2},  {1, 2}, "", nBins,  minEta, maxEta, true, false, label_MuEta2, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_MuEta2" + suffix,                             {dcData_Muon_HighDM_MuEta2, dcMC_Muon_HighDM_MuEta2}, {1, 2}, "", nBins,  minEta, maxEta, true, false, label_MuEta2, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_bestRecoZPt" + suffix,                         {dcData_Muon_LowDM_bestRecoZPt,  dcMC_Muon_LowDM_bestRecoZPt},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_bestRecoZPt, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_bestRecoZPt" + suffix,                        {dcData_Muon_HighDM_bestRecoZPt, dcMC_Muon_HighDM_bestRecoZPt}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_bestRecoZPt, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_bestRecoZM_50to250" + suffix,                  {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_bestRecoZM_50to250" + suffix,                 {dcData_Muon_HighDM_bestRecoZM,  dcMC_Muon_HighDM_bestRecoZM}, {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_bestRecoZM_0to400" + suffix,                   {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_bestRecoZM_0to400" + suffix,                  {dcData_Muon_HighDM_bestRecoZM,  dcMC_Muon_HighDM_bestRecoZM}, {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Normalization_bestRecoZM_50to250" + suffix,    {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Normalization_bestRecoZM_50to250" + suffix,   {dcData_Muon_HighDM_bestRecoZM,  dcMC_Muon_HighDM_Normalization_bestRecoZM}, {1, 2}, "", 40, 50.0, 250.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Normalization_bestRecoZM_0to400" + suffix,     {dcData_Muon_LowDM_bestRecoZM,   dcMC_Muon_LowDM_Normalization_bestRecoZM},  {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Normalization_bestRecoZM_0to400" + suffix,    {dcData_Muon_HighDM_bestRecoZM,  dcMC_Muon_HighDM_Normalization_bestRecoZM}, {1, 2}, "", 400, 0.0, 400.0, true, false, label_bestRecoZM, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_mtb" + suffix,                                 {dcData_Muon_LowDM_mtb,   dcMC_Muon_LowDM_mtb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_mtb" + suffix,                                {dcData_Muon_HighDM_mtb,  dcMC_Muon_HighDM_mtb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_mtb, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_ptb" + suffix,                                 {dcData_Muon_LowDM_ptb,   dcMC_Muon_LowDM_ptb},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ptb, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_ptb" + suffix,                                {dcData_Muon_HighDM_ptb,  dcMC_Muon_HighDM_ptb}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ptb, "Events"));
            
            // Loose selection
            vh.push_back(PHS("DataMC_Muon_LowDM_Loose_nj" + suffix,                            {dcData_Muon_LowDM_Loose_nj,   dcMC_Muon_LowDM_Loose_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Loose_nj" + suffix,                           {dcData_Muon_HighDM_Loose_nj,  dcMC_Muon_HighDM_Loose_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Loose_ht" + suffix,                            {dcData_Muon_LowDM_Loose_ht,   dcMC_Muon_LowDM_Loose_ht},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Loose_ht" + suffix,                           {dcData_Muon_HighDM_Loose_ht,  dcMC_Muon_HighDM_Loose_ht}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Loose_met" + suffix,                           {dcData_Muon_LowDM_Loose_met,  dcMC_Muon_LowDM_Loose_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Loose_met" + suffix,                          {dcData_Muon_HighDM_Loose_met, dcMC_Muon_HighDM_Loose_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            

            // Mid selection
            vh.push_back(PHS("DataMC_Muon_LowDM_Mid_nj" + suffix,                              {dcData_Muon_LowDM_Mid_nj,   dcMC_Muon_LowDM_Mid_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Mid_nj" + suffix,                             {dcData_Muon_HighDM_Mid_nj,  dcMC_Muon_HighDM_Mid_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Mid_ht" + suffix,                              {dcData_Muon_LowDM_Mid_ht,   dcMC_Muon_LowDM_Mid_ht},  {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Mid_ht" + suffix,                             {dcData_Muon_HighDM_Mid_ht,  dcMC_Muon_HighDM_Mid_ht}, {1, 2}, "",   nBins,  minPt, maxPt, true, false, label_ht, "Events"));
            vh.push_back(PHS("DataMC_Muon_LowDM_Mid_met" + suffix,                             {dcData_Muon_LowDM_Mid_met,  dcMC_Muon_LowDM_Mid_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));
            vh.push_back(PHS("DataMC_Muon_HighDM_Mid_met" + suffix,                            {dcData_Muon_HighDM_Mid_met, dcMC_Muon_HighDM_Mid_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithLL, "Events"));

            // dphi
            for (int i = 0; i < 4; i++)
            {
                std::string nameLowDM  = "DataMC_Muon_LowDM_dPhi"  + std::to_string(i+1) + suffix;
                std::string nameHighDM = "DataMC_Muon_HighDM_dPhi" + std::to_string(i+1) + suffix;
                vh.push_back(PHS(nameLowDM,  {dcVecData_Muon_LowDM_dPhi[i], dcVecMC_Muon_LowDM_dPhi[i]},    {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
                vh.push_back(PHS(nameHighDM, {dcVecData_Muon_HighDM_dPhi[i], dcVecMC_Muon_HighDM_dPhi[i]},  {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
            }
            
            // cut flow: name, DataCollection, cutLevels
            // Plotter::CutFlowSummary::CutFlowSummary(std::string n, Plotter::DataCollection ns, std::vector<std::string> cutLevels)
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Muon_LowDM_met"  + suffix, dcData_Muon_LowDM_met,  CutLevels_Data_Muon_LowDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Muon_HighDM_met" + suffix, dcData_Muon_HighDM_met, CutLevels_Data_Muon_HighDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Muon_LowDM_met"    + suffix, dcMC_Muon_LowDM_met,    CutLevels_MC_Muon_LowDM));
            cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Muon_HighDM_met"   + suffix, dcMC_Muon_HighDM_met,   CutLevels_MC_Muon_HighDM));
        }
        
        // begin photon section
        
        // baseline selections for photon control region
        SAT_Pass_lowDM  = ";SAT_Pass_lowDM_drPhotonCleaned"  + JetPtCut;
        SAT_Pass_highDM = ";SAT_Pass_highDM_drPhotonCleaned" + JetPtCut;

        if (doDataMCPhoton)
        {
            // begin loop over photon ID selections
            for (const auto& PhotonIDSelection : PhotonIDSelections) 
            {
                //suffix
                suffix = "_" + PhotonIDSelection + JetPtCut + eraTag;
    
                // Photon 
                
                // Data
                // Note: Apply passPhotonTrigger and Flag_eeBadScFilter to Data but not MC
                PDS dsData_Photon_LowDM("Data",  fileMap[PhotonDataset], "MET_pt<250;passPhotonTrigger;Flag_eeBadScFilter;" + PhotonIDSelection + SAT_Pass_lowDM  + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned,  "");
                PDS dsData_Photon_HighDM("Data", fileMap[PhotonDataset], "MET_pt<250;passPhotonTrigger;Flag_eeBadScFilter;" + PhotonIDSelection + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned,  "");
                // MC
                // all weights
                std::string PhotonWeights = "Stop0l_trigger_eff_Photon_pt;photonSF;BTagWeight" + PrefireWeight + puWeight;
                std::vector<std::vector<PDS>> StackMC_Photon_LowDM  = makeStackMC_Photon( "MET_pt<250;" + PhotonIDSelection + SAT_Pass_lowDM  + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned, PhotonWeights);
                std::vector<std::vector<PDS>> StackMC_Photon_HighDM = makeStackMC_Photon( "MET_pt<250;" + PhotonIDSelection + SAT_Pass_highDM + Flag_ecalBadCalibFilter + semicolon_HEMVeto_drPhotonCleaned, PhotonWeights);
                
                // n_jets
                PDC dcData_Photon_LowDM_nj(  "data",   "nJets_drPhotonCleaned", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_nj( "data",   "nJets_drPhotonCleaned", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_nj(    "stack",  "nJets_drPhotonCleaned", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_nj(   "stack",  "nJets_drPhotonCleaned", StackMC_Photon_HighDM);
                
                // HT
                PDC dcData_Photon_LowDM_ht(  "data",   "HT_drPhotonCleaned", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_ht( "data",   "HT_drPhotonCleaned", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_ht(    "stack",  "HT_drPhotonCleaned", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_ht(   "stack",  "HT_drPhotonCleaned", StackMC_Photon_HighDM);
                
                // met
                PDC dcData_Photon_LowDM_met(  "data",   "metWithPhoton", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_met( "data",   "metWithPhoton", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_met(    "stack",  "metWithPhoton", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_met(   "stack",  "metWithPhoton", StackMC_Photon_HighDM);
                
                // metphi
                PDC dcData_Photon_LowDM_metphi(  "data",   "metphiWithPhoton", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_metphi( "data",   "metphiWithPhoton", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_metphi(    "stack",  "metphiWithPhoton", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_metphi(   "stack",  "metphiWithPhoton", StackMC_Photon_HighDM);
                
                // photon pt
                PDC dcData_Photon_LowDM_PhotonPt(  "data",   "cutPhotonPt", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_PhotonPt( "data",   "cutPhotonPt", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_PhotonPt(    "stack",  "cutPhotonPt", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_PhotonPt(   "stack",  "cutPhotonPt", StackMC_Photon_HighDM);
                
                // photon eta
                PDC dcData_Photon_LowDM_PhotonEta(  "data",   "cutPhotonEta", {dsData_Photon_LowDM});
                PDC dcData_Photon_HighDM_PhotonEta( "data",   "cutPhotonEta", {dsData_Photon_HighDM});
                PDC dcMC_Photon_LowDM_PhotonEta(    "stack",  "cutPhotonEta", StackMC_Photon_LowDM);
                PDC dcMC_Photon_HighDM_PhotonEta(   "stack",  "cutPhotonEta", StackMC_Photon_HighDM);
                
                // dphi
                std::vector<PDC> dcVecData_Photon_LowDM_dPhi;
                std::vector<PDC> dcVecData_Photon_HighDM_dPhi;
                std::vector<PDC> dcVecMC_Photon_LowDM_dPhi;
                std::vector<PDC> dcVecMC_Photon_HighDM_dPhi;
                for (int i = 0; i < 4; i++)
                {
                    std::string var = "dPhiVec_drPhotonCleaned[" + std::to_string(i) + "]";
                    dcVecData_Photon_LowDM_dPhi.push_back(      PDC("data", var, {dsData_Photon_LowDM}));
                    dcVecData_Photon_HighDM_dPhi.push_back(     PDC("data", var, {dsData_Photon_HighDM}));
                    dcVecMC_Photon_LowDM_dPhi.push_back(        PDC("stack", var, StackMC_Photon_LowDM));
                    dcVecMC_Photon_HighDM_dPhi.push_back(       PDC("stack", var, StackMC_Photon_HighDM));
                }
                
                vh.push_back(PHS("DataMC_Photon_LowDM_nj" + suffix,          {dcData_Photon_LowDM_nj,   dcMC_Photon_LowDM_nj},   {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_nj" + suffix,         {dcData_Photon_HighDM_nj,  dcMC_Photon_HighDM_nj},  {1, 2}, "", maxJets,  minJets,  maxJets, true, false, label_nj, "Events"));
                vh.push_back(PHS("DataMC_Photon_LowDM_ht" + suffix,          {dcData_Photon_LowDM_ht,   dcMC_Photon_LowDM_ht},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_ht" + suffix,         {dcData_Photon_HighDM_ht,  dcMC_Photon_HighDM_ht}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_ht, "Events"));
                vh.push_back(PHS("DataMC_Photon_LowDM_met" + suffix,         {dcData_Photon_LowDM_met,  dcMC_Photon_LowDM_met},  {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithPhoton, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_met" + suffix,        {dcData_Photon_HighDM_met, dcMC_Photon_HighDM_met}, {1, 2}, "", nBins,  minPt, maxPt, true, false, label_metWithPhoton, "Events"));
                vh.push_back(PHS("DataMC_Photon_LowDM_metphi" + suffix,      {dcData_Photon_LowDM_metphi,     dcMC_Photon_LowDM_metphi},  {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithPhoton, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_metphi" + suffix,     {dcData_Photon_HighDM_metphi,    dcMC_Photon_HighDM_metphi}, {1, 2}, "", nBins,  minPhi, maxPhi, true, false, label_metphiWithPhoton, "Events"));
                vh.push_back(PHS("DataMC_Photon_LowDM_PhotonPt" + suffix,    {dcData_Photon_LowDM_PhotonPt,   dcMC_Photon_LowDM_PhotonPt},  {1, 2}, "", 36, 220.0, 2020.0, true, false, label_PhotonPt, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_PhotonPt" + suffix,   {dcData_Photon_HighDM_PhotonPt,  dcMC_Photon_HighDM_PhotonPt}, {1, 2}, "", 36, 220.0, 2020.0, true, false, label_PhotonPt, "Events"));
                vh.push_back(PHS("DataMC_Photon_LowDM_PhotonEta" + suffix,   {dcData_Photon_LowDM_PhotonEta,  dcMC_Photon_LowDM_PhotonEta},  {1, 2}, "", nBins,  minEta, maxEta, true, false, label_PhotonEta, "Events"));
                vh.push_back(PHS("DataMC_Photon_HighDM_PhotonEta" + suffix,  {dcData_Photon_HighDM_PhotonEta, dcMC_Photon_HighDM_PhotonEta}, {1, 2}, "", nBins,  minEta, maxEta, true, false, label_PhotonEta, "Events"));
                
                // dphi
                for (int i = 0; i < 4; i++)
                {
                    std::string nameLowDM  = "DataMC_Photon_LowDM_dPhi"  + std::to_string(i+1) + suffix;
                    std::string nameHighDM = "DataMC_Photon_HighDM_dPhi" + std::to_string(i+1) + suffix;
                    vh.push_back(PHS(nameLowDM,  {dcVecData_Photon_LowDM_dPhi[i], dcVecMC_Photon_LowDM_dPhi[i]},    {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
                    vh.push_back(PHS(nameHighDM, {dcVecData_Photon_HighDM_dPhi[i], dcVecMC_Photon_HighDM_dPhi[i]},  {1, 2}, "", nBins,  0.0, maxPhi, true, false, vec_label_dphi[i], "Events"));
                }

                // cut flow: name, DataCollection, cutLevels
                // Plotter::CutFlowSummary::CutFlowSummary(std::string n, Plotter::DataCollection ns, std::vector<std::string> cutLevels)
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Photon_LowDM_met"      + suffix,    dcData_Photon_LowDM_met,  CutLevels_Data_Photon_LowDM));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Photon_HighDM_met"     + suffix,    dcData_Photon_HighDM_met, CutLevels_Data_Photon_HighDM));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Photon_LowDM_met"        + suffix,    dcMC_Photon_LowDM_met,    CutLevels_MC_Photon_LowDM));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Photon_HighDM_met"       + suffix,    dcMC_Photon_HighDM_met,   CutLevels_MC_Photon_HighDM));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Photon_LowDM_All_met"  + suffix,    dcData_Photon_LowDM_met,  CutLevels_Data_Photon_LowDM_All));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_Data_Photon_HighDM_All_met" + suffix,    dcData_Photon_HighDM_met, CutLevels_Data_Photon_HighDM_All));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Photon_LowDM_All_met"    + suffix,    dcMC_Photon_LowDM_met,    CutLevels_MC_Photon_LowDM_All));
                cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_Photon_HighDM_All_met"   + suffix,    dcMC_Photon_HighDM_met,   CutLevels_MC_Photon_HighDM_All));
            }
            // end loop over photon ID selections
        }
        // end photon section
    }
    // end loop over jet pt cuts 


    //lambda is your friend
    //for electrons do not use muTrigWgt (it is 0.0 for electrons)
    auto makePDSMu     = [&](const std::string& label) {return PDS("DYJetsToLL "+label, fileMap["DYJetsToLL"], "", "muTrigWgt;bTagSF_EventWeightSimple_Central;_PUweightFactor"); };
    auto makePDSElec   = [&](const std::string& label) {return PDS("DYJetsToLL "+label, fileMap["DYJetsToLL"], "", "bTagSF_EventWeightSimple_Central;_PUweightFactor"); };
    auto makePDSPhoton = [&](const std::string& label, const std::string& sample="GJets", const std::string& cuts="passPhotonSelection;HTZinv>200") {return PDS("GJets "+label, fileMap[sample], cuts, "photonAcceptanceWeight;photonEfficiencyPtWeight;photonCrossSectionRatio"); };
    
    /*
    // acceptance
    // muons
    PDC dcMC_ngenMu(            "single", "ngenMu",                     {dsDY_mu});
    PDC dcMC_ngenMatchMu(       "single", "ngenMatchMu",                {dsDY_mu});
    PDC dcMC_ngenMuInAcc(       "single", "ngenMuInAcc",                {dsDY_mu});
    PDC dcMC_ngenMatchMuInAcc(  "single", "ngenMatchMuInAcc",           {dsDY_mu});
    PDC dcMC_genMuPt(           "single", "genMu(pt)",                  {dsDY_mu});
    PDC dcMC_genMuInAccPt(      "single", "genMuInAcc(pt)",             {dsDY_mu});
    PDC dcMC_genMatchMuInAccPt( "single", "genMatchMuInAcc(pt)",        {dsDY_mu});
    PDC dcMC_genMuEta(          "single", "genMu(eta)",                 {dsDY_mu});
    PDC dcMC_genMuInAccEta(     "single", "genMuInAcc(eta)",            {dsDY_mu});
    PDC dcMC_genMatchMuInAccEta("single", "genMatchMuInAcc(eta)",       {dsDY_mu});
    // electrons
    PDC dcMC_ngenElec(            "single", "ngenElec",                 {dsDY_elec});
    PDC dcMC_ngenMatchElec(       "single", "ngenMatchElec",            {dsDY_elec});
    PDC dcMC_ngenElecInAcc(       "single", "ngenElecInAcc",            {dsDY_elec});
    PDC dcMC_ngenMatchElecInAcc(  "single", "ngenMatchElecInAcc",       {dsDY_elec});
    PDC dcMC_genElecPt(           "single", "genElec(pt)",              {dsDY_elec});
    PDC dcMC_genElecInAccPt(      "single", "genElecInAcc(pt)",         {dsDY_elec});
    PDC dcMC_genMatchElecInAccPt( "single", "genMatchElecInAcc(pt)",    {dsDY_elec});
    PDC dcMC_genElecEta(          "single", "genElec(eta)",             {dsDY_elec});
    PDC dcMC_genElecInAccEta(     "single", "genElecInAcc(eta)",        {dsDY_elec});
    PDC dcMC_genMatchElecInAccEta("single", "genMatchElecInAcc(eta)",   {dsDY_elec});
    */

    // magic lambda functions... give it pt, eta, etc
    // leptons
    // acceptance = genInAcc / gen (acceptance / MC)
    auto makePDCElecAcc_single = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genElecInAcc("+var+")", makePDSElec("e acc")}, {"genElec("+var+")", makePDSElec("e gen")}}); };
    auto makePDCMuAcc_single   = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMuInAcc("+var+")",   makePDSMu("#mu acc")}, {"genMu("+var+")",   makePDSMu("#mu gen")}}); };
    auto makePDCElecAcc_ratio = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genElecInAcc("+var+")", makePDSElec("e acc over gen")}, {"genElec("+var+")", makePDSElec("e acc over gen")}}); };
    auto makePDCMuAcc_ratio   = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMuInAcc("+var+")",   makePDSMu("#mu acc over gen")}, {"genMu("+var+")",   makePDSMu("#mu acc over gen")}}); };
    // reco efficiency = genMatchInAcc / genInAcc (reco / acceptance)
    auto makePDCElecReco_single = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchElecInAcc("+var+")", makePDSElec("e reco")}, {"genElecInAcc("+var+")", makePDSElec("e acc")}}); };
    auto makePDCMuReco_single   = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchMuInAcc("+var+")",   makePDSMu("#mu reco")}, {"genMuInAcc("+var+")",   makePDSMu("#mu acc")}}); };
    auto makePDCElecReco_ratio = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchElecInAcc("+var+")", makePDSElec("e reco over acc")}, {"genElecInAcc("+var+")", makePDSElec("e reco over acc")}}); };
    auto makePDCMuReco_ratio   = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchMuInAcc("+var+")",   makePDSMu("#mu reco over acc")}, {"genMuInAcc("+var+")",   makePDSMu("#mu reco over acc")}}); };
    // iso efficiency = genMatchIsoInAcc / genMatchInAcc (iso / reco)
    auto makePDCElecIso_single  = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchIsoElecInAcc("+var+")", makePDSElec("e iso")}, {"genMatchElecInAcc("+var+")", makePDSElec("e reco")}}); };
    auto makePDCMuIso_single    = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchIsoMuInAcc("+var+")",   makePDSMu("#mu iso")}, {"genMatchMuInAcc("+var+")",   makePDSMu("#mu reco")}}); };
    auto makePDCElecIso_ratio  = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchIsoElecInAcc("+var+")", makePDSElec("e iso over reco")}, {"genMatchElecInAcc("+var+")", makePDSElec("e iso over reco")}}); };
    auto makePDCMuIso_ratio    = [&](const std::string& var, const std::string& style) {return PDC(style, {{"genMatchIsoMuInAcc("+var+")",   makePDSMu("#mu iso over reco")}, {"genMatchMuInAcc("+var+")",   makePDSMu("#mu iso over reco")}}); };
    // photons
    // Tag must be Gen or Reco
    // acceptance = gammaLVecTagEtaPt / gammaLVecTag (acc / tag)
    auto makePDCPhotonAcc_single   = [&](const std::string& tag, const std::string& var, const std::string& style) {return PDC(style, {{"gammaLVec"+tag+"Eta("+var+")", makePDSPhoton(tag+"Eta")},           {"gammaLVec"+tag+"("+var+")", makePDSPhoton(tag)}}); };
    auto makePDCPhotonAcc_ratio    = [&](const std::string& tag, const std::string& var, const std::string& style) {return PDC(style, {{"gammaLVec"+tag+"Eta("+var+")", makePDSPhoton(tag+"Eta over "+tag)}, {"gammaLVec"+tag+"("+var+")", makePDSPhoton(tag+"Eta over "+tag)}}); };
    // iso efficiency = gammaLVecTagIso / gammaLVecTagEtaPt (iso / acc)
    // iso is only for reco!
    auto makePDCPhotonIso_single   = [&](const std::string& tag, const std::string& var, const std::string& style) {return PDC(style, {{"gammaLVec"+tag+"Iso("+var+")", makePDSPhoton(tag+"Iso")},                   {"gammaLVec"+tag+"EtaPt("+var+")", makePDSPhoton(tag+"EtaPt")}}); };
    auto makePDCPhotonIso_ratio    = [&](const std::string& tag, const std::string& var, const std::string& style) {return PDC(style, {{"gammaLVec"+tag+"Iso("+var+")", makePDSPhoton(tag+"Iso over "+tag+"EtaPt")}, {"gammaLVec"+tag+"EtaPt("+var+")", makePDSPhoton(tag+"Iso over "+tag+"EtaPt")}}); };
    // Gen matched to Reco efficiency = gammaLVecGenMatched / gammaLVecGenEtaPt  (match / acc)
    // Reco matched to Gen efficiency = gammaLVecRecoMatched / gammaLVecRecoIso  (match / iso)
    auto makePDCPhotonMatch_single = [&](const std::string& tag, const std::string& var, const std::string& style) 
    {
        if (tag.compare("Gen") == 0)    return PDC(style, {{"gammaLVec"+tag+"EtaPtMatched("+var+")", makePDSPhoton(tag+"EtaPtMatched")}, {"gammaLVec"+tag+"EtaPt("+var+")", makePDSPhoton(tag+"EtaPt")}}); 
        if (tag.compare("Reco") == 0)   return PDC(style, {{"gammaLVec"+tag+"EtaPtMatched("+var+")", makePDSPhoton(tag+"EtaPtMatched")}, {"gammaLVec"+tag+"Iso("+var+")", makePDSPhoton(tag+"Iso")}});
        printf("ERROR: Tag must be Gen or Reco for Photon matching.\n");
    };
    auto makePDCPhotonMatch_ratio  = [&](const std::string& tag, const std::string& var, const std::string& style)
    {
        if (tag.compare("Gen") == 0)    return PDC(style, {{"gammaLVec"+tag+"EtaPtMatched("+var+")", makePDSPhoton(tag+"EtaPtMatched over "+tag+"EtaPt")}, {"gammaLVec"+tag+"EtaPt("+var+")", makePDSPhoton(tag+"EtaPtMatched over "+tag+"EtaPt")}}); 
        if (tag.compare("Reco") == 0)   return PDC(style, {{"gammaLVec"+tag+"EtaPtMatched("+var+")", makePDSPhoton(tag+"EtaPtMatched over "+tag+"Iso")},   {"gammaLVec"+tag+"Iso("+var+")", makePDSPhoton(tag+"EtaPtMatched over "+tag+"Iso")}}); 
        printf("ERROR: Tag must be Gen or Reco for Photon matching.\n");
    };
    

    // photons gen: gammaLVecGen
    // photons acc: gammaLVecGenEtaPt
    // photons reco: gammaLVecGenRecoMatched
    // photons iso: gammaLVecGenIso
    /*
    PDC dcMC_genPhotonPt(     "single", "gammaLVecGenEtaPt(pt)",      {dsPhoton});
    PDC dcMC_genPhotonEta(    "single", "gammaLVecGenEtaPt(eta)",     {dsPhoton});
    PDC dcMC_matchedPhotonPt( "single", "promptPhotons(pt)",        {dsPhoton});
    PDC dcMC_matchedPhotonEta("single", "promptPhotons(eta)",       {dsPhoton});
    */

    // tops
    PDC dcMC_T1tttt("single",  "genTops(pt)", {dsT1tttt_gluino1200_lsp800, dsT1tttt_gluino1500_lsp100, dsT1tttt_gluino2000_lsp100});

    // nj                                                                                                                                                                            
    PDC dcData_SingleMuon_nj("data",   "cntNJetsPt30Eta24Zinv", {dsData_SingleMuon});
    PDC dcMC_nj(             "stack",  "cntNJetsPt30Eta24Zinv", stack_MC);
    PDC dcwMC_nj(            "stack",  "cntNJetsPt30Eta24Zinv", stackw_MC);
    PDC dcwwMC_nj(           "stack",  "cntNJetsPt30Eta24Zinv", stackww_MC);
 
    // gen Z pt: genZPt
    PDC dcData_DY_gen_pt("data",  "genZPt", {dsData_SingleMuon});
    PDC dcMC_DY_gen_pt("stack",   "genZPt", stack_MC);
  
    // reco Z pt: bestRecoZPt
    PDC dcData_DY_reco_pt("data",  "bestRecoZPt", {dsData_SingleMuon});
    PDC dcMC_DY_reco_pt("stack",   "bestRecoZPt", stack_MC);

    // gen Z eta
    PDC dcData_DY_eta("data",  "genZEta", {dsData_SingleMuon});
    PDC dcMC_DY_eta("stack",   "genZEta", stack_MC);

    // gen Z phi
    PDC dcData_DY_phi("data",  "genZPhi", {dsData_SingleMuon});
    PDC dcMC_DY_phi("stack",   "genZPhi", stack_MC);

    // invariant mass
    PDC dcData_DY_mass("data",  "genZmass", {dsData_SingleMuon});
    PDC dcMC_DY_mass("stack",   "genZmass", stack_MC);

    // met                                                                                                                                                                           
    PDC dcData_SingleMuon_met("data",   "metWithLL", {dsData_SingleMuon});
    PDC dcMC_met(             "stack",  "metWithLL", stack_MC);
    PDC dcwMC_met(            "stack",  "metWithLL", stackw_MC);
    PDC dcwwMC_met(           "stack",  "metWithLL", stackww_MC);
    // ntops                                                                                                                                                                         
    PDC dcData_SingleMuon_nt("data",   "nTopCandSortedCntZinv", {dsData_SingleMuon});
    PDC dcMC_nt(             "stack",  "nTopCandSortedCntZinv", stack_MC);
    PDC dcwMC_nt(            "stack",  "nTopCandSortedCntZinv", stackw_MC);
    PDC dcwwMC_nt(           "stack",  "nTopCandSortedCntZinv", stackww_MC);
    // MT2                                                                                                                                                                           
    PDC dcData_SingleMuon_mt2("data",   "best_had_brJet_MT2Zinv", {dsData_SingleMuon});
    PDC dcMC_mt2(             "stack",  "best_had_brJet_MT2Zinv", stack_MC);
    PDC dcwMC_mt2(            "stack",  "best_had_brJet_MT2Zinv", stackw_MC);
    PDC dcwwMC_mt2(           "stack",  "best_had_brJet_MT2Zinv", stackww_MC);
    // nb                                                                                                                                  
    PDC dcData_SingleMuon_nb("data",   "cntCSVSZinv", {dsData_SingleMuon});
    PDC dcMC_nb(             "stack",  "cntCSVSZinv", stack_MC);
    PDC dcwMC_nb(            "stack",  "cntCSVSZinv", stackw_MC);
    PDC dcwwMC_nb(           "stack",  "cntCSVSZinv", stackww_MC);
    // ht                                                                                                                                                                            
    PDC dcData_SingleMuon_ht("data",   "HTZinv", {dsData_SingleMuon});
    PDC dcMC_ht(             "stack",  "HTZinv", stack_MC);
    PDC dcwMC_ht(            "stack",  "HTZinv", stackw_MC);
    PDC dcwwMC_ht(           "stack",  "HTZinv", stackww_MC);
    // gen mu pt
    // gen mu eta
    // mu1pt                                                                                                                                                                         
    PDC dcData_SingleMuon_mu1pt("data",   "cutMuPt1", {dsData_SingleMuon});
    PDC dcMC_mu1pt(             "stack",  "cutMuPt1", stack_MC);
    PDC dcwMC_mu1pt(            "stack",  "cutMuPt1", stackw_MC);
    // mu1eta                                                                                                                                                                         
    PDC dcData_SingleMuon_mu1eta("data",   "cutMuEta1", {dsData_SingleMuon});
    PDC dcMC_mu1eta(             "stack",  "cutMuEta1", stack_MC);
    PDC dcwMC_mu1eta(            "stack",  "cutMuEta1", stackw_MC);
    // mu2pt                                                                                                                                                                         
    PDC dcData_SingleMuon_mu2pt("data",   "cutMuPt2", {dsData_SingleMuon});
    PDC dcMC_mu2pt(             "stack",  "cutMuPt2", stack_MC);
    PDC dcwMC_mu2pt(            "stack",  "cutMuPt2", stackw_MC);
    // mu2eta                                                                                                                                                                         
    PDC dcData_SingleMuon_mu2eta("data",   "cutMuEta2", {dsData_SingleMuon});
    PDC dcMC_mu2eta(             "stack",  "cutMuEta2", stack_MC);
    PDC dcwMC_mu2eta(            "stack",  "cutMuEta2", stackw_MC);
    // mll                                                                                                                                                                           
    PDC dcData_SingleMuon_mll("data",   "bestRecoZM", {dsData_SingleMuon});
    PDC dcMC_mll(             "stack",  "bestRecoZM", stack_MC);
    PDC dcwMC_mll(            "stack",  "bestRecoZM", stackw_MC);
    // nsearchbins                                                                                                                                                                   
    PDC dcData_SingleMuon_nSearchBin("data",   "nSearchBin", {dsData_SingleMuon});
    PDC dcMC_nSearchBin(             "stack",  "nSearchBin", stack_MC);
    PDC dcwMC_nSearchBin(            "stack",  "nSearchBin", stackw_MC);
    PDC dcwwMC_nSearchBin(           "stack",  "nSearchBin", stackww_MC);

    PDC dcData_SingleMuon_mht("data",   "cleanMHt", {dsData_SingleMuon});
    PDC dcwwMC_mht(           "stack",  "cleanMHt", stackww_MC);
    PDC dcwMC_mht(            "stack",  "cleanMHt", stackw_MC);


    // electron cut levels
    std::vector<std::pair<std::string,std::string>> cutlevels_electrons = {
        {"nothing",                   ""},
        //{"nosel",                     "passNoiseEventFilterZinv"},
        //{"elecZinv",                  "passNoiseEventFilterZinv;passElecZinvSel"},
    };

    // muon cut levels; commented some cut leves to make less plots
    std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
        {"nothing",                   ""},
        //{"nosel",                     "passNoiseEventFilterZinv"},
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

    // for reference from Plotter.h
    // HistSummary(std::string l, std::vector<PDC> ns, std::pair<int, int> ratio, std::string cuts, int nb, double ll, double ul, bool log, bool norm, std::string xal, std::string yal, bool isRatio = true);
    // HistSummary(std::string l, std::vector<PDC> ns, std::pair<int, int> ratio, std::string cuts, std::vector<double> be, bool log, bool norm, std::string xal, std::string yal, bool isRatio = true);

    // structs: ordering is the same as the HistSummary constructor
    struct plotStruct
    {
        std::string particle;                   // Elec, Mu, etc.
        std::string measurement;                // Acc, RecoEff, IsoEff, etc.
        std::string variable;                   // Pt, Eta, Phi, Energy, Mass, etc.
        //std::string shortVar;                 // pt, eta, phi, E, M, etc.
        PDC dataCollection; // Data Collection returned by makePDCElecAcc, makePDCElecReco, etc.
        std::string style;                      // ratio, single, etc.
        int nbins;                              // number of bins
        double xMin;                            // x min for histo
        double xMax;                            // x max for histo
        bool logBool;                           // bool for log scale (true for pt)
        bool normBool;                          // bool for norm (false)
        std::string xLabel;                     // x axis label
        std::string yLabel;                     // y axis label
        bool ratioBool;                         // bool for ratio (true)
    };

    // vector of structs

    // electrons
    // measurements: Acc, Reco, Iso
    // variables: Pt, Energy, Mass Eta, Phi
    
    std::vector<plotStruct> plotParamsElec;
    
    // Acc
    plotParamsElec.push_back({"Elec", "Acc", "Pt",     makePDCElecAcc_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_ElecPt, label_acc_single, true});
    plotParamsElec.push_back({"Elec", "Acc", "Energy", makePDCElecAcc_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_ElecEnergy, label_acc_single, true});
    plotParamsElec.push_back({"Elec", "Acc", "Mass",   makePDCElecAcc_single("M","single"),   "single", nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_acc_single, true});
    plotParamsElec.push_back({"Elec", "Acc", "Eta",    makePDCElecAcc_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_ElecEta, label_acc_single, true});
    plotParamsElec.push_back({"Elec", "Acc", "Phi",    makePDCElecAcc_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_acc_single, true});
    plotParamsElec.push_back({"Elec", "Acc", "Pt",     makePDCElecAcc_ratio("pt","ratio"),    "ratio",  nBins, minPt, maxPt, false, false, label_ElecPt, label_acc_ratio, true});
    plotParamsElec.push_back({"Elec", "Acc", "Energy", makePDCElecAcc_ratio("E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy, false, false, label_ElecEnergy, label_acc_ratio, true});
    plotParamsElec.push_back({"Elec", "Acc", "Mass",   makePDCElecAcc_ratio("M","ratio"),     "ratio",  nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_acc_ratio, true});
    plotParamsElec.push_back({"Elec", "Acc", "Eta",    makePDCElecAcc_ratio("eta","ratio"),   "ratio",  nBins, minEta, maxEta, false, false, label_ElecEta, label_acc_ratio, true});
    plotParamsElec.push_back({"Elec", "Acc", "Phi",    makePDCElecAcc_ratio("phi","ratio"),   "ratio",  nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_acc_ratio, true});
    // Reco
    plotParamsElec.push_back({"Elec", "Reco", "Pt",     makePDCElecReco_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_ElecPt, label_matched_single, true});
    plotParamsElec.push_back({"Elec", "Reco", "Energy", makePDCElecReco_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_ElecEnergy, label_matched_single, true});
    plotParamsElec.push_back({"Elec", "Reco", "Mass",   makePDCElecReco_single("M","single"),   "single", nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_matched_single, true});
    plotParamsElec.push_back({"Elec", "Reco", "Eta",    makePDCElecReco_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_ElecEta, label_matched_single, true});
    plotParamsElec.push_back({"Elec", "Reco", "Phi",    makePDCElecReco_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_matched_single, true});
    plotParamsElec.push_back({"Elec", "Reco", "Pt",     makePDCElecReco_ratio("pt","ratio"),    "ratio",  nBins, minPt, maxPt, false, false, label_ElecPt, label_matched_ratio, true});
    plotParamsElec.push_back({"Elec", "Reco", "Energy", makePDCElecReco_ratio("E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy, false, false, label_ElecEnergy, label_matched_ratio, true});
    plotParamsElec.push_back({"Elec", "Reco", "Mass",   makePDCElecReco_ratio("M","ratio"),     "ratio",  nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_matched_ratio, true});
    plotParamsElec.push_back({"Elec", "Reco", "Eta",    makePDCElecReco_ratio("eta","ratio"),   "ratio",  nBins, minEta, maxEta, false, false, label_ElecEta, label_matched_ratio, true});
    plotParamsElec.push_back({"Elec", "Reco", "Phi",    makePDCElecReco_ratio("phi","ratio"),   "ratio",  nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_matched_ratio, true});
    // Iso
    plotParamsElec.push_back({"Elec", "Iso", "Pt",     makePDCElecIso_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_ElecPt, label_iso_single, true});
    plotParamsElec.push_back({"Elec", "Iso", "Energy", makePDCElecIso_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_ElecEnergy, label_iso_single, true});
    plotParamsElec.push_back({"Elec", "Iso", "Mass",   makePDCElecIso_single("M","single"),   "single", nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_iso_single, true});
    plotParamsElec.push_back({"Elec", "Iso", "Eta",    makePDCElecIso_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_ElecEta, label_iso_single, true});
    plotParamsElec.push_back({"Elec", "Iso", "Phi",    makePDCElecIso_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_iso_single, true});
    plotParamsElec.push_back({"Elec", "Iso", "Pt",     makePDCElecIso_ratio("pt","ratio"),    "ratio",  nBins, minPt, maxPt, false, false, label_ElecPt, label_iso_ratio, true});
    plotParamsElec.push_back({"Elec", "Iso", "Energy", makePDCElecIso_ratio("E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy, false, false, label_ElecEnergy, label_iso_ratio, true});
    plotParamsElec.push_back({"Elec", "Iso", "Mass",   makePDCElecIso_ratio("M","ratio"),     "ratio",  nBins, minMassElec, maxMassElec, false, false, label_ElecMass, label_iso_ratio, true});
    plotParamsElec.push_back({"Elec", "Iso", "Eta",    makePDCElecIso_ratio("eta","ratio"),   "ratio",  nBins, minEta, maxEta, false, false, label_ElecEta, label_iso_ratio, true});
    plotParamsElec.push_back({"Elec", "Iso", "Phi",    makePDCElecIso_ratio("phi","ratio"),   "ratio",  nBins, minPhi, maxPhi, false, false, label_ElecPhi, label_iso_ratio, true});
    
    // muons
    // measurements: Acc, Reco, Iso
    // variables: Pt, Energy, Mass Eta, Phi
    
    std::vector<plotStruct> plotParamsMu;
    
    // Acc
    plotParamsMu.push_back({"Mu", "Acc", "Pt",     makePDCMuAcc_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_MuPt, label_acc_single, true});
    plotParamsMu.push_back({"Mu", "Acc", "Energy", makePDCMuAcc_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_MuEnergy, label_acc_single, true});
    plotParamsMu.push_back({"Mu", "Acc", "Mass",   makePDCMuAcc_single("M","single"),   "single", nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_acc_single, true});
    plotParamsMu.push_back({"Mu", "Acc", "Eta",    makePDCMuAcc_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_MuEta, label_acc_single, true});
    plotParamsMu.push_back({"Mu", "Acc", "Phi",    makePDCMuAcc_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_MuPhi, label_acc_single, true});
    plotParamsMu.push_back({"Mu", "Acc", "Pt",     makePDCMuAcc_ratio("pt","ratio"),   "ratio",  nBins, minPt, maxPt, false, false, label_MuPt, label_acc_ratio, true});
    plotParamsMu.push_back({"Mu", "Acc", "Energy", makePDCMuAcc_ratio("E","ratio"),    "ratio",  nBins, minEnergy, maxEnergy, false, false, label_MuEnergy, label_acc_ratio, true});
    plotParamsMu.push_back({"Mu", "Acc", "Mass",   makePDCMuAcc_ratio("M","ratio"),    "ratio",  nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_acc_ratio, true});
    plotParamsMu.push_back({"Mu", "Acc", "Eta",    makePDCMuAcc_ratio("eta","ratio"),  "ratio",  nBins, minEta, maxEta, false, false, label_MuEta, label_acc_ratio, true});
    plotParamsMu.push_back({"Mu", "Acc", "Phi",    makePDCMuAcc_ratio("phi","ratio"),  "ratio",  nBins, minPhi, maxPhi, false, false, label_MuPhi, label_acc_ratio, true});
    // Reco
    plotParamsMu.push_back({"Mu", "Reco", "Pt",     makePDCMuReco_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_MuPt, label_matched_single, true});
    plotParamsMu.push_back({"Mu", "Reco", "Energy", makePDCMuReco_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_MuEnergy, label_matched_single, true});
    plotParamsMu.push_back({"Mu", "Reco", "Mass",   makePDCMuReco_single("M","single"),   "single", nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_matched_single, true});
    plotParamsMu.push_back({"Mu", "Reco", "Eta",    makePDCMuReco_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_MuEta, label_matched_single, true});
    plotParamsMu.push_back({"Mu", "Reco", "Phi",    makePDCMuReco_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_MuPhi, label_matched_single, true});
    plotParamsMu.push_back({"Mu", "Reco", "Pt",     makePDCMuReco_ratio("pt","ratio"),   "ratio",  nBins, minPt, maxPt, false, false, label_MuPt, label_matched_ratio, true});
    plotParamsMu.push_back({"Mu", "Reco", "Energy", makePDCMuReco_ratio("E","ratio"),    "ratio",  nBins, minEnergy, maxEnergy, false, false, label_MuEnergy, label_matched_ratio, true});
    plotParamsMu.push_back({"Mu", "Reco", "Mass",   makePDCMuReco_ratio("M","ratio"),    "ratio",  nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_matched_ratio, true});
    plotParamsMu.push_back({"Mu", "Reco", "Eta",    makePDCMuReco_ratio("eta","ratio"),  "ratio",  nBins, minEta, maxEta, false, false, label_MuEta, label_matched_ratio, true});
    plotParamsMu.push_back({"Mu", "Reco", "Phi",    makePDCMuReco_ratio("phi","ratio"),  "ratio",  nBins, minPhi, maxPhi, false, false, label_MuPhi, label_matched_ratio, true});
    // Iso
    plotParamsMu.push_back({"Mu", "Iso", "Pt",     makePDCMuIso_single("pt","single"),  "single", nBins, minPt, maxPt, true, false, label_MuPt, label_iso_single, true});
    plotParamsMu.push_back({"Mu", "Iso", "Energy", makePDCMuIso_single("E","single"),   "single", nBins, minEnergy, maxEnergy, true, false, label_MuEnergy, label_iso_single, true});
    plotParamsMu.push_back({"Mu", "Iso", "Mass",   makePDCMuIso_single("M","single"),   "single", nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_iso_single, true});
    plotParamsMu.push_back({"Mu", "Iso", "Eta",    makePDCMuIso_single("eta","single"), "single", nBins, minEta, maxEta, false, false, label_MuEta, label_iso_single, true});
    plotParamsMu.push_back({"Mu", "Iso", "Phi",    makePDCMuIso_single("phi","single"), "single", nBins, minPhi, maxPhi, false, false, label_MuPhi, label_iso_single, true});
    plotParamsMu.push_back({"Mu", "Iso", "Pt",     makePDCMuIso_ratio("pt","ratio"),   "ratio",  nBins, minPt, maxPt, false, false, label_MuPt, label_iso_ratio, true});
    plotParamsMu.push_back({"Mu", "Iso", "Energy", makePDCMuIso_ratio("E","ratio"),    "ratio",  nBins, minEnergy, maxEnergy, false, false, label_MuEnergy, label_iso_ratio, true});
    plotParamsMu.push_back({"Mu", "Iso", "Mass",   makePDCMuIso_ratio("M","ratio"),    "ratio",  nBins, minMassMu, maxMassMu, false, false, label_MuMass, label_iso_ratio, true});
    plotParamsMu.push_back({"Mu", "Iso", "Eta",    makePDCMuIso_ratio("eta","ratio"),  "ratio",  nBins, minEta, maxEta, false, false, label_MuEta, label_iso_ratio, true});
    plotParamsMu.push_back({"Mu", "Iso", "Phi",    makePDCMuIso_ratio("phi","ratio"),  "ratio",  nBins, minPhi, maxPhi, false, false, label_MuPhi, label_iso_ratio, true});

    // photons 
    // measurements: Acc, Reco, Iso
    // variables: Pt, Energy, Mass Eta, Phi
    
    std::vector<plotStruct> plotParamsPhoton;
    
    // Tag must be Gen or Reco
    std::vector<std::string> tags = {"Gen", "Reco"};
    for (std::string tag : tags)
    {
      // Acc
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Pt",     makePDCPhotonAcc_single(tag, "pt","single"),  "single", nBins, minPt, maxPt,                 true,  false, label_PhotonPt,     label_map[tag+"Acc_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Energy", makePDCPhotonAcc_single(tag, "E","single"),   "single", nBins, minEnergy, maxEnergy,         true,  false, label_PhotonEnergy, label_map[tag+"Acc_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Mass",   makePDCPhotonAcc_single(tag, "M","single"),   "single", nBins, minMassPhoton, maxMassPhoton, false, false, label_PhotonMass,   label_map[tag+"Acc_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Eta",    makePDCPhotonAcc_single(tag, "eta","single"), "single", nBins, minEta, maxEta,               false, false, label_PhotonEta,    label_map[tag+"Acc_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Phi",    makePDCPhotonAcc_single(tag, "phi","single"), "single", nBins, minPhi, maxPhi,               false, false, label_PhotonPhi,    label_map[tag+"Acc_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Pt",     makePDCPhotonAcc_ratio(tag, "pt","ratio"),    "ratio",  nBins, minPt, maxPt,                 false, false, label_PhotonPt,     label_map[tag+"Acc_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Energy", makePDCPhotonAcc_ratio(tag, "E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy,         false, false, label_PhotonEnergy, label_map[tag+"Acc_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Mass",   makePDCPhotonAcc_ratio(tag, "M","ratio"),     "ratio",  nBins, minMassPhoton, maxMassPhoton, false, false, label_PhotonMass,   label_map[tag+"Acc_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Eta",    makePDCPhotonAcc_ratio(tag, "eta","ratio"),   "ratio",  nBins, minEta, maxEta,               false, false, label_PhotonEta,    label_map[tag+"Acc_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Acc", "Phi",    makePDCPhotonAcc_ratio(tag, "phi","ratio"),   "ratio",  nBins, minPhi, maxPhi,               false, false, label_PhotonPhi,    label_map[tag+"Acc_ratio"], true});
      // Iso
      // only iso for reco
      if (tag.compare("Reco") == 0)
      {
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Pt",     makePDCPhotonIso_single(tag, "pt","single"),  "single", nBins, minPt, maxPt,                     true, false,  label_PhotonPt,     label_map[tag+"Iso_single"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Energy", makePDCPhotonIso_single(tag, "E","single"),   "single", nBins, minEnergy, maxEnergy,             true, false,  label_PhotonEnergy, label_map[tag+"Iso_single"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Mass",   makePDCPhotonIso_single(tag, "M","single"),   "single", nBins, minMassPhoton, maxMassPhoton,     false, false, label_PhotonMass,   label_map[tag+"Iso_single"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Eta",    makePDCPhotonIso_single(tag, "eta","single"), "single", nBins, minEta, maxEta,                   false, false, label_PhotonEta,    label_map[tag+"Iso_single"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Phi",    makePDCPhotonIso_single(tag, "phi","single"), "single", nBins, minPhi, maxPhi,                   false, false, label_PhotonPhi,    label_map[tag+"Iso_single"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Pt",     makePDCPhotonIso_ratio(tag, "pt","ratio"),    "ratio",  nBins, minPt, maxPt,                     false, false, label_PhotonPt,     label_map[tag+"Iso_ratio"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Energy", makePDCPhotonIso_ratio(tag, "E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy,             false, false, label_PhotonEnergy, label_map[tag+"Iso_ratio"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Mass",   makePDCPhotonIso_ratio(tag, "M","ratio"),     "ratio",  nBins, minMassPhoton, maxMassPhoton,     false, false, label_PhotonMass,   label_map[tag+"Iso_ratio"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Eta",    makePDCPhotonIso_ratio(tag, "eta","ratio"),   "ratio",  nBins, minEta, maxEta,                   false, false, label_PhotonEta,    label_map[tag+"Iso_ratio"], true});
        plotParamsPhoton.push_back({"Photon", tag+"Iso", "Phi",    makePDCPhotonIso_ratio(tag, "phi","ratio"),   "ratio",  nBins, minPhi, maxPhi,                   false, false, label_PhotonPhi,    label_map[tag+"Iso_ratio"], true});
      }
      // Match
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Pt",     makePDCPhotonMatch_single(tag, "pt","single"),  "single", nBins, minPt, maxPt,                 true, false,  label_PhotonPt,     label_map[tag+"Match_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Energy", makePDCPhotonMatch_single(tag, "E","single"),   "single", nBins, minEnergy, maxEnergy,         true, false,  label_PhotonEnergy, label_map[tag+"Match_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Mass",   makePDCPhotonMatch_single(tag, "M","single"),   "single", nBins, minMassPhoton, maxMassPhoton, false, false, label_PhotonMass,   label_map[tag+"Match_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Eta",    makePDCPhotonMatch_single(tag, "eta","single"), "single", nBins, minEta, maxEta,               false, false, label_PhotonEta,    label_map[tag+"Match_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Phi",    makePDCPhotonMatch_single(tag, "phi","single"), "single", nBins, minPhi, maxPhi,               false, false, label_PhotonPhi,    label_map[tag+"Match_single"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Pt",     makePDCPhotonMatch_ratio(tag, "pt","ratio"),    "ratio",  nBins, minPt, maxPt,                 false, false, label_PhotonPt,     label_map[tag+"Match_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Energy", makePDCPhotonMatch_ratio(tag, "E","ratio"),     "ratio",  nBins, minEnergy, maxEnergy,         false, false, label_PhotonEnergy, label_map[tag+"Match_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Mass",   makePDCPhotonMatch_ratio(tag, "M","ratio"),     "ratio",  nBins, minMassPhoton, maxMassPhoton, false, false, label_PhotonMass,   label_map[tag+"Match_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Eta",    makePDCPhotonMatch_ratio(tag, "eta","ratio"),   "ratio",  nBins, minEta, maxEta,               false, false, label_PhotonEta,    label_map[tag+"Match_ratio"], true});
      plotParamsPhoton.push_back({"Photon", tag+"Match", "Phi",    makePDCPhotonMatch_ratio(tag, "phi","ratio"),   "ratio",  nBins, minPhi, maxPhi,               false, false, label_PhotonPhi,    label_map[tag+"Match_ratio"], true});
    }
    
    if (doLeptons)
    {
        // loop over electron cut levels
        for(std::pair<std::string,std::string>& cut : cutlevels_electrons)
        {
            //ratios 
            for (plotStruct& p : plotParamsElec)
            {
                if (p.style.compare("single") == 0)
                    vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style+"_"+cut.first, {p.dataCollection}, {1, 2}, cut.second, p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
                else if (p.style.compare("ratio") == 0)
                    vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style+"_"+cut.first, {p.dataCollection}, {1, 1}, cut.second, p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
            }
            
        }

        // loop over muon cut levels
        for(std::pair<std::string,std::string>& cut : cutlevels_muon)
        {
            
            //ratios 
            for (plotStruct& p : plotParamsMu)
            {
                if (p.style.compare("single") == 0)
                    vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style+"_"+cut.first, {p.dataCollection}, {1, 2}, cut.second, p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
                else if (p.style.compare("ratio") == 0)
                    vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style+"_"+cut.first, {p.dataCollection}, {1, 1}, cut.second, p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
            }
            
            //Z things (DY)
            // vh.push_back(PHS("DataMC_DY_gen_pt_"                 +cut.first,  {dcData_DY_gen_pt, dcMC_DY_gen_pt},              {1, 2}, cut.second, 60, minPt, maxPt, true, false,   "gen Z Pt",   "Events"));
            // vh.push_back(PHS("DataMC_DY_reco_pt_"                +cut.first,  {dcData_DY_reco_pt, dcMC_DY_reco_pt},            {1, 2}, cut.second, 60, minPt, maxPt, true, false,   "reco Z Pt",  "Events"));
            // vh.push_back(PHS("DataMC_DY_eta_"                    +cut.first,  {dcData_DY_eta, dcMC_DY_eta},                    {1, 2}, cut.second, 60, -10, 10, true, false,   "gen Z Eta",  "Events"));
            // vh.push_back(PHS("DataMC_DY_phi_"                    +cut.first,  {dcData_DY_phi, dcMC_DY_phi},                    {1, 2}, cut.second, 60, -10, 10, true, false,   "gen Z Phi",  "Events"));
            // vh.push_back(PHS("DataMC_DY_mass_"                   +cut.first,  {dcData_DY_mass, dcMC_DY_mass},                  {1, 2}, cut.second, 60, 0, 200, true, false,    "gen Z Mass", "Events"));
            
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
                vh.push_back(PHS("DataMC_SingleMuon_mu1pt_"          +cut.first,  {dcData_SingleMuon_mu1pt, dcMC_mu1pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_MuPt1,    "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu2pt_"          +cut.first,  {dcData_SingleMuon_mu2pt, dcMC_mu2pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_MuPt2,    "Events"));
                ///vh.push_back(PHS("DataMC_SingleMuon_genMuEta_"       +cut.first,  {dcData_SingleMuon_genMuEta, dcMC_genMuEta},     {1, 2}, cut.second, 40, -10, 10, true, false,  label_genmueta, "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu1eta_"         +cut.first,  {dcData_SingleMuon_mu1eta, dcMC_mu1eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_MuEta1,   "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu2eta_"         +cut.first,  {dcData_SingleMuon_mu2eta, dcMC_mu2eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_MuEta2,   "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mll_"            +cut.first,  {dcData_SingleMuon_mll,   dcMC_mll},             {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll,      "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_nSearchBin_"     +cut.first,  {dcData_SingleMuon_nSearchBin, dcMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin",   "Events"));
                //Shape correction weights
                vh.push_back(PHS("DataMC_SingleMuon_met_Wgt_"        +cut.first,  {dcData_SingleMuon_met,   dcwMC_met},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_met,  "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_ht_Wgt_"         +cut.first,  {dcData_SingleMuon_ht,    dcwMC_ht},              {1, 2}, cut.second, 60, 0, 1500, true, false,  label_ht,   "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_nt_Wgt_"         +cut.first,  {dcData_SingleMuon_nt,    dcwMC_nt},              {1, 2}, cut.second, 5,  0, 5,    false, false,  label_nt,  "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mt2_Wgt_"        +cut.first,  {dcData_SingleMuon_mt2,   dcwMC_mt2},             {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mt2,  "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_nb_Wgt_"         +cut.first,  {dcData_SingleMuon_nb,    dcwMC_nb},              {1, 2}, cut.second, 10, 0, 10,   true, false,  label_nb,   "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_nj_Wgt_"         +cut.first,  {dcData_SingleMuon_nj,    dcwMC_nj},              {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,   "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu1pt_Wgt_"      +cut.first,  {dcData_SingleMuon_mu1pt, dcwMC_mu1pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_MuPt1, "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu2pt_Wgt_"      +cut.first,  {dcData_SingleMuon_mu2pt, dcwMC_mu2pt},           {1, 2}, cut.second, 50, 0, 1000, true, false,  label_MuPt2, "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu1eta_Wgt_"     +cut.first,  {dcData_SingleMuon_mu1eta, dcwMC_mu1eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_MuEta1, "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mu2eta_Wgt_"     +cut.first,  {dcData_SingleMuon_mu2eta, dcwMC_mu2eta},         {1, 2}, cut.second, 40, -10, 10, true, false,  label_MuEta2, "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_mll_Wgt_"        +cut.first,  {dcData_SingleMuon_mll,   dcwMC_mll},             {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll,  "Events"));
                vh.push_back(PHS("DataMC_SingleMuon_nSearchBin_Wgt_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", "Events"));
            }
        }
    }

    if (doWeights)
    {
        for(std::pair<std::string,std::string>& cut : cutlevels_muon)
        {
            vector<double> metBins = {0, 50, 100, 150, 200, 275, 300, 350, 400, 450, 2000};
            vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
            vh.push_back(PHS("DataMCw_SingleMuon_met_"        +cut.first,  {dcData_SingleMuon_met,        dcwMC_met},        {1, 2}, cut.second, 25, 0, 1500, true, false,  label_met,     "Events / 60 GeV"));
            vh.push_back(PHS("DataMCw_SingleMuon_rebin_met_"  +cut.first,  {dcData_SingleMuon_met,        dcwMC_met},        {1, 2}, cut.second, metBins,     true, false,  label_met,     "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_ht_"         +cut.first,  {dcData_SingleMuon_ht,         dcwMC_ht},         {1, 2}, cut.second, 60, 0, 1500, true, false,  label_ht,      "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_mht_"        +cut.first,  {dcData_SingleMuon_mht,        dcwMC_mht},        {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mht,     "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nt_"         +cut.first,  {dcData_SingleMuon_nt,         dcwMC_nt},         {1, 2}, cut.second, 5,  0, 5,    true, false,  label_nt,      "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_mt2_"        +cut.first,  {dcData_SingleMuon_mt2,        dcwMC_mt2},        {1, 2}, cut.second, 60, 0, 1500, true, false,  label_mt2,     "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt2_"  +cut.first,  {dcData_SingleMuon_mt2,        dcwMC_mt2},        {1, 2}, cut.second, mt2Bins,     true, false,  label_mt2,     "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nb_"         +cut.first,  {dcData_SingleMuon_nb,         dcwMC_nb},         {1, 2}, cut.second, 7, 0, 7,   true, false,  label_nb,        "Events / bin"));
            vh.push_back(PHS("DataMCw_SingleMuon_nj_"         +cut.first,  {dcData_SingleMuon_nj,         dcwMC_nj},         {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,      "Events"));
            vh.push_back(PHS("DataMCw_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin",  "Events"));

            //// Normalization weight applied, only for blnotag selections
            if(cut.first.rfind("blnotag") == (cut.first.size()-7))
            {
                // DataMC weights applied
                vh.push_back(PHS("DataMCww_SingleMuon_met_"   +cut.first,  {dcData_SingleMuon_met,   dcwwMC_met},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met,                    "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_ht_"    +cut.first,  {dcData_SingleMuon_ht,    dcwwMC_ht},    {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,                     "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_mht_"   +cut.first,  {dcData_SingleMuon_mht,   dcwwMC_mht},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht,                    "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nt_"    +cut.first,  {dcData_SingleMuon_nt,    dcwwMC_nt},    {1, 2}, cut.second, 5,  0, 5,    true, false,  label_nt,                     "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_mt2_"   +cut.first,  {dcData_SingleMuon_mt2,   dcwwMC_mt2},   {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2,                    "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nb_"    +cut.first,  {dcData_SingleMuon_nb,    dcwwMC_nb},    {1, 2}, cut.second, 10, 0, 10,   true, false,  label_nb,                     "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nj_"    +cut.first,  {dcData_SingleMuon_nj,    dcwwMC_nj},    {1, 2}, cut.second, 30, 0, 30,   true, false,  label_nj,                     "Events"));
                vh.push_back(PHS("DataMCww_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon_nSearchBin, dcwwMC_nSearchBin}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", "Events"));
            }
        }
    }

    // photons
    if (doPhotons)
    {
        for (plotStruct& p : plotParamsPhoton)
        {
            if (p.style.compare("single") == 0)
                vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style, {p.dataCollection}, {1, 2}, "", p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
            else if (p.style.compare("ratio") == 0)
                vh.push_back(PHS("MC_"+p.particle+p.measurement+p.variable+"_"+p.style, {p.dataCollection}, {1, 1}, "", p.nbins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel, p.ratioBool));
        }
        
            
    }

    //tops
    //vh.push_back(PHS("DataMC_T1tttt", {dcMC_T1tttt}, {1, 1}, "passNoiseEventFilterZinv", 60, minPt, maxPt, true, false, label_genTopPt, "Events"));

    // Znunu
    PDS dsDY_nunu("Z#rightarrow#nu#nu", fileMap["ZJetsToNuNu"], "passLeptVeto", "");
    // nothing
    PDS dsDY_nunu_njetnorm_TriggerCentral(          "Z#rightarrow#nu#nu Trigger Central ",        fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    PDS dsDY_nunu_njetnorm_TriggerUp(               "Z#rightarrow#nu#nu Trigger Up ",             fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    PDS dsDY_nunu_njetnorm_TriggerDown(             "Z#rightarrow#nu#nu Trigger Down ",           fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    PDS dsDY_nunu_njetnorm(                         "Z#rightarrow#nu#nu Njet+norm ",              fileMap["ZJetsToNuNu"], "passLeptVeto",    "");
    // only b jet scale factor
    PDS dsDY_nunu_njetnorm_TriggerCentral_scaled(   "Z#rightarrow#nu#nu Trigger scale Central ", fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    PDS dsDY_nunu_njetnorm_TriggerUp_scaled(        "Z#rightarrow#nu#nu Trigger scale Up ",      fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    PDS dsDY_nunu_njetnorm_TriggerDown_scaled(      "Z#rightarrow#nu#nu Trigger scale Down ",    fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    PDS dsDY_nunu_njetnorm_scaled(                  "Z#rightarrow#nu#nu Njet+norm scale ",       fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central");
    // bjet scale factor, shape factor and normalization
    PDS dsDY_nunu_njetnorm_TriggerCentral_weighted( "Z#rightarrow#nu#nu Trigger weight Central ", fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffMC");
    PDS dsDY_nunu_njetnorm_TriggerUp_weighted(      "Z#rightarrow#nu#nu Trigger weight Up ",      fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffUpMC");
    PDS dsDY_nunu_njetnorm_TriggerDown_weighted(    "Z#rightarrow#nu#nu Trigger weight Down ",    fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b;TriggerEffDownMC");
    PDS dsDY_nunu_njetnorm_weighted(                "Z#rightarrow#nu#nu Njet+norm weight ",       fileMap["ZJetsToNuNu"], "passLeptVeto",    "bTagSF_EventWeightSimple_Central;njWGJets;normWgt0b");
    PDC trigger_nSearchBin( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown}, {"nSearchBin",    dsDY_nunu_njetnorm}  });
    PDC trigger_nSearchBin_scaled( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown_scaled}, {"nSearchBin",    dsDY_nunu_njetnorm_scaled}  });
    PDC trigger_nSearchBin_weighted( "single", {{"nSearchBin",    dsDY_nunu_njetnorm_TriggerCentral_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerUp_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_TriggerDown_weighted}, {"nSearchBin",    dsDY_nunu_njetnorm_weighted}  });
    
    // Znunu
    auto makePDSZnunu       = [&](const std::string& label, const std::string& cuts = "", const std::string& weights = "") {return PDS("ZJetsToNuNu "+label, fileMap["ZJetsToNuNu" + yearTag], cuts, weights); };
    auto makePDCGJetsZnunu  = [&](const std::string& var,   const std::string& style, const std::string& label, const std::string& cuts) {return PDC(style, {{var, makePDSPhoton(label, "GJets", "passPhotonSelection;" + cuts)}, {var, makePDSZnunu(label, cuts)}}); };
    
    // study jet collections and jet cleaning
    // photon jet cleaning
    std::string met_cut = "";
    //std::string met_cut = "met>250";
    std::string photon_cut = "passPhotonSelection"; 
    std::string lepton_cut = "passLeptVeto"; 
    // use met cut if it is not an empty string
    if (! met_cut.empty() )
    {
      photon_cut += ";" + met_cut;
      lepton_cut += ";" + met_cut;
    }
    auto makePDSGJets = [&](const std::string& label, const std::string& cuts) {return PDS("GJets "+label, fileMap["GJets"], cuts, ""); };
    
    // lepton jet cleaning study
    // variables: HT, MET, METPHI, dPhi, n_j, n_t, n_b 
    auto makePDSDY = [&](const std::string& label, const std::string& cuts) {return PDS("DYJetsToLL "+label, fileMap["DYJetsToLL"], cuts, ""); };
    
    if (doDYAndZnunu)
    {
        struct simplePlotStruct
        {
            std::string variable;                                       // met, HT, etc.
            std::vector<PDC> dataCollectionVector;  // vector of Data Collections
            int nBins;                                                  // number of bins
            double xMin;                                                // x min for histo
            double xMax;                                                // x max for histo
            bool logBool;                                               // bool for log scale (true for pt)
            bool normBool;                                              // bool for norm (false)
            std::string xLabel;                                         // x axis label
            std::string yLabel;                                         // y axis label
        };
        
        // map of variable names to vector of data collections
        std::map<std::string, std::vector<PDC> > dataCollectionMap;
        // vetor of variable names
        std::vector<std::string> variables = {"cntNJetsPt20Eta24", "nTopCandSortedCnt", "cntCSVS", "HT"};
        // vector of pais with tags and labels
        std::vector< std::pair<std::string, std::string> > tagVector;
        //tagVector.emplace_back("NoVeto",          "all jets");
        tagVector.emplace_back("PFLeptonCleaned", "PF lepton cleaned jets");
        tagVector.emplace_back("_drLeptonCleaned", "DR lepton cleaned jets");
        // vector of DY selections
        std::vector<std::string> selectionVec = {"Elec", "Mu"};
        // vector of styles
        std::vector<string> styleVec = {"", "_ratio"};
        // map of jets 
        std::map<std::string, std::string> jetMap;
        
        // no pt or eta cuts
        //jetMap["NoVeto"]          = "jetsLVec";
        //jetMap["PFLeptonCleaned"] = "prodJetsNoLep_jetsLVec";
        //jetMap["_drLeptonCleaned"] = "jetsLVec_drLeptonCleaned";
        // with pt and eta cuts
        jetMap["NoVeto"]          = "jetsLVec_pt20eta24";
        jetMap["PFLeptonCleaned"] = "prodJetsNoLep_jetsLVec_pt20eta24";
        jetMap["_drLeptonCleaned"] = "jetsLVec_drLeptonCleaned_pt20eta24";

        // map of y axis limits
        std::map< std::string, std::vector<float> > YAxisLimits;
        YAxisLimits["jetphi"] = {    pow(10.0, 2),      pow(10.0, 5)};
        YAxisLimits["jeteta"] = {    pow(10.0, 2),      pow(10.0, 6)};
        YAxisLimits["metphi"] = {    pow(10.0, 2),      pow(10.0, 4)};
        YAxisLimits["dphi1"]  = {5 * pow(10.0, 1),      pow(10.0, 5)};
        YAxisLimits["dphi2"]  = {5 * pow(10.0, 1),      pow(10.0, 5)};
        YAxisLimits["dphi3"]  = {5 * pow(10.0, 1),      pow(10.0, 5)};

        // plot parameters
        std::vector<simplePlotStruct> plotParamsDY;

        // try combining the following sections into one loop
        // selection, then variable, then tag (standard pattern)
        // selection, then tag, then one line per variable (unique pattern)
        
        std::string selectionLL = "";
        std::string selectionNuNu = "passBaselineNoVeto";
        // fill map
        for (const auto& s : selectionVec) 
        {
            for (const auto& variable : variables)
            {
                // z nunu
                dataCollectionMap[variable + "_" + s].emplace_back( PDC("single", variable + "NoVeto", {makePDSZnunu("all jets", selectionNuNu)} ) );
                for (const auto& tag : tagVector)
                {
                    // note that DY to LL and Z to NuNu have different selections
                    // DY to LL has dilepton selection, while Z to NuNu does not
                    selectionLL = "passBaseline" + tag.first + ";pass" + s + "ZinvSel_lowpt";
                    // DY
                    dataCollectionMap[variable + "_" + s].emplace_back(           PDC("single", variable + tag.first, {makePDSDY(tag.second, selectionLL)} ) );
                    dataCollectionMap[variable + "_" + s+ "_ratio"].emplace_back( PDC("ratio",  variable + tag.first, {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                }
            }
            
            // z nunu
            dataCollectionMap["jetpt_" + s].emplace_back(  PDC("single", jetMap["NoVeto"] + "(pt)",    {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["jeteta_" + s].emplace_back( PDC("single", jetMap["NoVeto"] + "(eta)",   {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["jetphi_" + s].emplace_back( PDC("single", jetMap["NoVeto"] + "(phi)",   {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["jetE_" + s].emplace_back(   PDC("single", jetMap["NoVeto"] + "(E)",     {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["met_" + s].emplace_back(    PDC("single", "metWithLL",                  {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["metphi_" + s].emplace_back( PDC("single", "metphiWithLL",               {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["dphi1_" + s].emplace_back(  PDC("single", "dPhiVecNoVeto[0]",           {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["dphi2_" + s].emplace_back(  PDC("single", "dPhiVecNoVeto[1]",           {makePDSZnunu("all jets", selectionNuNu)} ) );
            dataCollectionMap["dphi3_" + s].emplace_back(  PDC("single", "dPhiVecNoVeto[2]",           {makePDSZnunu("all jets", selectionNuNu)} ) );
            for (const auto& tag : tagVector)
            {
                selectionLL = "passBaseline" + tag.first + ";pass" + s + "ZinvSel_lowpt";
                // DY
                dataCollectionMap["jetpt_" + s].emplace_back(  PDC("single", jetMap[tag.first] + "(pt)",    {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["jeteta_" + s].emplace_back( PDC("single", jetMap[tag.first] + "(eta)",   {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["jetphi_" + s].emplace_back( PDC("single", jetMap[tag.first] + "(phi)",   {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["jetE_" + s].emplace_back(   PDC("single", jetMap[tag.first] + "(E)",     {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["met_" + s].emplace_back(    PDC("single", "metWithLL",                   {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["metphi_" + s].emplace_back( PDC("single", "metphiWithLL",                {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["dphi1_" + s].emplace_back(  PDC("single", "dPhiVec" + tag.first + "[0]", {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["dphi2_" + s].emplace_back(  PDC("single", "dPhiVec" + tag.first + "[1]", {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["dphi3_" + s].emplace_back(  PDC("single", "dPhiVec" + tag.first + "[2]", {makePDSDY(tag.second, selectionLL)} ) );
                dataCollectionMap["jetpt_" + s + "_ratio"].emplace_back(  PDC("ratio", jetMap[tag.first] + "(pt)",    {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["jeteta_" + s + "_ratio"].emplace_back( PDC("ratio", jetMap[tag.first] + "(eta)",   {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["jetphi_" + s + "_ratio"].emplace_back( PDC("ratio", jetMap[tag.first] + "(phi)",   {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["jetE_" + s + "_ratio"].emplace_back(   PDC("ratio", jetMap[tag.first] + "(E)",     {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["met_" + s + "_ratio"].emplace_back(    PDC("ratio", "metWithLL",                   {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["metphi_" + s + "_ratio"].emplace_back( PDC("ratio", "metphiWithLL",                {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["dphi1_" + s + "_ratio"].emplace_back(  PDC("ratio", "dPhiVec" + tag.first + "[0]", {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["dphi2_" + s + "_ratio"].emplace_back(  PDC("ratio", "dPhiVec" + tag.first + "[1]", {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
                dataCollectionMap["dphi3_" + s + "_ratio"].emplace_back(  PDC("ratio", "dPhiVec" + tag.first + "[2]", {makePDSDY(tag.second, selectionLL), makePDSZnunu(tag.second, selectionNuNu)} ) );
            }

            selectionLL = "passBaseline_drLeptonCleaned;pass" + s + "ZinvSel_lowpt";
            dataCollectionMap["dr_" + s].emplace_back( PDC("single", "dR_jetsLVec_drLeptonCleaned", {makePDSDY("all jets", selectionLL)} ) );
           
            
            // fill plot parameters
            for (const auto& style : styleVec)
            {
                std::string y_axis_label = label_Events;
                if (style.compare("_ratio") == 0)
                {
                    y_axis_label = "DYJetsToLL / ZJetsToNuNu";
                }
             
                plotParamsDY.push_back({"jetpt_" + s + style,   dataCollectionMap["jetpt_" + s + style],              nBins, 0.0, 200.0, true, false, label_jetpt, y_axis_label});
                plotParamsDY.push_back({"jeteta_" + s + style,  dataCollectionMap["jeteta_" + s + style],             nBins, minEta, maxEta, true, false, label_jeteta, y_axis_label});
                plotParamsDY.push_back({"jetphi_" + s + style,  dataCollectionMap["jetphi_" + s + style],             nBins, minPhi, maxPhi, true, false, label_jetphi, y_axis_label});
                plotParamsDY.push_back({"jetE_" + s + style,    dataCollectionMap["jetE_" + s + style],               nBins, 0.0, 2000.0, true, false, label_jetE, y_axis_label});
                plotParamsDY.push_back({"nj_" + s + style,      dataCollectionMap["cntNJetsPt20Eta24_" + s + style],  20, 0, 20, true, false, label_nj, y_axis_label});
                plotParamsDY.push_back({"nt_" + s + style,      dataCollectionMap["nTopCandSortedCnt_" + s + style],  20, 0, 20, true, false, label_nt, y_axis_label});
                plotParamsDY.push_back({"nb_" + s + style,      dataCollectionMap["cntCSVS_" + s + style],            20, 0, 20, true, false, label_nb, y_axis_label});
                plotParamsDY.push_back({"ht_" + s + style,      dataCollectionMap["HT_" + s + style],                 nBins, 0.0, 2000.0, true, false, label_ht, y_axis_label});
                plotParamsDY.push_back({"met_" + s + style,     dataCollectionMap["met_" + s + style],                nBins, 0.0, 2000.0, true, false, label_met, y_axis_label});
                plotParamsDY.push_back({"metphi_" + s + style,  dataCollectionMap["metphi_" + s + style],             nBins, minPhi, maxPhi, true, false, label_metphi, y_axis_label});
                plotParamsDY.push_back({"dphi1_" + s + style,   dataCollectionMap["dphi1_" + s + style],              nBins, 0.0, maxPhi, true, false, label_dphi1, y_axis_label});
                plotParamsDY.push_back({"dphi2_" + s + style,   dataCollectionMap["dphi2_" + s + style],              nBins, 0.0, maxPhi, true, false, label_dphi2, y_axis_label});
                plotParamsDY.push_back({"dphi3_" + s + style,   dataCollectionMap["dphi3_" + s + style],              nBins, 0.0, maxPhi, true, false, label_dphi3, y_axis_label});
            }
                plotParamsDY.push_back({"dr_" + s,      dataCollectionMap["dr_" + s],                 nBins, 0.0, 1.0, true, false, label_dr, label_Events});
        }
        
        for (const auto& p : plotParamsDY)
        {
            // check for ratio before checking for y-axis limits
            // dr and ratio options
            std::vector<std::string> splitVar = SusyUtility::getVecFromString(p.variable, '_');
            std::string var = splitVar[0];
            if (p.variable.find("dr") != std::string::npos || p.variable.find("ratio") != std::string::npos)
            {
                vh.push_back(PHS("MC_DY_" + p.variable, {p.dataCollectionVector}, {1, 1}, "", p.nBins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel));
            }
            // user defined y-axis limits
            else if (YAxisLimits.find(var) != YAxisLimits.end())
            {
                vh.push_back(PHS("MC_DY_" + p.variable, {p.dataCollectionVector}, {2, 1}, "", p.nBins, p.xMin, p.xMax, YAxisLimits[var][0], YAxisLimits[var][1], p.logBool, p.normBool, p.xLabel, p.yLabel));
            }
            // default case
            else
            {
                vh.push_back(PHS("MC_DY_" + p.variable, {p.dataCollectionVector}, {2, 1}, "", p.nBins, p.xMin, p.xMax, p.logBool, p.normBool, p.xLabel, p.yLabel));
            }
        }
    }
    
    
    // photon cuts: passPhotonSelection;passLeptVeto;HTZinv>200

    if (doGJetsAndZnunu)
    {
        // Stack HT Plot
        // GJets_HT-200To400
        // GJets_HT-400To600
        // GJets_HT-600ToInf
        
        std::string style_stack = "stack";
        std::vector<std::vector<PDS>> Photon_HT_stack_cuts = {
                                                                                   {makePDSPhoton("200 < HT < 400", "GJets_HT-200To400", "HTZinv>200;metWithPhoton>250")}, 
                                                                                   {makePDSPhoton("400 < HT < 600", "GJets_HT-400To600", "HTZinv>200;metWithPhoton>250")},
                                                                                   {makePDSPhoton("600 > HT",       "GJets_HT-600ToInf", "HTZinv>200;metWithPhoton>250")}
                                                                                 };
        std::vector<std::vector<PDS>> Photon_HT_stack_BaselineLowDM  = {
                                                                                             {makePDSPhoton("200 < HT < 400", "GJets_HT-200To400", "SAT_Pass_lowDMZinv")}, 
                                                                                             {makePDSPhoton("400 < HT < 600", "GJets_HT-400To600", "SAT_Pass_lowDMZinv")},
                                                                                             {makePDSPhoton("600 > HT",       "GJets_HT-600ToInf", "SAT_Pass_lowDMZinv")}
                                                                                           };
        std::vector<std::vector<PDS>> Photon_HT_stack_BaselineHighDM = {
                                                                                             {makePDSPhoton("200 < HT < 400", "GJets_HT-200To400", "SAT_Pass_highDMZinv")}, 
                                                                                             {makePDSPhoton("400 < HT < 600", "GJets_HT-400To600", "SAT_Pass_highDMZinv")},
                                                                                             {makePDSPhoton("600 > HT",       "GJets_HT-600ToInf", "SAT_Pass_highDMZinv")}
                                                                                           };
  
        PDC dc_GJets_ht_cuts(     "stack", "HTZinv",  Photon_HT_stack_cuts);
        PDC dc_GJets_ht_BaselineLowDM(  "stack", "HTZinv",  Photon_HT_stack_BaselineLowDM);
        PDC dc_GJets_ht_BaselineHighDM( "stack", "HTZinv",  Photon_HT_stack_BaselineHighDM);
        PDC dcMC_ZNuNu_ht_cuts(     "data",  "HTZinv",  {makePDSZnunu("HT > 200", "HTZinv>200;metWithPhoton>250")});
        PDC dcMC_ZNuNu_ht_BaselineLowDM( "data",  "HTZinv",  {makePDSZnunu("HT > 200",  "SAT_Pass_lowDMZinv")});
        PDC dcMC_ZNuNu_ht_BaselineHighDM( "data",  "HTZinv",  {makePDSZnunu("HT > 200", "SAT_Pass_highDMZinv")});
        vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_ht_" + style_stack, {dc_GJets_ht_cuts, dcMC_ZNuNu_ht_cuts},         {1, 2}, "", 100, 0.0, 2000.0, true, false, label_ht, label_Events));
        vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_ht_"   + style_stack, {dc_GJets_ht_BaselineLowDM,  dcMC_ZNuNu_ht_BaselineLowDM}, {1, 2}, "", 100, 0.0, 2000.0, true, false, label_ht, label_Events));
        vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_ht_"  + style_stack, {dc_GJets_ht_BaselineHighDM, dcMC_ZNuNu_ht_BaselineHighDM}, {1, 2}, "", 100, 0.0, 2000.0, true, false, label_ht, label_Events));
        
        std::vector<std::string> styles = {"single", "ratio"};
        // Z#rightarrow#nu#nu
        // ratio = GJets / ZJetsToNuNu
        for (const auto & style : styles)
        {
            // denominator index (2 for single, 1 for ratio)
            int d = 0;
            bool log_scale = false;
            std::string y_axis_label;
            std::string legend_label;
            if (style.compare("single") == 0)
            {
                d = 2;
                log_scale = true;
                y_axis_label = label_Events;
                legend_label = "MC";
            }
            else if (style.compare("ratio") == 0)
            {
                d = 1;
                log_scale = false;
                y_axis_label = "GJets / ZJetsToNuNu";
                legend_label = "over ZJetsToNuNu";
            }
            // use metWithPhoton instead of met (it had the photon pt added to it)
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_met_"       + style, {makePDCGJetsZnunu("metWithPhoton",                style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 80, minPt,  maxPt,  log_scale, false, label_met, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_metphi_"    + style, {makePDCGJetsZnunu("metphiWithPhoton",             style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 80, minPhi, maxPhi, log_scale, false, label_metphi, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_ht_"        + style, {makePDCGJetsZnunu("HTZinv",                       style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 80, 0.0, 2000.0,    log_scale, false, label_ht,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_nj_"        + style, {makePDCGJetsZnunu("cntNJetsPt20Eta24Zinv",        style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 10, 0, 10,          log_scale, false, label_nj,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_nb_"        + style, {makePDCGJetsZnunu("cntCSVSZinv",                  style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 10, 0, 10,          log_scale, false, label_nb,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_nt_"        + style, {makePDCGJetsZnunu("nTopCandSortedCntZinv",        style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 10, 0, 10,          log_scale, false, label_nt,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_ht200_met250_dr_"        + style, {makePDCGJetsZnunu("dR_jetsLVec_drPhotonCleaned",  style, legend_label, "HTZinv>200;metWithPhoton>250")}, {1, d}, "", 80, 0, 10.0,        log_scale, false, label_dr,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_met_"      + style, {makePDCGJetsZnunu("metWithPhoton",                style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 80, minPt,  maxPt,  log_scale, false, label_met, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_metphi_"   + style, {makePDCGJetsZnunu("metphiWithPhoton",             style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 80, minPhi, maxPhi, log_scale, false, label_metphi, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_ht_"       + style, {makePDCGJetsZnunu("HTZinv",                       style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 80, 0.0, 2000.0,    log_scale, false, label_ht,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_nj_"       + style, {makePDCGJetsZnunu("cntNJetsPt20Eta24Zinv",        style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 10, 0, 10,          log_scale, false, label_nj,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_nb_"       + style, {makePDCGJetsZnunu("cntCSVSZinv",                  style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 10, 0, 10,          log_scale, false, label_nb,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_nt_"       + style, {makePDCGJetsZnunu("nTopCandSortedCntZinv",        style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 10, 0, 10,          log_scale, false, label_nt,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineLowDM_dr_"       + style, {makePDCGJetsZnunu("dR_jetsLVec_drPhotonCleaned",  style, legend_label, "SAT_Pass_lowDMZinv")},        {1, d}, "", 80, 0, 10.0,        log_scale, false, label_dr,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_met_"     + style, {makePDCGJetsZnunu("metWithPhoton",                style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 80, minPt,  maxPt,  log_scale, false, label_met, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_metphi_"  + style, {makePDCGJetsZnunu("metphiWithPhoton",             style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 80, minPhi, maxPhi, log_scale, false, label_metphi, y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_ht_"      + style, {makePDCGJetsZnunu("HTZinv",                       style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 80, 0.0, 2000.0,    log_scale, false, label_ht,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_nj_"      + style, {makePDCGJetsZnunu("cntNJetsPt20Eta24Zinv",        style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 10, 0, 10,          log_scale, false, label_nj,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_nb_"      + style, {makePDCGJetsZnunu("cntCSVSZinv",                  style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 10, 0, 10,          log_scale, false, label_nb,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_nt_"      + style, {makePDCGJetsZnunu("nTopCandSortedCntZinv",        style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 10, 0, 10,          log_scale, false, label_nt,  y_axis_label));
            vh.push_back(PHS("MC_GJets_ZJetsToNuNu_BaselineHighDM_dr_"      + style, {makePDCGJetsZnunu("dR_jetsLVec_drPhotonCleaned",  style, legend_label, "SAT_Pass_highDMZinv")},       {1, d}, "", 80, 0, 10.0,        log_scale, false, label_dr,  y_axis_label));
        }
    }

    // search bins and validation bins
    PDC dcMC_ZNuNu_nSearchBin_LowDM("data",             "nSearchBinLowDM",            {makePDSZnunu("Search Bin Low DM",              "SAT_Pass_lowDM"           + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_low_dm;BTagWeight"          + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_nSearchBin_HighDM("data",            "nSearchBinHighDM",           {makePDSZnunu("Search Bin High DM",             "SAT_Pass_highDM"          + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_high_dm;BTagWeight"         + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_nValidationBin_LowDM("data",         "nValidationBinLowDM",        {makePDSZnunu("Validation Bin Low DM",          "SAT_Pass_lowDM"           + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_nValidationBin_LowDM_HighMET("data", "nValidationBinLowDMHighMET", {makePDSZnunu("Validation Bin Low DM High MET", "SAT_Pass_lowDM_mid_dPhi"  + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_nValidationBin_HighDM("data",        "nValidationBinHighDM",       {makePDSZnunu("Validation Bin High DM",         "SAT_Pass_highDM_mid_dPhi" + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    // MET using validation bins for cutflow
    PDC dcMC_ZNuNu_met_LowDM("data",                    "MET_pt",                     {makePDSZnunu("MET Low DM",                     "SAT_Pass_lowDM"           + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_met_LowDM_HighMET("data",            "MET_pt",                     {makePDSZnunu("MET Low DM High MET",            "SAT_Pass_lowDM_mid_dPhi"  + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    PDC dcMC_ZNuNu_met_HighDM("data",                   "MET_pt",                     {makePDSZnunu("MET High DM",                    "SAT_Pass_highDM_mid_dPhi" + Flag_ecalBadCalibFilter + semicolon_HEMVeto, "Stop0l_trigger_eff_MET_loose_baseline;BTagWeight"  + PrefireWeight + puWeight)});
    
    if (doSearchBins)
    {
        vh.push_back(PHS("ZNuNu_nSearchBin_LowDM" + eraTag,             {dcMC_ZNuNu_nSearchBin_LowDM},             {1, 1}, "", max_sb_low_dm - min_sb_low_dm,                      min_sb_low_dm,          max_sb_low_dm,          false, false,  "Search Bin Low DM", "Events", true));
        vh.push_back(PHS("ZNuNu_nSearchBin_HighDM" + eraTag,            {dcMC_ZNuNu_nSearchBin_HighDM},            {1, 1}, "", max_sb_high_dm - min_sb_high_dm,                    min_sb_high_dm,         max_sb_high_dm,         false, false,  "Search Bin High DM", "Events", true));
        vh.push_back(PHS("ZNuNu_nValidationBin_LowDM" + eraTag,         {dcMC_ZNuNu_nValidationBin_LowDM},         {1, 1}, "", max_vb_low_dm - min_vb_low_dm,                      min_vb_low_dm,          max_vb_low_dm,          false, false,  "Validation Bin Low DM", "Events", true));
        vh.push_back(PHS("ZNuNu_nValidationBin_LowDM_HighMET" + eraTag, {dcMC_ZNuNu_nValidationBin_LowDM_HighMET}, {1, 1}, "", max_vb_low_dm_high_met - min_vb_low_dm_high_met,    min_vb_low_dm_high_met, max_vb_low_dm_high_met, false, false,  "Validation Bin Low DM High MET", "Events", true));
        vh.push_back(PHS("ZNuNu_nValidationBin_HighDM" + eraTag,        {dcMC_ZNuNu_nValidationBin_HighDM},        {1, 1}, "", max_vb_high_dm - min_vb_high_dm,                    min_vb_high_dm,         max_vb_high_dm,         false, false,  "Validation Bin High DM", "Events", true));
    
        // cut flow: name, DataCollection, cutLevels
        // Plotter::CutFlowSummary::CutFlowSummary(std::string n, Plotter::DataCollection ns, std::vector<std::string> cutLevels)
        cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_ZNuNu_met_LowDM",         dcMC_ZNuNu_met_LowDM,           CutLevels_MC_ZNuNu_LowDM));
        cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_ZNuNu_met_LowDM_HighMET", dcMC_ZNuNu_met_LowDM_HighMET,   CutLevels_MC_ZNuNu_LowDM_HighMET));
        cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("CutFlow_MC_ZNuNu_met_HighDM",        dcMC_ZNuNu_met_HighDM,          CutLevels_MC_ZNuNu_HighDM));
    }
    

    set<AFS> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, sbEra, year);

    if (verbose)
    {
        std::cout << "Creating Plotter: Plotter plotter(vh, vvf, fromTuple, filename, nFiles, startFile, nEvts);" << std::endl;
        printf("    fromTuple: %s\n", fromTuple ? "true" : "false"); fflush(stdout);
        printf("    filename: %s\n", filename.c_str());              fflush(stdout);
        printf("    nFiles: %d\n", nFiles);                          fflush(stdout);
        printf("    startFile: %d\n", startFile);                    fflush(stdout);
        printf("    nEvts: %d\n", nEvts);                            fflush(stdout);
    }
  
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
