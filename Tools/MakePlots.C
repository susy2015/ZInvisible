#include "Plotter.h"
#include "samples.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SusyUtility.h"
#include "TMath.h"

#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>

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
        {"verbose",          no_argument, 0, 'v'},
        {"filename",   required_argument, 0, 'I'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
        {"era",        required_argument, 0, 'Y'},
        {"syst",       required_argument, 0, 'S'}
    };
    bool runOnCondor        = false;
    bool doPlots = true;
    bool doSave = true;
    bool doTuple = true;
    bool fromTuple = true;
    bool verbose = false;
    int startFile   = 0;
    int nFiles      = -1;
    int nEvts       = -1;
    std::string filename = "histoutput.root";
    std::string dataSets = "";
    std::string sampleloc = AnaSamples::fileDir;
    std::string plotDir = "plots";
    std::string era  = "";
    std::string year = "";
    std::string suffix = "";
    std::string JESsuffix = "";
    std::string METUnClustsuffix = "";
    std::string JetPtsuffix = "";
    std::string MetPtsuffix = "";
    while((opt = getopt_long(argc, argv, "pstfcvI:D:N:M:E:P:Y:S:", long_options, &option_index)) != -1)
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

            case 'Y':
                era = optarg;
                break;

            case 'S':
                suffix = std::string("_") + optarg;
                if(suffix.find("_JES") == 0)
                {
                    JESsuffix = suffix;
                    JetPtsuffix = "_jesTotal" + suffix.substr(4);
                    MetPtsuffix = "_jesTotal" + suffix.substr(4);
                }
                if(suffix.find("_METUnClust") == 0)
                {
                    METUnClustsuffix = suffix;
                    MetPtsuffix = "_unclustEn" + suffix.substr(11);
                }
                break;
        }
    }

    // get year from era
    // if era is 2018_PreHEM, year is 2018
    if (era.length() >= 4)
    {
        year = era.substr(0, 4);
    }
    std::string earTag = "_" + era;
    std::string yearTag   = "_" + year;

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
    AnaSamples::SampleSet        SS_2016("sampleSets_PostProcessed_2016.cfg", runOnCondor, 1);
    AnaSamples::SampleSet        SS_2017("sampleSets_PostProcessed_2017.cfg", runOnCondor, 1);
    AnaSamples::SampleSet        SS_2018("sampleSets_PostProcessed_2018.cfg", runOnCondor, 1);

    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples)
    AnaSamples::SampleCollection SC_2016("sampleCollections_2016.cfg", SS_2016);
    AnaSamples::SampleCollection SC_2017("sampleCollections_2017.cfg", SS_2017);
    AnaSamples::SampleCollection SC_2018("sampleCollections_2018.cfg", SS_2018);

    // Warning: keep years together when you add them to sampleList:
    std::vector<sampleStruct> sampleList;
    sampleList.push_back({SS_2016, SC_2016, "2016"});
    sampleList.push_back({SS_2017, SC_2017, "2017"});
    sampleList.push_back({SS_2018, SC_2018, "2018"});

    map<string, vector<AFS>> fileMap;

    std::cout << "dataSets:          " << dataSets << std::endl;
    std::cout << "yearTag:           " << yearTag << std::endl;
    std::cout << "suffix:            " << suffix << std::endl;
    std::cout << "JES suffix:        " << JESsuffix << std::endl;
    std::cout << "METUnClust suffix: " << METUnClustsuffix << std::endl;
    std::cout << "JetPt suffix:      " << JetPtsuffix << std::endl;
    std::cout << "MetPt suffix:      " << MetPtsuffix << std::endl;
    for (const auto& sample : sampleList)
    {
        AnaSamples::SampleSet           ss = sample.sample_set;
        AnaSamples::SampleCollection    sc = sample.sample_collection;
        std::string                     sy = sample.sample_year;
        sc.getSampleLabels(dataSets); // This is a weird hack, but it prevents a strange but inconsequential bug -- JSW
        if(ss[dataSets] != ss.null())
        {
            fileMap[dataSets] = {ss[dataSets]};
            for(const auto& colls : ss[dataSets].getCollections())
                fileMap[colls] = {ss[dataSets]};
        }
        else if(sc[dataSets] != sc.null())
        {
            fileMap[dataSets] = {sc[dataSets]};
            int i = 0;
            for(const auto& fs : sc[dataSets])
                fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }

    std::vector<PHS> vh;

    std::vector<Plotter::Scanner> scanners;

    std::string weights = "Stop0l_evtWeight";

    std::vector<std::string> cut_vec;
    cut_vec.push_back("Pass_CaloMETRatio"+suffix);

    std::string cuts = boost::join(cut_vec, ";");

    std::map<std::string, std::string> var_name_map = {
        {"Jet_pt", "Jet_pt"+JetPtsuffix},
        {"Jet_dPhiMET", "Jet_dPhiMET"+suffix},
        {"Jet_sortedIdx", "Jet_sortedIdx"+suffix},
        {"Jet_nsortedIdx", "Jet_nsortedIdx"+suffix},
        {"Jet_btagStop0l", "Jet_btagStop0l"+JESsuffix},
        {"Jet_Stop0l", "Jet_Stop0l"+JESsuffix},
        {"MET_pt", "MET_pt"+MetPtsuffix},
        {"MET_phi", "MET_phi"+MetPtsuffix},
        {"ResolvedTop_Stop0l", "ResolvedTop_Stop0l"+JESsuffix},
        {"Pass_exHEMVeto30", "Pass_exHEMVeto30"+suffix},
        {"Pass_QCDCR", "Pass_QCDCR"+suffix},
        {"Pass_QCDCR_highDM", "Pass_QCDCR_highDM"+suffix},
        {"Pass_QCDCR_lowDM", "Pass_QCDCR_lowDM"+suffix},
        {"Pass_Baseline", "Pass_Baseline"+suffix},
        {"Pass_highDM", "Pass_highDM"+suffix},
        {"Pass_lowDM", "Pass_lowDM"+suffix},
        {"Stop0l_HT", "Stop0l_HT"+JESsuffix},
        {"Stop0l_Mtb", "Stop0l_Mtb"+suffix},
        {"Stop0l_Ptb", "Stop0l_Ptb"+JESsuffix},
        {"Stop0l_ISRJetPt", "Stop0l_ISRJetPt"+JESsuffix},
        {"Stop0l_METSig", "Stop0l_METSig"+suffix},
        {"Stop0l_nTop", "Stop0l_nTop"+JESsuffix},
        {"Stop0l_nW", "Stop0l_nW"+JESsuffix},
        {"Stop0l_nResolved", "Stop0l_nResolved"+JESsuffix},
        {"Stop0l_nbtags", "Stop0l_nbtags"+JESsuffix},
        {"Stop0l_nJets", "Stop0l_nJets"+JESsuffix},
		{"Stop0l_ResTopWeight", "Stop0l_ResTopWeight"+JESsuffix},
        {"ResolvedTopCandidate_genMatch", "ResolvedTopCandidate"+JESsuffix+"_genMatch"}
    };

    std::set<std::string> vars = {
        "run", "event",
        "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_puId", "Jet_dPhiMET",
        "Jet_sortedIdx", "Jet_nsortedIdx", "Jet_btagStop0l", "Jet_Stop0l",
        "MET_pt", "MET_phi",
		"rpseudo", "best_dphi", "best_index",
        "SB_Stop0l", "FatJet_pt", "FatJet_Stop0l", "ResolvedTop_Stop0l",
        "Pass_exHEMVeto30",
        "Pass_QCDCR", "Pass_QCDCR_highDM", "Pass_QCDCR_lowDM",
        "Pass_Baseline", "Pass_highDM", "Pass_lowDM",
        "Stop0l_HT", "Stop0l_Mtb", "Stop0l_Ptb", "Stop0l_ISRJetPt",
        "Stop0l_METSig", "Stop0l_nTop", "Stop0l_nW", "Stop0l_nResolved",
        "Stop0l_nbtags", "Stop0l_nSoftb", "Stop0l_nJets",
        "Stop0l_evtWeight",
        "nValidationBinLowDM", "nValidationBinLowDMHighMET", "nValidationBinHighDM",
        "Pass_trigger_MET",
        "Stop0l_trigger_eff_MET_loose_baseline", "Stop0l_trigger_eff_MET_loose_baseline_down", "Stop0l_trigger_eff_MET_loose_baseline_up",
        "Stop0l_trigger_eff_MET_low_dm", "Stop0l_trigger_eff_MET_low_dm_down", "Stop0l_trigger_eff_MET_low_dm_up",
        "Stop0l_trigger_eff_MET_high_dm", "Stop0l_trigger_eff_MET_high_dm_down", "Stop0l_trigger_eff_MET_high_dm_up",
        "Stop0l_trigger_eff_MET_loose_baseline_QCD", "Stop0l_trigger_eff_MET_loose_baseline_QCD_down", "Stop0l_trigger_eff_MET_loose_baseline_QCD_up",
        "Stop0l_trigger_eff_MET_low_dm_QCD", "Stop0l_trigger_eff_MET_low_dm_QCD_down", "Stop0l_trigger_eff_MET_low_dm_QCD_up",
        "Stop0l_trigger_eff_MET_high_dm_QCD", "Stop0l_trigger_eff_MET_high_dm_QCD_down", "Stop0l_trigger_eff_MET_high_dm_QCD_up"
    };

    /*std::set<std::string> vars = {"run", "event",
        "puWeight", "puWeight_Up", "puWeight_Down",
        "Jet_pt", "Jet_pt_jesTotalUp", "Jet_pt_jesTotalDown", "Jet_pt_jerUp", "Jet_pt_jerDown",
        "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_puId",
        "Jet_dPhiMET", "Jet_dPhiMET_JESUp", "Jet_dPhiMET_JESDown", "Jet_dPhiMET_METUnClustUp", "Jet_dPhiMET_METUnClustDown",
        "Jet_sortedIdx", "Jet_sortedIdx_JESUp", "Jet_sortedIdx_JESDown", "Jet_sortedIdx_METUnClustUp", "Jet_sortedIdx_METUnClustDown",
        "Jet_nsortedIdx", "Jet_nsortedIdx_JESUp", "Jet_nsortedIdx_JESDown", "Jet_nsortedIdx_METUnClustUp", "Jet_nsortedIdx_METUnClustDown",
        "Jet_btagStop0l", "Jet_btagStop0l_JESUp", "Jet_btagStop0l_JESDown",
        "Jet_Stop0l", "Jet_Stop0l_JESUp", "Jet_Stop0l_JESDown",
        "MET_pt", "MET_pt_jesTotalUp", "MET_pt_unclustEnUp", "MET_pt_jesTotalDown", "MET_pt_unclustEnDown",
        "MET_phi", "MET_phi_jesTotalUp", "MET_phi_unclustEnUp", "MET_phi_jesTotalDown", "MET_phi_unclustEnDown",
        "SB_Stop0l", "FatJet_Stop0l",
        "ResolvedTop_Stop0l", "ResolvedTop_Stop0l_JESUp", "ResolvedTop_Stop0l_JESDown",
        "Pass_exHEMVeto30", "Pass_exHEMVeto30_JESUp", "Pass_exHEMVeto30_JESDown", "Pass_exHEMVeto30_METUnClustUp", "Pass_exHEMVeto30_METUnClustDown",
        "Pass_QCDCR", "Pass_QCDCR_JESUp", "Pass_QCDCR_JESDown", "Pass_QCDCR_METUnClustUp", "Pass_QCDCR_METUnClustDown",
        "Pass_QCDCR_highDM", "Pass_QCDCR_highDM_JESUp", "Pass_QCDCR_highDM_JESDown", "Pass_QCDCR_highDM_METUnClustUp", "Pass_QCDCR_highDM_METUnClustDown",
        "Pass_QCDCR_lowDM", "Pass_QCDCR_lowDM_JESUp", "Pass_QCDCR_lowDM_JESDown", "Pass_QCDCR_lowDM_METUnClustUp", "Pass_QCDCR_lowDM_METUnClustDown",
        "Pass_Baseline", "Pass_Baseline_JESUp", "Pass_Baseline_JESDown", "Pass_Baseline_METUnClustUp", "Pass_Baseline_METUnClustDown",
        "Pass_highDM", "Pass_highDM_JESUp", "Pass_highDM_JESDown", "Pass_highDM_METUnClustUp", "Pass_highDM_METUnClustDown",
        "Pass_lowDM", "Pass_lowDM_JESUp", "Pass_lowDM_JESDown", "Pass_lowDM_METUnClustUp", "Pass_lowDM_METUnClustDown",
        "Stop0l_HT", "Stop0l_HT_JESUp", "Stop0l_HT_JESDown",
        "Stop0l_Mtb", "Stop0l_Mtb_JESUp", "Stop0l_Mtb_JESDown", "Stop0l_Mtb_METUnClustUp", "Stop0l_Mtb_METUnClustDown",
        "Stop0l_Ptb", "Stop0l_Ptb_JESUp", "Stop0l_Ptb_JESDown",
        "Stop0l_ISRJetPt", "Stop0l_ISRJetPt_JESUp", "Stop0l_ISRJetPt_JESDown",
        "Stop0l_METSig", "Stop0l_METSig_JESUp", "Stop0l_METSig_JESDown", "Stop0l_METSig_METUnClustUp", "Stop0l_METSig_METUnClustDown",
        "Stop0l_nTop", "Stop0l_nTop_JESUp", "Stop0l_nTop_JESDown",
        "Stop0l_nW", "Stop0l_nW_JESUp", "Stop0l_nW_JESDown",
        "Stop0l_nResolved", "Stop0l_nResolved_JESUp", "Stop0l_nResolved_JESDown",
        "nResolvedTopCandidate", "nResolvedTopCandidate_JESUp", "nResolvedTopCandidate_JESDown",
        "Stop0l_nbtags", "Stop0l_nbtags_JESUp", "Stop0l_nbtags_JESDown",
        "Stop0l_nSoftb",
        "Stop0l_nJets", "Stop0l_nJets_JESUp", "Stop0l_nJets_JESDown",
        "Stop0l_evtWeight",
        "nValidationBinLowDM", "nValidationBinLowDMHighMET", "nValidationBinHighDM",
        "Pass_trigger_MET",
        "Stop0l_trigger_eff_MET_loose_baseline", "Stop0l_trigger_eff_MET_loose_baseline_down", "Stop0l_trigger_eff_MET_loose_baseline_up",
        "Stop0l_trigger_eff_MET_low_dm", "Stop0l_trigger_eff_MET_low_dm_down", "Stop0l_trigger_eff_MET_low_dm_up",
        "Stop0l_trigger_eff_MET_high_dm", "Stop0l_trigger_eff_MET_high_dm_down", "Stop0l_trigger_eff_MET_high_dm_up",
        "Stop0l_trigger_eff_MET_loose_baseline_QCD", "Stop0l_trigger_eff_MET_loose_baseline_QCD_down", "Stop0l_trigger_eff_MET_loose_baseline_QCD_up",
        "Stop0l_trigger_eff_MET_low_dm_QCD", "Stop0l_trigger_eff_MET_low_dm_QCD_down", "Stop0l_trigger_eff_MET_low_dm_QCD_up",
        "Stop0l_trigger_eff_MET_high_dm_QCD", "Stop0l_trigger_eff_MET_high_dm_QCD_down", "Stop0l_trigger_eff_MET_high_dm_QCD_up"
    };*/

    std::set<std::string> MCvars = {
        "Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_partonFlavour", "GenJet_hadronFlavour",
		"rjet", "best_genindex", "best_dR",
        "genWeight",
        "LHE_HTIncoming",
        "SB_SF", "SB_SFerr",
        "puWeight", "puWeight_Up", "puWeight_Down",
        "BTagWeight", "BTagWeight_Up", "BTagWeight_Down",
        "ISRWeight", "ISRWeight_Up", "ISRWeight_Down",
        "pdfWeight_Up", "pdfWeight_Down",
		"Stop0l_ResTopWeight", "Stop0l_ResTopWeight_Up", "Stop0l_ResTopWeight_Dn",
		"Stop0l_DeepAK8_SFWeight",
		"Stop0l_DeepAK8_SFWeight_total_up", "Stop0l_DeepAK8_SFWeight_total_dn",
		"Stop0l_DeepAK8_SFWeight_top_up", "Stop0l_DeepAK8_SFWeight_top_dn",
		"Stop0l_DeepAK8_SFWeight_w_up", "Stop0l_DeepAK8_SFWeight_w_dn",
		"Stop0l_DeepAK8_SFWeight_veto_up", "Stop0l_DeepAK8_SFWeight_veto_dn",
        "nWplusLep", "nWminusLep", "TopPt", "TbarPt", "LHEScaleWeight"
    };

    if (year == "2017")
    {
        vars.insert("Flag_ecalBadCalibFilterV2");
        MCvars.insert("17BtoEpuWeight");
        MCvars.insert("17BtoEpuWeight_Up");
        MCvars.insert("17BtoEpuWeight_Down");
        MCvars.insert("17FpuWeight");
        MCvars.insert("17FpuWeight_Up");
        MCvars.insert("17FpuWeight_Down");
    }

    if (year != "2018")
    {
        MCvars.insert("PrefireWeight");
        MCvars.insert("PrefireWeight_Up");
        MCvars.insert("PrefireWeight_Down");
    }

    std::map<std::string, std::vector<std::string>> region_to_cuts = {
        {"QCDCR", {"Pass_QCDCR"+suffix}},
        {"Baseline", {"Pass_Baseline"+suffix}},
        {"vLowDM", {"Pass_lowDM"+suffix}},
        {"vLowDMHighMET", {"Pass_EventFilter"+suffix,
                           "Pass_JetID"+suffix,
                           "Pass_LeptonVeto"+suffix,
                           "Pass_NJets30"+suffix,
                           "Pass_MET"+suffix,
                           "Pass_HT"+suffix,
                           "Pass_dPhiMETMedDM"+suffix,
                           "Stop0l_nTop" + JESsuffix + "==0",
                           "Stop0l_nW" + JESsuffix + "==0",
                           "Stop0l_nResolved" + JESsuffix + "==0",
                           "Stop0l_Mtb" + suffix + "<175",
                           "Stop0l_ISRJetPt" + JESsuffix + ">=200",
                           "Stop0l_METSig" + suffix + ">10"}},
        {"vHighDM", {"Pass_EventFilter" + suffix,
                     "Pass_JetID" + suffix,
                     "Pass_LeptonVeto" + suffix,
                     "Pass_NJets30" + suffix,
                     "Pass_MET" + suffix,
                     "Pass_HT" + suffix,
                     "Stop0l_nbtags" + JESsuffix + ">=1",
                     "Stop0l_nJets" + JESsuffix + ">=5",
                     "!Pass_dPhiMETHighDM" + suffix,
                     "Pass_dPhiMETLowDM" + suffix}}};

    MCvars.insert(vars.begin(), vars.end());

    //std::set<std::string> QCDvars = {"bootstrapWeight", "bootstrapWeight_smaller", "Stop0l_smearWeight", "uniqueID"};
    std::set<std::string> QCDvars = {"uniqueID"};
    QCDvars.insert(MCvars.begin(), MCvars.end());

    std::set<std::string> sQCDvars = {"Stop0l_smearWeight", "uniqueID"};
    sQCDvars.insert(MCvars.begin(), MCvars.end());

	std::set<std::string> Topvars = {"Stop0l_topptWeight", "Stop0l_topMGPowWeight", "Stop0l_topptOnly", "Stop0l_topptOnly_Up", "Stop0l_topptOnly_Down"};
	Topvars.insert(MCvars.begin(), MCvars.end());

    std::string region_cuts;
    for (string region : {"QCDCR", "Baseline", "vLowDM", "vLowDMHighMET", "vHighDM"})
    {
        region_cuts = boost::join(region_to_cuts[region], ";");
        PDS dsData        = PDS("Data",         fileMap["Data_MET_2016"        ], cuts + ";" + region_cuts, "");
        //PDS dsDataPreHEM  = PDS("Data_PreHEM",  fileMap["Data_MET_2018_PreHEM" ], cuts + ";" + region_cuts, "");
        //PDS dsDataPostHEM = PDS("Data_PostHEM", fileMap["Data_MET_2018_PostHEM"], cuts + ";" + region_cuts, "");
        PDS dsDataBE      = PDS("Data_BE",      fileMap["Data_MET_2017_BE"     ], cuts + ";" + region_cuts, "");
        PDS dsDataF       = PDS("Data_F",       fileMap["Data_MET_2017_F"      ], cuts + ";" + region_cuts, "");
        PDS dsData18      = PDS("Data18",       fileMap["Data_MET_2018"        ], cuts + ";" + region_cuts, "");
        PDS dstt          = PDS("ttbar",        fileMap["TTbar"       + yearTag], cuts + ";" + region_cuts, weights);
        PDS dstt2         = PDS("ttbarnohad",   fileMap["TTbarNoHad"  + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsttHT        = PDS("ttbarHT",      fileMap["TTbarHT"     + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsWJets       = PDS("Wjets",        fileMap["WJetsToLNu"  + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsZnunu       = PDS("Znunu",        fileMap["ZJetsToNuNu" + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsrare        = PDS("rare",         fileMap["Rare"        + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsQCD         = PDS("QCD",          fileMap["QCD"         + yearTag], cuts + ";" + region_cuts, weights);
        PDS dsQCD_s       = PDS("QCD_smear",    fileMap["QCD_Smear"   + yearTag], cuts + ";" + region_cuts, weights);

        scanners.push_back(Plotter::Scanner(region, vars, {dsData, dsData18, dsDataBE, dsDataF}, var_name_map));
        scanners.push_back(Plotter::Scanner(region, MCvars, {dsWJets, dsZnunu, dsrare}, var_name_map));
        scanners.push_back(Plotter::Scanner(region, Topvars, {dstt, dstt2, dsttHT}, var_name_map));
        scanners.push_back(Plotter::Scanner(region, MCvars, {dsQCD}, var_name_map));
        scanners.push_back(Plotter::Scanner(region, sQCDvars, {dsQCD_s}, var_name_map));
    }

    set<AFS> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, year, var_name_map);

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
    plotter.setScanners(scanners);
    //plotter.setCutFlows(cutFlowSummaries);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
