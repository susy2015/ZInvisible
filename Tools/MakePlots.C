#include "Plotter.h"
#include "samples.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SusyUtility.h"
#include "TMath.h"

#include <getopt.h>
#include <iostream>
#include <algorithm>

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
        {"luminosity", required_argument, 0, 'L'},
        {"refLumi",    required_argument, 0, 'R'},
        {"era",        required_argument, 0, 'Y'}
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
    double lumi     = -1.0;
    std::string refLumi = "";
    std::string filename = "histoutput.root";
    std::string dataSets = "";
    std::string sampleloc = AnaSamples::fileDir;
    std::string plotDir = "plots";
    std::string era  = "";
    std::string year = "";
    while((opt = getopt_long(argc, argv, "pstfcvI:D:N:M:E:P:L:R:Y:", long_options, &option_index)) != -1)
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

        case 'L':
            lumi = atof(optarg);
            break;
        
        case 'R':
            refLumi = optarg;
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
    // year and periods
    std::string eraTag    = "_" + era; 
    std::string yearTag   = "_" + year; 
    std::string periodTag = ""; 
    // HEM veto for 2018 periods C and C
    std::string HEMVeto = "";
    std::string semicolon_HEMVeto = "";
    // Note: Don't apply Flag_ecalBadCalibFilter to 2016, but apply it to 2017 and 2018
    std::string Flag_ecalBadCalibFilter = "";
    // PrefireWeight
    std::string PrefireWeight = "";
	std::string ISRWeight = "";
    std::string puWeight                = ";puWeight";

    // lumi for Plotter
    if (era.compare("2016") == 0)
    {
        PrefireWeight   = ";PrefireWeight";
		ISRWeight = ";ISRWeight";
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
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2018_PreHEM") == 0)
    {
        periodTag               = "_PreHEM";
        Flag_ecalBadCalibFilter = ";Flag_ecalBadCalibFilter";
    }
    else if (era.compare("2018_PostHEM") == 0)
    {
        // HEM vetos: use ";veto_name" so that it can be appended to cuts
        HEMVeto                           = "SAT_Pass_HEMVeto";
        semicolon_HEMVeto                 = ";" + HEMVeto;
        periodTag                         = "_PostHEM";
        Flag_ecalBadCalibFilter           = ";Flag_ecalBadCalibFilter";
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
    
    map<string, vector<AFS>> fileMap;

	std::cout << "dataSets: " << dataSets << std::endl;
	std::cout << "yearTag: " << yearTag << std::endl;
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

<<<<<<< HEAD

    //vector<Plotter::HistSummary> vh;
    vector<PHS> vh;

	std::vector<Plotter::Scanner> scanners;

	string weights = "BTagWeight;genWeight" + PrefireWeight + puWeight;

	string cuts = "Pass_CaloMETRatio;Pass_LeptonVeto;Flag_goodVertices;Flag_HBHENoiseFilter;Flag_HBHENoiseIsoFilter;Flag_BadPFMuonFilter;Flag_globalSuperTightHalo2016Filter;Flag_eeBadScFilter;Pass_JetID;Pass_NJets20;Pass_MET;Pass_HT" + semicolon_HEMVeto;

	std::set<std::string> vars = {"run", "event",
	   	"PV_npvsGood",
	   	"Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_puId", "Jet_dPhiMET", "Jet_order",
	   	"MET_pt", "MET_phi",
	   	"Pass_QCDCR", "Pass_QCDCR_lowDM", "Pass_QCDCR_highDM",
		"Pass_LeptonVeto", "Pass_CaloMETRatio", "Pass_JetID", "Pass_EventFilter",
	   	"Pass_trigger_MET",
		"Jet_neEmEF", "Jet_neHEF", "Jet_chEmEF", "Jet_chHEF", "Jet_muEF",
	   	"Stop0l_HT", "Stop0l_Mtb", "Stop0l_Ptb", "Stop0l_ISRJetPt", "Stop0l_METSig",
	   	"Stop0l_nTop", "Stop0l_nW", "Stop0l_nResolved", "Stop0l_nbtags", "Stop0l_nSoftb",
		"Stop0l_evtWeight",
		"SAT_Pass_HEMVeto20", "SAT_Pass_HEMVeto30", "Pass_HEMVeto20", "Pass_HEMVeto30", "Pass_exHEMVeto20", "Pass_exHEMVeto30",
		"nSearchBinLowDM", "nSearchBinHighDM", "nSearchBinHighDMLoose",
		"nValidationBinLowDM", "nValidationBinLowDMHighMET", "nValidationBinHighDM",      
		"nBottoms", "nSoftBottoms", "nMergedTops", "nJets", "nWs", "nResolvedTops", "HT", "ptb", "mtb", "ISRJetPt",
		"Flag_BadChargedCandidateFilter", "Flag_BadChargedCandidateSummer16Filter", "Flag_ecalBadCalibFilter",
		"Flag_goodVertices", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
		"Flag_BadPFMuonFilter", "Flag_globalSuperTightHalo2016Filter", "Flag_eeBadScFilter",
		"Stop0l_trigger_eff_MET_loose_baseline", "Stop0l_trigger_eff_MET_low_dm", "Stop0l_trigger_eff_MET_high_dm",
		"Stop0l_trigger_eff_MET_low_dm_QCD", "Stop0l_trigger_eff_MET_high_dm_QCD",
		"Pass_Baseline", "Pass_highDM", "Pass_lowDM"
	};

	if (era == "2017")
	{
		vars.insert("Flag_ecalBadCalibFilterV2");
	}

	std::string datacuts = cuts + ";!Pass_Baseline"; // Do this to keep the signal region blinded

	PDS dsData  = PDS("Data", fileMap["Data_MET"+eraTag], datacuts, "");
	PDS dstt    = PDS("ttbar", fileMap["TTbar"+yearTag], cuts, weights + ISRWeight);
	PDS dstt2   = PDS("ttbarnohad", fileMap["TTbarNoHad"+yearTag], cuts, weights + ISRWeight);
	PDS dsWJets = PDS("Wjets", fileMap["WJetsToLNu"+yearTag], cuts, weights);
	PDS dsZnunu = PDS("Znunu", fileMap["ZJetsToNuNu"+yearTag], cuts, weights);
	PDS dsrare  = PDS("rare", fileMap["Rare"+yearTag], cuts, weights);
	PDS dsQCD   = PDS("QCD", fileMap["QCD"+yearTag], cuts, weights);
	PDS dsQCD_s = PDS("QCD_smear", fileMap["QCD_smear"+yearTag], cuts, weights + ";Stop0l_evtWeight");
	// This weight is a dirty hack for the smeared QCD.  The smearing process alters the effective cross section, and this is
	// encoded in "Stop0l_evtWeight", which normally only includes the cross section divided by the number of generated events
	// Because the Plotter framework gets the cross section and number of generated events (and luminosity) from the config files,
	// I usually don't use "Stop0l_evtWeight", but in this case I need to.  In order to avoid double-counting the cross section
	// and number of generated events, I set those to 1 in the config file

	string tag = "QCDCR";
	scanners.push_back(Plotter::Scanner(tag, vars, {dsData, dstt, dstt2, dsWJets, dsZnunu, dsrare, dsQCD, dsQCD_s}));

	/*
	std::string SRcuts = "Pass_Baseline" + semicolon_HEMVeto;

	PDS srQCD   = PDS("QCD", fileMap["QCD"+yearTag], SRcuts, weights);

	scanners.push_back(Plotter::Scanner("SR", vars, {srQCD}));
	*/

    set<AFS> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, year);

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
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
