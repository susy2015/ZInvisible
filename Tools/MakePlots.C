#include "Plotter.h"
#include "samples.h"
#include "RegisterFunctions.h"
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
        {"verbose",          no_argument, 0, 'v'},
        {"filename",   required_argument, 0, 'I'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
        {"luminosity", required_argument, 0, 'L'},
        {"sbEra",      required_argument, 0, 'S'},
        {"year",       required_argument, 0, 'Y'}
    };

    bool runOnCondor        = false;
    bool doPlots = true;
    bool doSave = true;
    bool doTuple = true;
    bool fromTuple = true;
    bool verbose = false;
    string filename = "histoutput.root", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi      = -1.0;
    std::string sbEra = "SB_v1_2017";
    std::string year = "";
    while((opt = getopt_long(argc, argv, "pstfcvI:D:N:M:E:P:L:S:Y:", long_options, &option_index)) != -1)
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

        case 'S':
            sbEra = optarg;
            break;
        
        case 'Y':
            year = optarg;
            break;
        }
    }
    
    // datasets
    // year and periods
    std::string yearTag   = "_" + year; 
    // HEM veto for 2018 periods C and C
	std::string HEMcut = "";
    // PrefireWeight
    std::string PrefireWeight = "";

    // lumi for Plotter
    if (year.compare("2016") == 0)
    {
        PrefireWeight   = ";PrefireWeight";
		HEMcut = "";
    }
    else if (year.compare("2017") == 0)
    {
        PrefireWeight   = ";PrefireWeight";
		HEMcut = "";
    }
    else if (year.compare("2018") == 0)
    {
		HEMcut = "";
    }
    else if (year.compare("2018_PreHEM") == 0)
    {
		HEMcut = "";
        yearTag         = "_2018"; 
    }
    else if (year.compare("2018_PostHEM") == 0)
    {
		HEMcut = ";Pass_exHEMVeto20";
        yearTag                 = "_2018"; 
    }
    else
    {
        std::cout << "Please enter 2016, 2017, 2018, 2018_PreHEM or 2018_PostHEM for the year using the -Y option." << std::endl;
        exit(1);
    }

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

    std::cout << "input filename: " << filename << std::endl;
    std::cout << "Sample location: " << sampleloc << std::endl;

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

	for (const auto& sample : sampleList)
	{
		AnaSamples::SampleSet           ss = sample.sample_set;
		AnaSamples::SampleCollection    sc = sample.sample_collection;
		std::string                     sy = sample.sample_year; 
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


    //vector<Plotter::HistSummary> vh;
    vector<PHS> vh;


	std::string basecuts = "Pass_CaloMETRatio;Pass_JetID;";

	std::vector<Plotter::Scanner> scanners;

	string weights = "puWeight;BTagWeight"+PrefireWeight;

	string cuts = basecuts + "Pass_QCDCR" + HEMcut;

	std::set<std::string> vars = {"run", "event", "nJets", "nJets30", "PV_npvsGood", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_Stop0l", "Jet_jetId", "Jet_puId", "MET_pt", "MET_phi",
	   	"Pass_QCDCR_lowDM", "Pass_QCDCR_highDM", "Pass_trigger_MET",
		"Stop0l_trigger_eff_MET_loose_baseline", "Stop0l_trigger_eff_MET_low_dm", "Stop0l_trigger_eff_MET_high_dm",
		"Stop0l_trigger_eff_MET_low_dm_QCD", "Stop0l_trigger_eff_MET_high_dm_QCD"};

	std::set<std::string> MCvars = {"run", "event", "nJets", "nJets30", "PV_npvsGood", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_Stop0l", "Jet_jetId", "Jet_puId", "MET_pt", "MET_phi",
	   	"Pass_QCDCR_lowDM", "Pass_QCDCR_highDM", "Pass_trigger_MET",
		"Stop0l_trigger_eff_MET_loose_baseline", "Stop0l_trigger_eff_MET_low_dm", "Stop0l_trigger_eff_MET_high_dm",
		"Stop0l_trigger_eff_MET_low_dm_QCD", "Stop0l_trigger_eff_MET_high_dm_QCD",
		"Jet_corr_JEC", "Jet_jecUncertTotal", "Jet_corr_JER", "Jet_pt_jerUp", "Jet_pt_jerDown"};

	PDS dsData  = PDS("Data", fileMap["Data_MET"+yearTag], cuts, "");
	PDS dstt    = PDS("t#bar{t}", fileMap["TTbarAll"+yearTag], cuts, weights + ";ISRWeight");
	PDS dsWJets = PDS("Wjets", fileMap["WJetsToLNu"+yearTag], cuts, weights);
	PDS dsZnunu = PDS("Znunu", fileMap["ZJetsToNuNu"+yearTag], cuts, weights);
	PDS dsrare  = PDS("rare", fileMap["Rare"+yearTag], cuts, weights);
	PDS dsQCD   = PDS("QCD", fileMap["QCD"+yearTag], cuts, weights);

	string tag = "QCDCR";
	scanners.push_back(Plotter::Scanner(tag, vars, {dsData}));
	scanners.push_back(Plotter::Scanner(tag, MCvars, {dstt, dsWJets, dsZnunu, dsrare, dsQCD}));

    set<AFS> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, sbEra);

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
