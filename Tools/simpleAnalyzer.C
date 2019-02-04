#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SATException.h"

#include "TopTaggerTools/Tools/include/HistoContainer.h"

#include "CleanedJets.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "BTagCorrector.h"
#include "TTbarCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"
#include "customize.h"
//#include "RegisterFunctions.h"

#include "TopTaggerResults.h"
#include "Constituent.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"


std::set<AnaSamples::FileSummary> getFS(const std::string& dataSets, const bool& isCondor)
{
    // follow this syntax; order matters for your arguments
    //SampleSet::SampleSet(std::string file, bool isCondor, double lumi)
    AnaSamples::SampleSet        ss("sampleSets.cfg", isCondor, AnaSamples::luminosity);
    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples)
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    std::map<std::string, std::vector<AnaSamples::FileSummary>> fileMap;
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
    std::set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    return vvf;
}

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

double SF_13TeV(double top_pt){

    return exp(0.0615-0.0005*top_pt);

}

bool filterEvents(NTupleReader& tr)
{
    const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const double& met = tr.getVar<double>("met");

    return jetsLVec.size() >= 4 && jetsLVec[3].Pt() > 30;// && met > 250;
}

int main(int argc, char* argv[])
{

    std::string jetVecLabel           = "jetsLVec";

    int opt;
    int option_index = 0;

    static struct option long_options[] = {
        {"condor",              no_argument, 0, 'c'},
        {"TTbar weight",        no_argument, 0, 't'},
        {"no event weighting",  no_argument, 0, 'd'},
        {"run stealth version", no_argument, 0, 's'},
        {"dataSets",      required_argument, 0, 'D'},
        {"numFiles",      required_argument, 0, 'N'},
        {"startFile",     required_argument, 0, 'M'},
        {"numEvts",       required_argument, 0, 'E'},
        {"output",        required_argument, 0, 'O'}
    };

    bool runOnCondor = false, enableTTbar = false, doWgt = true, runStealth = false;
    int nFiles = -1, startFile = 0, nEvts = -1;
    std::string dataSets = "Signal_T2tt_mStop850_mLSP100", filename = "example.root";

    while((opt = getopt_long(argc, argv, "ctdsD:N:M:E:O:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'c':
            runOnCondor = true;
            std::cout << "Configured for condor compatibility." << std::endl;
            break;

        case 't':
            enableTTbar = true;
            std::cout << "Enabled TTbar event weighting." << std::endl;
            break;

        case 'd':
            doWgt = false;
            std::cout << "No Event weighting." << std::endl;
            break;

        case 's':
            runStealth = true;
            std::cout << "Running stealth verison" << std::endl;
            break;

        case 'D':
            dataSets = optarg;
            std::cout << "Running over the " << dataSets << " data sets." << std::endl;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            std::cout << "Running over " << nFiles << " files." << std::endl;
            break;

        case 'M':
            startFile = int(atoi(optarg));
            std::cout << "Starting on file #" << startFile << std::endl;
            break;

        case 'E':
            nEvts = int(atoi(optarg));
            std::cout << "Events: " << nEvts << std::endl;
            break;

        case 'O':
            filename = optarg;
            std::cout << "Filename: " << filename << std::endl;
        }
    }

    // printing option for jet study
    bool verbose = false;

    std::string sampleloc = AnaSamples::fileDir;
    //if running on condor override all optional settings
    if(runOnCondor)
    {
        std::cout << "Run on condor is true." << std::endl;
        char thistFile[128];
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
        std::cout << "Filename modified for use with condor: " << filename << std::endl;
        sampleloc = "condor";
    }
    else
    {
        std::cout << "Run on condor is false." << std::endl;
    }

    TH1::AddDirectory(false);

    bool savefile = true;
    if(filename == "-"){
        savefile = false;
        std::cout << "Histogram file will not be saved." << std::endl;
    }
    
    const auto& myFS = getFS(dataSets, runOnCondor);

    int events = 0, pevents = 0;

    HistoContainer<NTupleReader> hists0Lep("Lep0"), hists1Lep("Lep1"), histsTTbar("ttbar"), histsTTbarNob("ttbarNob"), histsTTbarLep("ttbarLep"), histsQCD("QCD"), histsQCDb("QCDb"), histsPhoton("photon"), histsDilepton("dilepton");

    TRandom* trand = new TRandom3();

    try
    {
        for(auto& fs : myFS)
        {
            std::cout << "Tag: " << fs.tag << std::endl;

            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);

            std::cout << "File: " << fs.filePath << std::endl;
            std::cout << "Tree: " << fs.treePath << std::endl;

            BaselineVessel myBLV("", "");
            BaselineVessel blvZinv("Zinv");
            BaselineVessel blvNoVeto("NoVeto");
            BaselineVessel blvPFLeptonCleaned("PFLeptonCleaned");
            BaselineVessel blvDRLeptonCleaned("DRLeptonCleaned");
            BaselineVessel blvDRPhotonCleaned("DRPhotonCleaned");
            CleanedJets cleanedJets;
            plotterFunctions::GeneratePhotonEfficiency generatePhotonEfficiency;
            plotterFunctions::Gamma gamma;
            plotterFunctions::GenerateWeight weights;
            plotterFunctions::LepInfo lepInfo;
            plotterFunctions::BasicLepton basicLepton;
            plotterFunctions::TriggerInfo triggerInfo(false, false);
            BTagCorrector bTagCorrector("allINone_bTagEff.root", "", false, fs.tag);
            TTbarCorrector ttbarCorrector(false, "");
            ISRCorrector ISRcorrector("allINone_ISRJets.root","","");
            Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");

            NTupleReader tr(t);            
            tr.registerFunction(gamma);
            tr.registerFunction(basicLepton);
            tr.registerFunction(generatePhotonEfficiency);
            tr.registerFunction(cleanedJets);
            tr.registerFunction(myBLV);
            tr.registerFunction(lepInfo);
            tr.registerFunction(weights);
            tr.registerFunction(blvZinv);
            tr.registerFunction(blvNoVeto);
            tr.registerFunction(blvPFLeptonCleaned);
            tr.registerFunction(blvDRLeptonCleaned);
            tr.registerFunction(blvDRPhotonCleaned);
            tr.registerFunction(ISRcorrector);

            double fileWgt = fs.getWeight();

            const int printInterval = 1000;
            int printNumber = 0;
            //////////////////////////////
            // --- Start Event Loop --- //
            //////////////////////////////
            while(tr.getNextEvent())
            {
                events++;

                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                ////////////////////////////////
                // --- Jet Cleaning Study --- //
                ////////////////////////////////
                std::map<std::string, std::vector<TLorentzVector> > objMap;
                
                const auto& cutElecVec               = tr.getVec<TLorentzVector>("cutElecVec");
                const auto& cutMuVec                 = tr.getVec<TLorentzVector>("cutMuVec");
                const auto& jetsLVec                 = tr.getVec<TLorentzVector>("jetsLVec_pt20eta24");
                const auto& jetsLVec_PFLeptonCleaned = tr.getVec<TLorentzVector>("prodJetsNoLep_jetsLVec_pt20eta24");
                const auto& jetsLVec_DRLeptonCleaned = tr.getVec<TLorentzVector>("jetsLVec_drLeptonCleaned_pt20eta24");
                const auto& n_jets                   = tr.getVar<int>("cntNJetsPt20Eta24");

                std::vector<std::string> names = {"cutElecVec", "cutMuVec", "jetsLVec", "jetsLVec_DRLeptonCleaned", "jetsLVec_PFLeptonCleaned"};
                objMap["cutElecVec"] = cutElecVec;
                objMap["cutMuVec"] = cutMuVec;
                objMap["jetsLVec"] = jetsLVec;
                objMap["jetsLVec_DRLeptonCleaned"] = jetsLVec_DRLeptonCleaned;
                objMap["jetsLVec_PFLeptonCleaned"] = jetsLVec_PFLeptonCleaned;

                int min_jets = 10;

                if (n_jets > min_jets)
                {
                    pevents++;
                    printf("- n_jets (cntNJetsPt20Eta24) = %d ", n_jets);
                    if (verbose) printf("\n");
                    for (const auto& name : names)
                    {
                        if (verbose)
                        {
                            printf("- Collection: %s; %d objects\n", name.c_str(), objMap[name].size());
                            int i = 0;
                            for (const auto& obj : objMap[name])
                            {
                                printf("LVec %d: (pt, eta, phi, E) = (%f, %f, %f, %f)\n", i, obj.Pt(), obj.Eta(), obj.Phi(), obj.E());
                                i++;
                            }
                        }
                        else
                        {
                            printf("| n_%s = %d ", name.c_str(), objMap[name].size());
                        }
                    }
                    if (!verbose) printf("\n");
                }
            }
            /////////////////////////////
            // --- Stop Event Loop --- //
            /////////////////////////////
        }
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const TTException e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const SATException e)
    {
        std::cout << e << std::endl;
        return 0;
    }

    std::cout << "Processed " << events << " events. " << pevents << " passed selection." << std::endl;

    if(savefile)
    {
        std::cout << "Saving root file..." << std::endl;

        TFile *f = new TFile(filename.c_str(),"RECREATE");
        if(f->IsZombie())
        {
            std::cout << "Cannot create " << filename << std::endl;
            throw "File is zombie";
        }

        hists0Lep.save(f);
        hists1Lep.save(f);
        histsQCD.save(f);
        histsQCDb.save(f);
        histsTTbar.save(f);
        histsTTbarNob.save(f);
        histsTTbarLep.save(f);
        histsPhoton.save(f);
        histsDilepton.save(f);

        f->Write();
        f->Close();
    }
}
