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
    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples) : ss_(samples)
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
    

    //std::cout << "Sample location: " << sampleloc << std::endl;

    //AnaSamples::SampleSet        ss("sampleSets.cfg", AnaSamples::luminosity, runOnCondor);
    //AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    //if(dataSets.find("Data") != std::string::npos){
    //   std::cout << "This looks like a data n-tuple. No weighting will be applied." << std::endl;
    //   doWgt = false;
    //}

    //if(dataSets.find("TT") != std::string::npos){
    //   std::cout << "This looks like a TTbar sample. Applying TTbar weighting" << std::endl;
    //   enableTTbar = true;
    //}

    //std::cout << "Dataset: " << dataSets << std::endl;

    int events = 0, pevents = 0;

    HistoContainer<NTupleReader> hists0Lep("Lep0"), hists1Lep("Lep1"), histsTTbar("ttbar"), histsTTbarNob("ttbarNob"), histsTTbarLep("ttbarLep"), histsQCD("QCD"), histsQCDb("QCDb"), histsPhoton("photon"), histsDilepton("dilepton");

    TRandom* trand = new TRandom3();

    try
    {
        //auto& fs = ss[dataSets];
        //for(auto& fs : sc[dataSets])
        for(auto& fs : myFS)
        {
            std::cout << "Tag: " << fs.tag << std::endl;

            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);

            std::cout << "File: " << fs.filePath << std::endl;
            std::cout << "Tree: " << fs.treePath << std::endl;
            //std::cout << "sigma*lumi: " << fs.getWeight() << std::endl;

            //BaselineVessel myBLV(*static_cast<NTupleReader*>(nullptr), "TopTag", "");
            //plotterFunctions::AliasStealthVars setUpStealth;
            //plotterFunctions::PrepareTopCRSelection prepTopCR;
            //plotterFunctions::PrepareTopVars prepareTopVars;
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
            //tr.registerFunction(pileup);
            //tr.registerFunction(filterEvents);
            //tr.registerFunction(prepTopCR);
            //tr.registerFunction(prepareTopVars);
            //tr.registerFunction(triggerInfo);
            //tr.registerFunction(bTagCorrector);
            //tr.registerFunction(ttbarCorrector);
            tr.registerFunction(ISRcorrector);

            double fileWgt = fs.getWeight();
            //std::cout << "Weight: " << fileWgt << std::endl;

            const int printInterval = 1000;
            int printNumber = 0;
            while(tr.getNextEvent())
            {
                events++;

//                tr.printTupleMembers();
//                return 0;

                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                //////////////////////////
                //  Jet Cleaning Study  //
                //////////////////////////
                std::map<std::string, std::vector<TLorentzVector> > objMap;
                
                const auto& cutElecVec               = tr.getVec<TLorentzVector>("cutElecVec");
                const auto& cutMuVec                 = tr.getVec<TLorentzVector>("cutMuVec");
                const auto& jetsLVec                 = tr.getVec<TLorentzVector>("jetsLVec");
                const auto& jetsLVec_PFLeptonCleaned = tr.getVec<TLorentzVector>("prodJetsNoLep_jetsLVec");
                const auto& jetsLVec_DRLeptonCleaned = tr.getVec<TLorentzVector>("jetsLVec_drLeptonCleaned");
                const auto& n_jets                   = tr.getVar<int>("cntNJetsPt20Eta24");

                std::vector<std::string> names = {"cutElecVec", "cutMuVec", "jetsLVec", "jetsLVec_DRLeptonCleaned", "jetsLVec_PFLeptonCleaned"};
                objMap["cutElecVec"] = cutElecVec;
                objMap["cutMuVec"] = cutMuVec;
                objMap["jetsLVec"] = jetsLVec;
                objMap["jetsLVec_DRLeptonCleaned"] = jetsLVec_DRLeptonCleaned;
                objMap["jetsLVec_PFLeptonCleaned"] = jetsLVec_PFLeptonCleaned;

                bool verbose = true;
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
                



                //////////////////////////////////////
                //  Old Version of simpleAnalyzer.C //
                //////////////////////////////////////
                
                /*

                const double& met    = tr.getVar<double>("met");
                const double& metphi = tr.getVar<double>("metphi");
                TLorentzVector MET;
                MET.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

                const bool&   passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter");
                const bool&   passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&   passSingleLep30      = tr.getVar<bool>("passSingleLep30");
                const bool&   passSingleMu40       = tr.getVar<bool>("passSingleMu40");
                const bool&   passLeptonVeto       = tr.getVar<bool>("passLeptVeto");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhis");
                const double& ht                   = tr.getVar<double>("HT");

                const int&    nbCSV                = tr.getVar<int>("cntCSVS");

                const bool& passMuTrigger     = tr.getVar<bool>("passMuTrigger");
                const bool& passElecTrigger   = tr.getVar<bool>("passElecTrigger");
                const bool& passMETMHTTrigger = tr.getVar<bool>("passMETMHTTrigger");
                const bool& passSearchTrigger = tr.getVar<bool>("passSearchTrigger");
                const bool& passHighHtTrigger = tr.getVar<bool>("passHighHtTrigger");
                const bool& passPhotonTrigger = tr.getVar<bool>("passPhotonTrigger");

                const bool& passDoubleLep      = tr.getVar<bool>("passDoubleLep");

                const std::vector<TLorentzVector>& cutMuVec = tr.getVec<TLorentzVector>("cutMuVec");
                const std::vector<double>& cutMuMTlepVec = tr.getVec<double>("cutMuMTlepVec");
                const std::vector<TLorentzVector>& cutElecVec = tr.getVec<TLorentzVector>("cutElecVec");
                const std::vector<double>& cutElecMTlepVec = tr.getVec<double>("cutElecMTlepVec");

                const double isData = !tr.checkBranch("genDecayLVec");

                double eWeight = fileWgt;

                double muTrigEff = 1.0;

                if(!isData && doWgt)
                {
                    const double& puWF               = tr.getVar<double>("_PUweightFactor");
                    const double& bTagWF             = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
                    if(enableTTbar)
                    {
                        const double& ttbarWF            = tr.getVar<double>("TTbarWF");
                        eWeight *= ttbarWF;
                    }
                    const double& triggerWF          = tr.getVar<double>("TriggerEffMC");
                
                    muTrigEff = tr.getVar<double>("muTrigWgt");
                
                    eWeight *= puWF * bTagWF;
                }

                //check on overall event weight
                //if(eWeight > 50.0 || eWeight < 1/50.0) continue;
                

                const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>(jetVecLabel);
                const std::vector<double>& recoJetsBtag     = tr.getVec<double>("recoJetsBtag_0");

                //Find lepton (here it is assumed there is exactly 1 lepton)
                TLorentzVector lepton;
                double mTLep = 999.9;
                for(unsigned int i = 0; i < cutMuVec.size(); ++i)
                {
                    if(cutMuVec[i].Pt() > 20)
                    {
                        lepton = cutMuVec[i];
                        mTLep = cutMuMTlepVec[i];
                        break;
                    }
                }
                for(unsigned int i = 0; i < cutElecVec.size(); ++i)
                {
                    if(cutElecVec[i].Pt() > 20)
                    {
                        lepton = cutElecVec[i];
                        mTLep = cutElecMTlepVec[i];
                        break;
                    }
                }

                tr.registerDerivedVar("lepton", lepton);
                tr.registerDerivedVar("mTLep", mTLep);

                //                                    minAbsEta, maxAbsEta, minPt, maxPt
                const AnaConsts::AccRec pt45Eta24Arr = {-1,         2.4,      45,   -1  };

                int cntNJetsPt45Eta24 = AnaFunctions::countJets(jetsLVec,            pt45Eta24Arr);
                int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

                const double& HT = tr.getVar<double>("HT");

                const bool& passPhoton200 = tr.getVar<bool>("passPhoton200");

                // calculate passBLep
                bool passBLep = false;
                bool passLepTtag = false;
                for(int i = 0; i < jetsLVec.size(); i++)
                {
                    //Is this a b-tagged jet (loose wp?)?
                    if(recoJetsBtag[i] < 0.8) continue;

                    double lepTopMass = (lepton + jetsLVec[i]).M();
                    if(lepTopMass > 30 && lepTopMass < 180)
                    {
                        passLepTtag = true;
                    }

                    if(jetsLVec[i].DeltaR(lepton) < 1.5)
                    {
                        passBLep = true;
                    }
                }

                double deltaPhiLepMET = fabs(lepton.DeltaPhi(MET));

                //High HT QCD control sample
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 1000)
                    )
                {
                    histsQCD.fill(tr, eWeight, trand);
                }

                //High HT QCD control sample
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 1000)
                    )
                {
                    histsQCDb.fill(tr, eWeight, trand);
                }

                //photon control sample
                if( (!isData || passPhotonTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && passPhoton200
                    && ht > 400
                    )
                {
                    histsPhoton.fill(tr, eWeight, trand);
                }

                //dilepton control sample
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && cntNJetsPt30Eta24 >= 4
                    && passDoubleLep                    
                    )
                {
                    histsDilepton.fill(tr, eWeight * muTrigEff, trand);
                }

                //semileptonic ttbar enriched control sample MET triggered
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passSingleLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && passBLep
                    && passLepTtag
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 250)
                    && (met > 250)
                    )
                {
                    histsTTbar.fill(tr, eWeight, trand);
                }

                //semileptonic ttbar enriched control sample MET triggered
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passSingleLep20
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 250)
                    && (met > 250)
                    )
                {
                    histsTTbarNob.fill(tr, eWeight, trand);
                }

                //semileptonic ttbar enriched control sample Mu triggered
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && passSingleMu40
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && passBLep
                    && passLepTtag
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 200)
                    && (met > 50)
                    )
                {
                    //trigger weight really matters for single mu trigger
                    
                    histsTTbarLep.fill(tr, eWeight*muTrigEff, trand);
                }

                //Stealth Event Selection - 0 Lepton
                if( !isData  //lets not accidently unblind the stealth SR
                    && passNoiseEventFilter 
                    && passLeptonVeto
                    && cntNJetsPt45Eta24 >= 6 
                    && (ht > 500)
                    && nbCSV >= 1                //Atleast 1 medium B-Jet
                    )
                {
                    hists0Lep.fill(tr, eWeight, trand);
                }

                //Stealth Event Selection - 1 Lepton
                if( !isData  //lets not accidently unblind the stealth SR
                    && passNoiseEventFilter 
                    && passSingleLep30
                    && cntNJetsPt30Eta24 >= 6
                    && nbCSV >= 1               //Atleast 1 medium B-Jet
                    && passLepTtag
                    )
                {
                    hists1Lep.fill(tr, eWeight, trand);
                }

                */
            }
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