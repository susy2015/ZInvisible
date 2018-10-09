#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "../../SusyAnaTools/Tools/SATException.h"

#include "../../TopTaggerTools/Tools/include/HistoContainer.h"

#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "BTagCorrector.h"
#include "TTbarCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"
#include "customize.h"

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


std::vector<TLorentzVector> GetHadTopLVec(const std::vector<TLorentzVector>& genDecayLVec, const std::vector<int>& genDecayPdgIdVec, const std::vector<int>& genDecayIdxVec, const std::vector<int>& genDecayMomIdxVec)
{
//    std::cout << "Called GetHadTopLVec" << std::endl;
    std::vector<TLorentzVector> tLVec;
//    std::cout << "genDecayPdgIDs: ";
    for(unsigned it=0; it<genDecayLVec.size(); it++)
    {
        int pdgId = genDecayPdgIdVec.at(it);
//        std::cout << " " << pdgId;
        if(abs(pdgId)==6) //pdg(6), top quark
        {

//            std::cout << "Found a genTop" << std::endl;
            for(unsigned ig=0; ig<genDecayLVec.size(); ig++)
            {
                if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ) //Is this a daughter particle of the genTop?
                {
//                    std::cout << "Found a genTop decay product" << std::endl;
                    int pdgId = genDecayPdgIdVec.at(ig);
                    if(abs(pdgId)==24) //pdg(24), W boson //Let's look at the W boson that is coming off of the top
                    {
//                        std::cout << "Found a genW from a genTop" << std::endl;
                        int flag = 0;
                        for(unsigned iq=0; iq<genDecayLVec.size(); iq++)
                        {
                            if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) //Is this a daughter particle of the W?
                            {
                                int pdgid = genDecayPdgIdVec.at(iq);
                                if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++; //pdg(11), electron; pdg(13), muon; pdg(15), tau
                            }
                        }
                        if(!flag) tLVec.push_back(genDecayLVec.at(it)); //If the W didn't have a lepton daughter product then let's include in the list of hadronic tops.
                    }
                }
            }//dau. loop
        }//top cond
//        std::cout << std::endl;
    }//genloop
    return tLVec;
}

std::vector<TLorentzVector> GetLepTopLVec(const std::vector<TLorentzVector>& genDecayLVec, const std::vector<int>& genDecayPdgIdVec, const std::vector<int>& genDecayIdxVec, const std::vector<int>& genDecayMomIdxVec)
{
//    std::cout << "Called GetLepTopLVec" << std::endl;
    std::vector<TLorentzVector> tLVec;
    for(unsigned it=0; it<genDecayLVec.size(); it++)
    {
        int pdgId = genDecayPdgIdVec.at(it);
        if(abs(pdgId)==6) //pdg(6), top quark
        {
            for(unsigned ig=0; ig<genDecayLVec.size(); ig++)
            {
                if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ) //Is this a daughter particle of the genTop?
                {
                    int pdgId = genDecayPdgIdVec.at(ig);
                    if(abs(pdgId)==24) //pdg(24), W boson //Let's look at the W boson that is coming off of the top
                    {
                        int flag = 0;
                        for(unsigned iq=0; iq<genDecayLVec.size(); iq++)
                        {
                            if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) //Is this a daughter particle of the W?
                            {
                                int pdgid = genDecayPdgIdVec.at(iq);
                                if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++; //pdg(11), electron; pdg(13), muon; pdg(15), tau
                            }
                        }
                        if(flag) tLVec.push_back(genDecayLVec.at(it)); //If the W has a lepton daughter product then let's include it in the list of leptonic tops.
                    }
                }
            }//dau. loop
        }//top cond
    }//genloop
    return tLVec;
}

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

float SF_13TeV(float top_pt){

    return exp(0.0615-0.0005*top_pt);

}

bool filterEvents(NTupleReader& tr)
{
    const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const float& met = tr.getVar<float>("met");

    return true; //Let's do a little more work for a better count. jetsLVec.size() >= 4 && jetsLVec[3].Pt() > 30;// && met > 250;
}

int main(int argc, char* argv[])
{

    std::string jetVecLabel           = "jetsLVec";

    int opt;
    int option_index = 0;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"TTbar weight",       no_argument, 0, 't'},
        {"Stored weight",      no_argument, 0, 's'},
        {"Supress Noise Filter", no_argument, 0, 'f'},
        {"Override Stored weight", no_argument, 0, 'r'},
        {"No TTbar or Btag corrections", no_argument, 0, 'a'},
        {"Use prodjetNolep branches",no_argument, 0, 'z'},
        {"no event weighting", no_argument, 0, 'd'},
        {"jecUnc up",          no_argument, 0, 'u'},
        {"jecUnc down",        no_argument, 0, 'b'}, 
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
        {"numEvts",      required_argument, 0, 'T'},
        {"output",       required_argument, 0, 'O'}
    };

    bool runOnCondor = false, enableTTbar = false, doWgt = true, enableStored = false, overrideStored = false, noCorr = false, useAltBranch = false, noNoiseFilter = false, onlyTTbarAllHad = false, onlyTTbarSingleLep = false, onlyTTbarDiLep = false;
    int nFiles = -1, startFile = 0, nEvts = -1, tF = -1;
    std::string dataSets = "Signal_T2tt_mStop850_mLSP100", filename = "example.root";

    int JECSys = 0;

    while((opt = getopt_long(argc, argv, "ctsfrazdubD:N:M:E:T:O:", long_options, &option_index)) != -1)
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

        case 's':
            enableStored = true;
            std::cout << "Enabled stored event weighting." << std::endl;
            break;

        case 'f':
            noNoiseFilter = true;
            std::cout << "Do not use noise filter." << std::endl;
            break;

        case 'r':
            overrideStored = true;
            enableStored = false;
            std::cout << "Use of stored event weight is disabled." << std::endl;
            break;

        case 'a':
            noCorr = true;
            std::cout << "Disabled Btag and TTbar corrections." << std::endl;
            break;

       case 'z':
            useAltBranch = true;
            std::cout << "Using Alternate Branches." << std::endl;
            break;

        case 'd':
            doWgt = false;
            std::cout << "No Event weighting." << std::endl;
            break;

        case 'u':
            JECSys = 1;
            break;

        case 'b':
            JECSys = -1;
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

        case 'T':
            tF = int(atoi(optarg));
            if(tF == 0){
                onlyTTbarAllHad = true;
                std::cout << "Include only TTbarAllHad = True" << std::endl;
            }
            if(tF == 1){
                onlyTTbarSingleLep = true;
                std::cout << "Include only TTbarSingleLep = True" << std::endl;
            }
            if(tF == 2){
                onlyTTbarDiLep = true;
                std::cout << "Include only TTbarDiLep = True" << std::endl;
            }
            break;

        case 'O':
            filename = optarg;
            std::cout << "Filename: " << filename << std::endl;

        }
    }

    if(JECSys == 1) std::cout << "JEC uncertainty up." << std::endl;
    if(JECSys == -1) std::cout << "JEC uncertainty down." << std::endl;

    std::string sampleloc = AnaSamples::fileDir;
    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
        std::cout << "Filename modified for use with condor: " << filename << std::endl;
        sampleloc = "condor";
    }

    TH1::AddDirectory(false);

    bool savefile = true;
    if(filename == "-"){
        savefile = false;
        std::cout << "Histogram file will not be saved." << std::endl;
    }

    std::cout << "Sample location: " << sampleloc << std::endl;

    AnaSamples::SampleSet        ss("sampleSets.txt", runOnCondor, AnaSamples::luminosity);
    AnaSamples::SampleCollection sc("sampleCollections.txt", ss);

    //if(dataSets.find("Data") != std::string::npos){
    //   std::cout << "This looks like a data n-tuple. No weighting will be applied." << std::endl;
    //   doWgt = false;
    //}

    if(dataSets.find("TT") != std::string::npos){
       std::cout << "This looks like a TTbar sample. Applying TTbar weighting" << std::endl;
       enableTTbar = true;
    }

    if(dataSets.find("Pt15to7000") != std::string::npos && !overrideStored){
       std::cout << "This looks like an HT spanning QCD sample. Applying stored weight" << std::endl;
       enableStored = true;
    }

    std::cout << "Dataset: " << dataSets << std::endl;

    int events = 0, eventsPhoton = 0, eventsQCD = 0, eventsLowHTQCD = 0, eventsTaggedLowHTQCD = 0, events1L = 0, events2L = 0, eventsTTbar1l = 0, eventsTTbar2l = 0, eventsTTbarNol = 0, fevents = 0;

    HistoContainer<NTupleReader> hists0Lep("Lep0"), hists1Lep("Lep1"), histsTTbar("ttbar"), histsTTbarLep("ttbarLep"), histsQCD("QCD"), histsLowHTQCD("lowHTQCD"), histsPhoton("photon"), histsDilepton("dilepton"), histsSimpleSemiLept("simpleSemiLep"), histsTTbar1l("histsTTbar1l"), histsTTbar2l("histsTTbar2l"), histsTTbar1lnoMET("histsTTbar1lnoMET"), histsTTbar2lnoMET("histsTTbar2lnoMET"), histsTTbarNol("histsTTbarNol");

    TRandom* trand = new TRandom3();

    int qcN = 0, qcNL = 0, qcNLJ = 0, qcNLJH = 0;
    int tcN = 0, tcNL = 0, tcNLB = 0, tcNLBJ = 0, tcNLBJM = 0;
    int lcN = 0, lcNL = 0, lcNLJ = 0, lcNLJH = 0;

    try
    {
        //for(auto& fs : sc[dataSets])
        auto& fs = ss[dataSets];
        {
            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);

            std::cout << "File: " << fs.filePath << std::endl;
            std::cout << "Tree: " << fs.treePath << std::endl;
            //std::cout << "sigma*lumi: " << fs.getWeight() << std::endl;

            //BaselineVessel myBLV(*static_cast<NTupleReader*>(nullptr), "TopTag", "");
            //plotterFunctions::SystematicPrep sysPrep;
            plotterFunctions::PrepareTopCRSelection prepTopCR(JECSys);
            plotterFunctions::PrepareTopVars prepareTopVars("TopTagger.cfg",JECSys);
            plotterFunctions::TriggerInfo triggerInfo(false, false);

            //BTagCorrector bTagCorrector("allINone_bTagEff.root", "", false);
            TTbarCorrector ttbarCorrector(false, "");
            ISRCorrector ISRcorrector("allINone_ISRJets.root","","");
            Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");

            NTupleReader tr(t);

            tr.setConvertFloatingPointVectors(true, true, true);
            tr.setConvertFloatingPointScalars(true, true, true);

            //Add some 'alias's
            //tr.addAlias("recoJetsCSVv2","recoJetsBtag_0");
            //tr.addVectorAlias<double,float>("recoJetsBtag_0","recoJetsCSVv2");
            //tr.setVectorAlias<double,float>("recoJetsBtag_0", "recoJetsCSVv2");

            if(useAltBranch){
                std::cout << "Using prodjetNolep branches" << std::endl;
                tr.setPrefix("prodJetsNoLep_");
            }

            tr.registerFunction(filterEvents);
            tr.registerFunction(prepTopCR);
            //tr.registerFunction(sysPrep);
            tr.registerFunction(prepareTopVars);
            tr.registerFunction(triggerInfo);
            //std::cout << "We are going to register bTagCorrector" << std::endl;
            //tr.registerFunction(bTagCorrector);
            tr.registerFunction(ttbarCorrector);
            tr.registerFunction(ISRcorrector);
            tr.registerFunction(pileup);
            //std::cout << "Finished registering functions" << std::endl;
            float fileWgt = fs.getWeight();

            const int printInterval = 1000;
            int printNumber = 0;
            //std::cout << "Start looping over events" << std::endl;
            while(tr.getNextEvent())
            {
                //std::cout << "Entered loop" << std::endl;
                events++;

                //int run = tr.getVar<int>("run");
                //int lumi = tr.getVar<int>("lumi");
                //long event = tr.getVar<long>("event",true);
 
                //std::cout << "Let's print run:lumi:event" << std::endl;
                //std::cout << run << "::" << lumi << "::" << event << std::endl;
                //std::cout << "Let's print the event numbers" << std::endl;

//                tr.printTupleMembers();
//                return 0;

                //std::cout << "Let's print the event numbers" << std::endl;

                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                //std::cout << "We're done printing the event number" << std::endl;

                const float& met    = tr.getVar<float>("met");
                const float& metphi = tr.getVar<float>("metphi");

                //const float& highestDisc = tr.getVar<float>("highestDisc");
                //std::cout << "Highest Discriminator: " << highestDisc << std::endl;

                TLorentzVector MET;
                MET.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

                //std::cout << "Let's calculate the selections" << std::endl;

                //Filter by MC event type

                const std::vector<TLorentzVector>& genDecayLVec = (tr.checkBranch("genDecayLVec") ? tr.getVec<TLorentzVector>("genDecayLVec") : std::vector<TLorentzVector>{});
//                std::cout << "Size of genDecayLVec: " << genDecayLVec.size() << std::endl; 
                const std::vector<int>& genDecayPdgIdVec = (tr.checkBranch("genDecayPdgIdVec") ? tr.getVec<int>("genDecayPdgIdVec") : std::vector<int>{});
                const std::vector<int>& genDecayIdxVec = (tr.checkBranch("genDecayIdxVec") ? tr.getVec<int>("genDecayIdxVec") : std::vector<int>{});
                const std::vector<int>& genDecayMomIdxVec = (tr.checkBranch("genDecayMomIdxVec") ? tr.getVec<int>("genDecayMomIdxVec") : std::vector<int>{});
                
                std::vector<TLorentzVector> hadTop = GetHadTopLVec(genDecayLVec,genDecayPdgIdVec,genDecayIdxVec,genDecayMomIdxVec);                
                std::vector<TLorentzVector> lepTop = GetLepTopLVec(genDecayLVec,genDecayPdgIdVec,genDecayIdxVec,genDecayMomIdxVec);                
                
                //if(lepTop.size()+hadTop.size() > 0) std::cout << "Total tops: " << hadTop.size()+lepTop.size() << ", Had Tops: " << hadTop.size() << ", Lep Tops: " << lepTop.size() << std::endl;
                bool bTTbarAllHad = lepTop.size() == 0 && hadTop.size() == 2;
                bool bTTbarSingleLep = lepTop.size() == 1 && hadTop.size() == 1;
                bool bTTbarDiLep = lepTop.size() == 2 && hadTop.size() == 0;

                if((onlyTTbarAllHad && !bTTbarAllHad) || (onlyTTbarSingleLep && !bTTbarSingleLep) || (onlyTTbarDiLep && !bTTbarDiLep)) continue; // We only want to include the specific ttbar sample
                fevents++;

                const bool&   passNoiseEventFilter = (tr.getVar<bool>("passNoiseEventFilter") || noNoiseFilter);
                const bool&   passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&   passLep20            = tr.getVar<bool>("passLep20");
                const bool&   passSingleLep30      = tr.getVar<bool>("passSingleLep30");
                const bool&   passSingleMu30       = tr.getVar<bool>("passSingleMu30");
                const bool&   passLeptonVeto       = tr.getVar<bool>("passLeptVeto");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhis");
                const float& ht                   = tr.getVar<float>("HT");

                const int&    nbCSV                = tr.getVar<int>("cntCSVS");

                const bool& passMuTrigger     = tr.getVar<bool>("passMuTrigger");
                const bool& passElecTrigger   = tr.getVar<bool>("passElecTrigger");
                const bool& passMETMHTTrigger = tr.getVar<bool>("passMETMHTTrigger");
                const bool& passSearchTrigger = tr.getVar<bool>("passSearchTrigger");
                const bool& passHighHtTrigger = tr.getVar<bool>("passHighHtTrigger");
                const bool& passPhotonTrigger = tr.getVar<bool>("passPhotonTrigger");

                const bool& passfloatLep      = tr.getVar<bool>("passfloatLep");

                const std::vector<TLorentzVector>& cutMuVec = tr.getVec<TLorentzVector>("cutMuVec");
                const std::vector<float>& cutMuMTlepVec = tr.getVec<float>("cutMuMTlepVec");
                const std::vector<TLorentzVector>& cutElecVec = tr.getVec<TLorentzVector>("cutElecVec");
                const std::vector<float>& cutElecMTlepVec = tr.getVec<float>("cutElecMTlepVec");

                const TopTaggerResults *ttr_                 =  tr.getVar<TopTaggerResults*>("ttrMVA",true);

                const float isData = !tr.checkBranch("genDecayLVec");

                float eWeight = fileWgt;

                float muTrigEff = 1.0;
                const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>(jetVecLabel);

                if(!isData && doWgt)
                {
                    const float& puWF               = tr.getVar<float>("_PUweightFactor");
                    //std::cout << "Calculate btag WF" << std::endl;
                    //const float& bTagWF             = tr.getVar<float>("bTagSF_EventWeightSimple_Central");
                    const float& stored_weight      = tr.getVar<float>("stored_weight");
                    if(enableTTbar & !noCorr)
                    {
                        const float& ttbarWF            = tr.getVar<float>("TTbarWF");
                        eWeight *= ttbarWF;
                    }
                    if(enableStored)
                    {
                        //std::cout << "weight: " << stored_weight << " jet pT: " << jetsLVec[0].Pt() << std::endl;
                        eWeight *= stored_weight;
                        //std::cout << "eWeight after stored weight: " << eWeight;
                    }
                    const float& triggerWF          = tr.getVar<float>("TriggerEffMC");

                    muTrigEff = tr.getVar<float>("muTrigWgt");

                    eWeight *= puWF;

                    //if(!noCorr) eWeight *= bTagWF;
                }
                //std::cout << " eWeight passed to fillHisto: " << eWeight << std::endl;

                //check on overall event weight
                //if(eWeight > 50.0 || eWeight < 1/50.0) continue;
                

                const std::vector<float>& recoJetsBtag     = tr.getVec<float>("recoJetsBtag_0");
                //std::cout << "Do we get this far?" << std::endl;
//                const std::vector<float>& recoJetsBtag     = tr.getVec<float>("recoJetsCSVv2");

                //Find lepton (here it is assumed there is exactly 1 lepton)
                TLorentzVector lepton;
                float mTLep = 999.9;
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

                const float& HT = tr.getVar<float>("HT");

                const bool& passPhoton200 = tr.getVar<bool>("passPhoton200");

                // calculate passBLep
                bool passBLep = false;
                bool passLepTtag = false;
                for(int i = 0; i < jetsLVec.size(); i++)
                {
                    //Is this a b-tagged jet (loose wp?)?
                    if(recoJetsBtag[i] < 0.8) continue;

                    float lepTopMass = (lepton + jetsLVec[i]).M();
                    if(lepTopMass > 30 && lepTopMass < 180)
                    {
                        passLepTtag = true;
                    }

                    if(jetsLVec[i].DeltaR(lepton) < 1.5)
                    {
                        passBLep = true;
                    }
                }

                float deltaPhiLepMET = fabs(lepton.DeltaPhi(MET));

                //std::cout << "Selection Critera " << run << ":" << lumi << ":" << event << ", passNoiseFilter: " << (passNoiseEventFilter ? "TRUE" : "FALSE") << ", Jets cnt: " << cntNJetsPt30Eta24 << ", HT: " << ht << ", MET: " << met << std::endl;

                //High HT QCD control sample
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 1000)
                    )
                {
                    histsQCD.fill(tr, eWeight, trand);
                    eventsQCD++;
                }

                //Look at acceptance
                if(passNoiseEventFilter) qcN++;
                if(passNoiseEventFilter&&passLeptonVeto) qcNL++;
                if(passNoiseEventFilter&&passLeptonVeto&&(cntNJetsPt30Eta24 >= 4)) qcNLJ++;
                if(passNoiseEventFilter&&passLeptonVeto&&(cntNJetsPt30Eta24 >= 4)&&(ht>1000)) qcNLJH++;
          

                //Low HT QCD control sample (For Fake study)
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 250)
                    )
                {
                    histsLowHTQCD.fill(tr, eWeight, trand);
                    eventsLowHTQCD++;
                    if(ttr_->getTops().size() > 0) eventsTaggedLowHTQCD++;
                }

                if(passNoiseEventFilter) lcN++;
                if(passNoiseEventFilter&&passLeptonVeto) lcNL++;
                if(passNoiseEventFilter&&passLeptonVeto&&(cntNJetsPt30Eta24 >= 4)) lcNLJ++;
                if(passNoiseEventFilter&&passLeptonVeto&&(cntNJetsPt30Eta24 >= 4)&&(ht>250)) lcNLJH++;

                //Simple SemiLeptonic criteria (just MC)
                if( (!isData)
                    && passNoiseEventFilter
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 250)
                    )
                {
                    histsSimpleSemiLept.fill(tr, eWeight, trand);
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
                    eventsPhoton++;
                }

                //dilepton control sample
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && cntNJetsPt30Eta24 >= 4
                    && passfloatLep                    
                    )
                {
                    histsDilepton.fill(tr, eWeight * muTrigEff, trand);
                    events2L++;
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
                    events1L++;
                }

                //TTbar1Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbar1l.fill(tr, eWeight, trand);
                    eventsTTbar1l++;
                }

                if(passNoiseEventFilter) tcN++;
                if(passNoiseEventFilter&&passLep20) tcNL++;
                if(passNoiseEventFilter&&passLep20&&(nbCSV >= 1)) tcNLB++;
                if(passNoiseEventFilter&&passLep20&&(nbCSV >= 1)&&(cntNJetsPt30Eta24 >= 4)) tcNLBJ++;
                if(passNoiseEventFilter&&passLep20&&(nbCSV >= 1)&&(cntNJetsPt30Eta24 >= 4)&&(met > 150)) tcNLBJM++;

                //TTbar2Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbar2l.fill(tr, eWeight, trand);
                    eventsTTbar2l++;
                }

                //TTbar No Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbarNol.fill(tr, eWeight, trand);
                    eventsTTbarNol++;
                }

                //TTbar1Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    )
                {
                    histsTTbar1lnoMET.fill(tr, eWeight, trand);
                }

                //TTbar2Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    )
                {
                    histsTTbar2lnoMET.fill(tr, eWeight, trand);
                }

                //semileptonic ttbar enriched control sample Mu triggered
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && passSingleMu30
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
                    histsTTbarLep.fill(tr, eWeight * muTrigEff, trand);
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

    std::cout << "Processed " << fevents<< " out of " << events << " events. " << std::endl;
    std::cout << eventsQCD << " events passed highHT QCD selection." << std::endl;
    std::cout << eventsLowHTQCD << " events passed lowHT QCD selection." << std::endl;
    std::cout << eventsPhoton << " events passed single Photon selection." << std::endl;
    std::cout << events1L << " events passed semileptonic TTbar selection." << std::endl;
    std::cout << events2L << " events passed di-mu selection." << std::endl;

    std::cout << eventsTTbar1l << " events passed the TTbar1l selecetion." << std::endl;
    std::cout << eventsTTbar2l << " events passed the TTbar2l selecetion." << std::endl;

    std::cout << eventsTTbarNol << " events passed the TTbarNol selecetion." << std::endl;

    std::cout << eventsTaggedLowHTQCD << " tagged events passed lowHT QCD selection." << std::endl;

    std::cout << std::endl << "CR acceptance (QCD):" << std::endl;
    std::cout << qcN << " events pass noise filter" << std::endl;
    std::cout << qcNL << " events pass noise filter and lepton veto" << std::endl;
    std::cout << qcNLJ << " events pass noise filter, lepton veto, and nJets" << std::endl;
    std::cout << qcNLJH << " events pass noise filter, lepton veto, nJets, and HT > 1000" << std::endl;

    std::cout << std::endl << "CR acceptance (lowHTQCD):" << std::endl;
    std::cout << lcN << " events pass noise filter" << std::endl;
    std::cout << lcNL << " events pass noise filter and lepton veto" << std::endl;
    std::cout << lcNLJ << " events pass noise filter, lepton veto, and nJets" << std::endl;
    std::cout << lcNLJH << " events pass noise filter, lepton veto, nJets, and HT > 250" << std::endl;

    std::cout << std::endl << "CR acceptance (TTbar1L):" << std::endl;
    std::cout << tcN << " events pass noise filter" << std::endl;
    std::cout << tcNL << " events pass noise filter and lepton requirement" << std::endl;
    std::cout << tcNLB << " events pass noise filter, lepton requirement, and b requirement" << std::endl;
    std::cout << tcNLBJ << " events pass noise filter, lepton requirement, b requirement, and nJets" << std::endl;
    std::cout << tcNLBJM << " events pass noise filter, lepton requirement, b requirement, nJets, and MET > 150" << std::endl;

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
        histsLowHTQCD.save(f);
        histsSimpleSemiLept.save(f);
        histsTTbar.save(f);
        histsTTbarLep.save(f);
        histsPhoton.save(f);
        histsDilepton.save(f);
        histsTTbar1l.save(f);
        histsTTbar2l.save(f);
        histsTTbarNol.save(f);
        histsTTbar1lnoMET.save(f);
        histsTTbar2lnoMET.save(f);

        f->Write();
        f->Close();
    }
}
