#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "BTagCorrector.h"
#include "TTbarCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"
#include "customize.h"


#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"

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

int main(int argc, char* argv[])
{

    std::string jetVecLabel           = "jetsLVec";

    int opt;
    int option_index = 0;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"TTbar weight",       no_argument, 0, 't'},
        {"no event weighting", no_argument, 0, 'd'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
        {"output",       required_argument, 0, 'O'}
    };

    bool runOnCondor = false, enableTTbar = false, doWgt = true;
    int nFiles = -1, startFile = 0, nEvts = -1;
    std::string dataSets = "Signal_T2tt_mStop850_mLSP100", filename = "example.root";

    while((opt = getopt_long(argc, argv, "ctdD:N:M:E:O:", long_options, &option_index)) != -1)
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

    TFile *f;

    if(savefile){
        f = new TFile(filename.c_str(),"NEW");
        if(f->IsZombie()){
            std::cout << "Cannot create " << filename << std::endl;
            return 0;
        }

        f->cd();
    }

    TH1* hMET = new TH1D("MET","MET",100,0, 1000);
    TH1* hNJets = new TH1D("nJets","nJets",21,-0.5, 20.5);
    TH1* hNVertices = new TH1D("nVertices","nVertices",61,-0.5, 60.5);
    TH1* hTopMass = new TH1D("TopMass","TopMass", 100, 0, 300);
    TH1* hTopP = new TH1D("TopP","TopP", 100 , 0, 1000);
    TH1* hTopPt = new TH1D("TopPt","TopPt", 100, 0, 1000);
    TH1* hDiTopMass = new TH1D("DiTopMass","DiTopMass", 100, 0, 1500);



    std::cout << "Sample location: " << sampleloc << std::endl;

    AnaSamples::SampleSet        ss(sampleloc, AnaSamples::luminosity);
    AnaSamples::SampleCollection sc(ss);

    if(dataSets.find("Data") != std::string::npos){
       std::cout << "This looks like a data n-tuple. No weighting will be applied." << std::endl;
       doWgt = false;
    }

    if(dataSets.find("TT") != std::string::npos){
       std::cout << "This looks like a TTbar sample. Applying TTbar weighting" << std::endl;
       enableTTbar = true;
    }

    std::cout << "Dataset: " << dataSets << std::endl;

    int events = 0, pevents = 0;


    try
    {
        //for(auto& fs : sc[dataSets])
        auto& fs = ss[dataSets];
        {
            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);

            std::cout << "Tree: " << fs.treePath << std::endl;
            //std::cout << "sigma*lumi: " << fs.getWeight() << std::endl;

            
            BaselineVessel myBLV(*static_cast<NTupleReader*>(nullptr), "TopTag", "");
            plotterFunctions::PrepareTopVars prepareTopVars;
            plotterFunctions::TriggerInfo triggerInfo(false, false);
            
            BTagCorrector bTagCorrector("allINone_bTagEff.root", "", false);
            TTbarCorrector ttbarCorrector(false, "");
            ISRCorrector ISRcorrector("allINone_ISRJets.root","","");   
            Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");

            NTupleReader tr(t);
            tr.registerFunction(myBLV);
            tr.registerFunction(prepareTopVars);
            tr.registerFunction(triggerInfo);
            tr.registerFunction(bTagCorrector);
            tr.registerFunction(ttbarCorrector);
            tr.registerFunction(ISRcorrector);
            tr.registerFunction(pileup);

            double fileWgt = fs.getWeight();


            while(tr.getNextEvent())
            {
                events++;
              
//                tr.printTupleMembers();
//                return 0;

                if(nEvts > 0 && events > nEvts) break;
                if(events % 100 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

                const double& met    = tr.getVar<double>("met");

                const bool&   passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilterTopTag");
                const bool&   passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&   passBJets            = tr.getVar<bool>("passBJetsTopTag");
                const bool&   passnJets            = tr.getVar<bool>("passnJetsTopTag");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhisTopTag");
                const double& ht                   = tr.getVar<double>("HTTopTag");
                const int&    vtxSize              = tr.getVar<int>("vtxSize");

                double eWeight = 1;

                if(doWgt){
                    const double& puWF               = tr.getVar<double>("_PUweightFactor");
                    const double& bTagWF             = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
                    if(enableTTbar){
                        const double& ttbarWF            = tr.getVar<double>("TTbarWF");         
                        eWeight *= ttbarWF;
                    }
                    const double& triggerWF          = tr.getVar<double>("TriggerEffMC");

                    const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");          
                    const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

                    eWeight *= puWF * bTagWF * triggerWF;

                }
//                std::cout << "Event Weight: " << eWeight << std::endl;

                int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt30Eta24Arr);

                const std::vector<TLorentzVector>& vTops        = tr.getVec<TLorentzVector>("vTopsNewMVA");

//                std::cout << "Noise: " << passNoiseEventFilter << ", SingleLep: " << passSingleLep20 << ", bjets: " << passBJets << ", njets: " << passnJets << ", phis: " << passdPhis << ", ht: " << ht << ", met: " << met << std::endl;
//                continue;

                if( passNoiseEventFilter 
                 && passSingleLep20
                 && passBJets
                 && passnJets
                 && passdPhis
                 && (ht >300)
                 && (met > 250))
                {

                    pevents++;
                    //std::cout << "Trigger weight: " << triggerWF << std::endl;

                    hMET->Fill(met, eWeight);
                    hNJets->Fill(cntNJetsPt30Eta24, eWeight);
                    hNVertices->Fill(vtxSize,eWeight);


                    if(vTops.size() > 0){

                        for(int tidx = 0; tidx < vTops.size(); tidx++){
                            hTopMass->Fill(vTops[tidx].M(),eWeight);
                            hTopP->Fill(vTops[tidx].Rho(),eWeight);
                            hTopPt->Fill(vTops[tidx].Perp(),eWeight);
                        }

                        if(vTops.size() == 2){
                            TLorentzVector diTop = vTops[0] + vTops[1];
                            hDiTopMass->Fill(diTop.M(),eWeight);
                        }
                    }
                    //std::cout << "MET: " << met << ", puWF: " << puWF << ", bTagWF: " << bTagWF << ", ttbarWF: " << ttbarWF << std::endl;
                }
            }
        }
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }

    std::cout << "Processed " << events << " events. " << pevents << " passed selection." << std::endl;

    if(savefile){
        std::cout << "Saving root file..." << std::endl;

        f->cd();

        hMET->Write();
        hNJets->Write();
        hNVertices->Write();
 
        hTopMass->Write();
        hTopP->Write();
        hTopPt->Write();
        hDiTopMass->Write();

        f->Write();

        f->Close();
    }
}
