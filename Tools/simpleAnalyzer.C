#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "../../SusyAnaTools/Tools/SATException.h"
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

bool filterEvents(NTupleReader& tr)
{
    const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const double& met = tr.getVar<double>("met");

    return jetsLVec.size() >= 4 && met > 250;
}

class HistoContainer
{
private:
    std::vector<TH1*> histos;

    TH1* bookHisto(const std::string& name, int nBins, double min, double max)
    {
        TH1* hptr = new TH1D(name.c_str(), name.c_str(),nBins, min, max);
        histos.push_back(hptr);
        return hptr;
    }

public:
    TH1* hMET;
    TH1* hNJets;
    TH1* hNVertices;
    TH1* hTopMass;
    TH1* hTopP;
    TH1* hTopPt;
    TH1* hDiTopMass;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopPt, *bestTopMass, *bestTopEta;

    HistoContainer()
    {
        hMET       = bookHisto("MET",100,0, 1000);
        hNJets     = bookHisto("nJets",21,-0.5, 20.5);
        hNVertices = bookHisto("nVertices",61,-0.5, 60.5);
        hTopMass   = bookHisto("TopMass", 100, 0, 300);
        hTopP      = bookHisto("TopP", 100 , 0, 1000);
        hTopPt     = bookHisto("TopPt", 100, 0, 1000);
        hDiTopMass = bookHisto("DiTopMass", 100, 0, 1500);

        bestTopPt   = bookHisto("bestTopPt",   100,  0, 1000);
        bestTopMass = bookHisto("bestTopMass", 100,  0, 500);
        bestTopEta  = bookHisto("bestTopEta",  100, -5, 5);
        bestTopCandPt   = bookHisto("bestTopCandPt",   100,  0, 1000);
        bestTopCandMass = bookHisto("bestTopCandMass", 100,  0, 500);
        bestTopCandEta  = bookHisto("bestTopCandEta",  100, -5, 5);
    }

    void save(const std::string& filename)
    {
        TFile *f;

        f = new TFile(filename.c_str(),"RECREATE");
        if(f->IsZombie()){
            std::cout << "Cannot create " << filename << std::endl;
            throw "File is zombie";
        }

        f->cd();

        for(TH1* hist : histos) hist->Write();

        f->Write();

        f->Close();
    }
};

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

    HistoContainer hists;

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
            tr.registerFunction(filterEvents);
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

                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() % 100 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

                const double& met    = tr.getVar<double>("met");

                const bool&   passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilterTopTag");
                const bool&   passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&   passBJets            = tr.getVar<bool>("passBJetsTopTag");
                const bool&   passnJets            = tr.getVar<bool>("passnJetsTopTag");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhisTopTag");
                const double& ht                   = tr.getVar<double>("HTTopTag");
                const int&    vtxSize              = tr.getVar<int>("vtxSize");

                double eWeight = fileWgt;

                if(doWgt){
                    const double& puWF               = tr.getVar<double>("_PUweightFactor");
                    const double& bTagWF             = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
                    if(enableTTbar){
                        const double& ttbarWF            = tr.getVar<double>("TTbarWF");         
                        eWeight *= ttbarWF;
                    }
                    const double& triggerWF          = tr.getVar<double>("TriggerEffMC");

                    //const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");          
                    //const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

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

                    hists.hMET->Fill(met, eWeight);
                    hists.hNJets->Fill(cntNJetsPt30Eta24, eWeight);
                    hists.hNVertices->Fill(vtxSize,eWeight);

                    //SF plots
                    if(tr.getVar<double>("bestTopMass") > 0.0)
                    {
                        const TLorentzVector& bestCandLV = tr.getVar<TLorentzVector>("bestTopMassLV");
                        hists.bestTopPt->Fill(bestCandLV.Pt(), eWeight);
                        hists.bestTopMass->Fill(bestCandLV.M(), eWeight);
                        hists.bestTopEta->Fill(bestCandLV.Eta(), eWeight);
                        
                        if(tr.getVar<bool>("bestTopMassTopTag"))
                        {
                            hists.bestTopCandPt->Fill(bestCandLV.Pt(), eWeight);
                            hists.bestTopCandMass->Fill(bestCandLV.M(), eWeight);
                            hists.bestTopCandEta->Fill(bestCandLV.Eta(), eWeight);
                        }
                    }

                    if(vTops.size() > 0)
                    {

                        for(int tidx = 0; tidx < vTops.size(); tidx++)
                        {
                            hists.hTopMass->Fill(vTops[tidx].M(),eWeight);
                            hists.hTopP->Fill(vTops[tidx].Rho(),eWeight);
                            hists.hTopPt->Fill(vTops[tidx].Perp(),eWeight);
                        }

                        if(vTops.size() == 2)
                        {
                            TLorentzVector diTop = vTops[0] + vTops[1];
                            hists.hDiTopMass->Fill(diTop.M(),eWeight);
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
        hists.save(filename);
    }
}
