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

#include "TopTaggerResults.h"
#include "Constituent.h"

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

    return jetsLVec.size() >= 4;//met > 250;
}

class HistoContainer
{
private:
    std::vector<TH1*> histos;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H(name.c_str(), name.c_str(), args...);
        hptr->Sumw2();
        histos.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1* hMET;
    TH1* hNJets;
    TH1* hNBJets;
    TH1* hNVertices;
    TH1* hTopMass;
    TH1* hTopP;
    TH1* hTopPt;
    TH1* hDiTopMass;
    TH1 *topPt, *topMass, *topEta;
    TH1 *topCandPt, *topCandMass, *topCandEta;
    TH1 *genTopPt, *genTopMass, *genTopEta;
    TH1 *genTopMatchPt, *genTopMatchMass, *genTopMatchEta;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopCandSumPt, *bestTopCandSumMass, *bestTopCandSumEta;
    TH1 *bestTopGenPt, *bestTopGenMass, *bestTopGenEta;
    TH1 *bestTopNotGenPt, *bestTopNotGenMass, *bestTopNotGenEta;
    TH1 *bestTopPt, *bestTopMass, *bestTopEta;
    TH1 *fakerateMET, *fakerateNj, *fakerateNb;

    TH1 *massTemplateTop, *massTemplateNotTop;

    TH1 *allSumPt, *bestSumPt, *genSumPt;

    TH2 *topCandMassByPt, *massTemplateTopByPt, *massTemplateNotTopByPt;
    TH2 *bestTopCandSumMassByPt;
    TH2 *bestTopCandSumMassRecoMatchByPt;
    TH2 *massTemplateGen0MatchByPt;
    TH2 *massTemplateGen1MatchByPt;
    TH2 *massTemplateGen2MatchByPt;
    TH2 *massTemplateGen3MatchByPt;

    HistoContainer()
    {
        hMET       = bookHisto<TH1D>("MET",100,0, 1000);
        hNJets     = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
        hNBJets    = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
        hNVertices = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);
        hTopMass   = bookHisto<TH1D>("TopMass", 100, 0, 300);
        hTopP      = bookHisto<TH1D>("TopP", 100 , 0, 1000);
        hTopPt     = bookHisto<TH1D>("TopPt", 100, 0, 1000);
        hDiTopMass = bookHisto<TH1D>("DiTopMass", 100, 0, 1500);

        topPt   = bookHisto<TH1D>("topPt",   100,  0, 1000);
        topMass = bookHisto<TH1D>("topMass", 100,  0, 500);
        topEta  = bookHisto<TH1D>("topEta",  100, -5, 5);
        topCandPt   = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
        topCandMass = bookHisto<TH1D>("topCandMass", 100,  0, 500);
        topCandEta  = bookHisto<TH1D>("topCandEta",  100, -5, 5);

        genTopPt   = bookHisto<TH1D>("genTopPt",   100,  0, 1000);
        genTopMass = bookHisto<TH1D>("genTopMass", 100,  0, 500);
        genTopEta  = bookHisto<TH1D>("genTopEta",  100, -5, 5);
        genTopMatchPt   = bookHisto<TH1D>("genTopMatchPt",   100,  0, 1000);
        genTopMatchMass = bookHisto<TH1D>("genTopMatchMass", 100,  0, 500);
        genTopMatchEta  = bookHisto<TH1D>("genTopMatchEta",  100, -5, 5);

        bestTopPt   = bookHisto<TH1D>("bestTopPt",   100,  0, 1000);
        bestTopMass = bookHisto<TH1D>("bestTopMass", 100,  0, 500);
        bestTopEta  = bookHisto<TH1D>("bestTopEta",  100, -5, 5);
        bestTopCandPt   = bookHisto<TH1D>("bestTopCandPt",   100,  0, 1000);
        bestTopCandMass = bookHisto<TH1D>("bestTopCandMass", 100,  0, 500);
        bestTopCandEta  = bookHisto<TH1D>("bestTopCandEta",  100, -5, 5);
        bestTopCandSumPt   = bookHisto<TH1D>("bestTopCandSumPt",   100,  0, 1000);
        bestTopCandSumMass = bookHisto<TH1D>("bestTopCandSumMass", 100,  0, 500);
        bestTopCandSumEta  = bookHisto<TH1D>("bestTopCandSumEta",  100, -5, 5);
        bestTopGenPt   = bookHisto<TH1D>("bestTopGenPt",   100,  0, 1000);
        bestTopGenMass = bookHisto<TH1D>("bestTopGenMass", 100,  0, 500);
        bestTopGenEta  = bookHisto<TH1D>("bestTopGenEta",  100, -5, 5);
        bestTopNotGenPt   = bookHisto<TH1D>("bestTopNotGenPt",   100,  0, 1000);
        bestTopNotGenMass = bookHisto<TH1D>("bestTopNotGenMass", 100,  0, 500);
        bestTopNotGenEta  = bookHisto<TH1D>("bestTopNotGenEta",  100, -5, 5);

        fakerateMET = bookHisto<TH1D>("fakerateMET", 100,0, 1000);
        fakerateNj  = bookHisto<TH1D>("fakerateNj", 21,-0.5, 20.5);
        fakerateNb  = bookHisto<TH1D>("fakerateNb", 21,-0.5, 20.5);

        massTemplateTop = bookHisto<TH1D>("massTemplateTop", 100,  0, 500);
        massTemplateNotTop = bookHisto<TH1D>("massTemplateBG", 100,  0, 500);

        allSumPt = bookHisto<TH1D>("allSumPt", 100,  0, 1000);
        bestSumPt = bookHisto<TH1D>("bestSumPt", 100,  0, 1000);
        genSumPt = bookHisto<TH1D>("genSumPt", 100,  0, 1000);

        topCandMassByPt = bookHisto<TH2D>("topCandMassByPt", 100,  0, 500, 100, 0, 1000);
        bestTopCandSumMassByPt = bookHisto<TH2D>("bestTopCandSumMassByPt", 100,  0, 500, 100, 0, 1000);
        bestTopCandSumMassRecoMatchByPt = bookHisto<TH2D>("bestTopCandSumMassRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateTopByPt = bookHisto<TH2D>("massTemplateTopByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateNotTopByPt = bookHisto<TH2D>("massTemplateBGByPt", 100,  0, 500, 100, 0, 1000);

        massTemplateGen0MatchByPt = bookHisto<TH2D>("massTemplateGen0MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen1MatchByPt = bookHisto<TH2D>("massTemplateGen1MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen2MatchByPt = bookHisto<TH2D>("massTemplateGen2MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen3MatchByPt = bookHisto<TH2D>("massTemplateGen3MatchByPt", 100,  0, 500, 100, 0, 1000);

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

            std::cout << "File: " << fs.filePath << std::endl;
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

                const double& met    = tr.getVar<double>("met");
                const double& metphi = tr.getVar<double>("metphi");

                const bool&   passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilterTopTag");
                const bool&   passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&   passBJets            = tr.getVar<bool>("passBJetsTopTag");
                const bool&   passnJets            = tr.getVar<bool>("passnJetsTopTag");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhisTopTag");
                const double& ht                   = tr.getVar<double>("HTTopTag");
                const int&    vtxSize              = tr.getVar<int>("vtxSize");
                const int&    cntCSVS              = tr.getVar<int>("cntCSVS");
                const TopTaggerResults* ttr        = tr.getVar<TopTaggerResults*>("ttrMVA");

                const std::vector<TLorentzVector>& cutMuVec = tr.getVec<TLorentzVector>("cutMuVec");
                const std::vector<TLorentzVector>& cutElecVec = tr.getVec<TLorentzVector>("cutElecVec");

                double eWeight = fileWgt;

                if(doWgt)
                {
                    const double& puWF               = tr.getVar<double>("_PUweightFactor");
                    const double& bTagWF             = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
                    if(enableTTbar)
                    {
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
//                 && passSingleLep20
                 && passBJets
                 && passnJets
//                 && passdPhis
                    && (ht > 1000)//300)
//                 && (met > 250)
                    )
                {

                    pevents++;
                    //std::cout << "Trigger weight: " << triggerWF << std::endl;

                    hists.hMET->Fill(met, eWeight);
                    hists.hNJets->Fill(cntNJetsPt30Eta24, eWeight);
                    hists.hNBJets->Fill(cntCSVS, eWeight);
                    hists.hNVertices->Fill(vtxSize,eWeight);

                    const std::vector<TLorentzVector>& genTops = tr.getVec<TLorentzVector>("genTops");
                    const std::vector<TLorentzVector>& genTopsRecoMatch = tr.getVec<TLorentzVector>("vTopsGenMatchTriNewMVA");

                    //plots for gen efficiency 
                    for(const TLorentzVector& genTop : genTops)
                    {
                        hists.genTopPt->Fill(genTop.Pt(), eWeight);
                        hists.genTopMass->Fill(genTop.M(), eWeight);
                        hists.genTopEta->Fill(genTop.Eta(), eWeight);
                    }

                    for(const TLorentzVector& genTop : genTopsRecoMatch)
                    {
                        hists.genTopMatchPt->Fill(genTop.Pt(), eWeight);
                        hists.genTopMatchMass->Fill(genTop.M(), eWeight);
                        hists.genTopMatchEta->Fill(genTop.Eta(), eWeight);
                    }

                    //fakerate histograms 
                    const auto& vTopsNCandNewMVA = tr.getVec<int>("vTopsNCandNewMVA");
                    const auto& vTopsMatchNewMVA = tr.getVec<int>("vTopsMatchNewMVABool");
                    for(unsigned int i = 0; i < vTopsNCandNewMVA.size(); ++i)
                    {
                        if(vTopsNCandNewMVA[i] == 3 && !vTopsMatchNewMVA[i])
                        {
                            hists.fakerateMET->Fill(met, eWeight);
                            hists.fakerateNj->Fill(cntNJetsPt30Eta24, eWeight);
                            hists.fakerateNb->Fill(cntCSVS, eWeight);
                            break;
                        }
                    }

                    //SF plots
                    if(tr.getVar<double>("bestTopMass") > 0.0)
                    {
                        const TLorentzVector& bestCandLV = tr.getVar<TLorentzVector>("bestTopMassLV");
                        hists.bestTopCandPt->Fill(bestCandLV.Pt(), eWeight);
                        hists.bestTopCandMass->Fill(bestCandLV.M(), eWeight);
                        hists.bestTopCandEta->Fill(bestCandLV.Eta(), eWeight);
                        
                        if(tr.getVar<bool>("bestTopMassTopTag"))
                        {
                            hists.bestTopPt->Fill(bestCandLV.Pt(), eWeight);
                            hists.bestTopMass->Fill(bestCandLV.M(), eWeight);
                            hists.bestTopEta->Fill(bestCandLV.Eta(), eWeight);
                        }

                        if(tr.getVar<bool>("bestTopMassGenMatch"))
                        {
                            hists.bestTopGenPt->Fill(bestCandLV.Pt(), eWeight);
                            hists.bestTopGenMass->Fill(bestCandLV.M(), eWeight);
                            hists.bestTopGenEta->Fill(bestCandLV.Eta(), eWeight);
                        }
                        else
                        {
                            hists.bestTopNotGenPt->Fill(bestCandLV.Pt(), eWeight);
                            hists.bestTopNotGenMass->Fill(bestCandLV.M(), eWeight);
                            hists.bestTopNotGenEta->Fill(bestCandLV.Eta(), eWeight);
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

                    for(auto& top : ttr->getTops())
                    {
                        hists.massTemplateTop->Fill(top->p().M(), eWeight);
                        hists.massTemplateTopByPt->Fill(top->p().M(), top->p().Pt(), eWeight);

                        hists.topPt->Fill(top->p().Pt(), eWeight);
                        hists.topMass->Fill(top->p().M(), eWeight);
                        hists.topEta->Fill(top->p().Eta(), eWeight);

                    }

                    //Find best candiate
                    
                    //Find b jets
                    std::vector<const Constituent*> bjets;
                    for(const auto& constituent : ttr->getConstituents())
                    {
                        if(constituent.getBTagDisc() > 0.8484)
                        {
                            bjets.push_back(&constituent);
                        }
                    }

                    //Find lepton (here it is assumed there is exactly 1 lepton)
                    TLorentzVector lepton;
                    for(const auto& lep : cutMuVec)
                    {
                        if(lep.Pt() > 20)
                        {
                            lepton = lep;
                            break;
                        }
                    }
                    for(const auto& lep : cutElecVec)
                    {
                        if(lep.Pt() > 20)
                        {
                            lepton = lep;
                            break;
                        }
                    }

                    //met TLorentz vector
                    TLorentzVector MET;
                    MET.SetPtEtaPhiM(met, 0.0, metphi, 0);

                    double bestSumPt = 99999.999;
                    const TopObject* bestCand = nullptr;
                    for(auto& topCand : ttr->getTopCandidates())
                    {
                        switch(topCand.getGenTopMatches().size())
                        {
                        case 0:
                            hists.massTemplateGen0MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                            break;
                        case 1:
                            hists.massTemplateGen1MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                            break;
                        case 2:
                            hists.massTemplateGen2MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                            break;
                        case 3:
                            hists.massTemplateGen3MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                            break;
                        }

                        const auto& constituents = topCand.getConstituents();
                        for(const auto& bjet : bjets)
                        {
                            //Check that the b-jet is not inside the top candidate 
                            if(std::find(constituents.begin(), constituents.end(), bjet) == constituents.end())
                            {
                                double sumPt = (bjet->p() + MET + topCand.p() + lepton).Pt();
                                hists.allSumPt->Fill(sumPt, eWeight);
                                if(topCand.getBestGenTopMatch() != nullptr) hists.genSumPt->Fill(sumPt, eWeight);
                                if(sumPt < bestSumPt)
                                {
                                    bestSumPt = sumPt;
                                    bestCand = &topCand;
                                }
                            }
                        }

                        hists.topCandPt->Fill(topCand.p().Pt(), eWeight);
                        hists.topCandMass->Fill(topCand.p().M(), eWeight);
                        hists.topCandEta->Fill(topCand.p().Eta(), eWeight);

                        hists.topCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);

                        if(topCand.getDiscriminator() > std::min(0.97, 0.8 + 0.0005*topCand.p().Pt()))
                        {
                            //hists.massTemplateTop->Fill(topCand.p().M(), eWeight);
                            //hists.massTemplateTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                        }
                        else
                        {
                            hists.massTemplateNotTop->Fill(topCand.p().M(), eWeight);
                            hists.massTemplateNotTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                        }
                    }

                    if(bestCand)
                    {
                        hists.bestTopCandSumPt->Fill(bestCand->p().Pt(), eWeight);
                        hists.bestTopCandSumMass->Fill(bestCand->p().M(), eWeight);
                        hists.bestTopCandSumEta->Fill(bestCand->p().Eta(), eWeight);
                        hists.bestTopCandSumMassByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);
                        hists.bestSumPt->Fill(bestSumPt, eWeight);

                        for(const auto& topPtr : ttr->getTops()) 
                        {
                            if(topPtr == bestCand) 
                            {
                                hists.bestTopCandSumMassRecoMatchByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);
                                break;
                            }
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
