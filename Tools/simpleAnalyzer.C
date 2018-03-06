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

#include "Math/VectorUtil.h"
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
    std::string csName_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H((csName_ + name).c_str(), (csName_ + name).c_str(), args...);
        hptr->Sumw2();
        histos.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1* hMET;
    TH1* hNJets;
    TH1* hNBJets;
    TH1* hNVertices;
    TH1 *topPt, *topMass, *topEta;
    TH1 *topCandPt, *topCandMass, *topCandEta;
    TH1 *topPtGenMatch, *topMassGenMatch, *topEtaGenMatch;
    TH1 *topCandPtGenMatch, *topCandMassGenMatch, *topCandEtaGenMatch;
    TH1 *genTopPt, *genTopMass, *genTopEta;
    TH1 *genTopMatchPt, *genTopMatchMass, *genTopMatchEta;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopCandSumPt, *bestTopCandSumMass, *bestTopCandSumEta;
    TH1 *bestTopGenPt, *bestTopGenMass, *bestTopGenEta;
    TH1 *bestTopNotGenPt, *bestTopNotGenMass, *bestTopNotGenEta;
    TH1 *bestTopPt, *bestTopMass, *bestTopEta;
    TH1 *randomTopCandPt,   *randomTopCandMass,   *randomTopCandEta;
    TH1 *randomTopPt, *randomTopMass, *randomTopEta;
    TH2 *randomTopCandMassByPt, *randomTopMassByPt;
    TH1 *fakerateMET, *fakerateNj, *fakerateNb;
    TH1 *fakerateMET2, *fakerateNj2, *fakerateNb2, *fakerateNvert2;

    TH1 *massTemplateTop, *massTemplateNotTop;

    TH2 *topCandMassByPt, *massTemplateTopByPt, *massTemplateNotTopByPt;
    TH2 *bestTopCandSumMassByPt;
    TH2 *bestTopCandSumMassRecoMatchByPt;
    TH2 *massTemplateGen0MatchByPt;
    TH2 *massTemplateGen1MatchByPt;
    TH2 *massTemplateGen2MatchByPt;
    TH2 *massTemplateGen3MatchByPt;
    TH2 *massTemplateGen0MatchRecoMatchByPt;
    TH2 *massTemplateGen1MatchRecoMatchByPt;
    TH2 *massTemplateGen2MatchRecoMatchByPt;
    TH2 *massTemplateGen3MatchRecoMatchByPt;

    TH1 *hdPhiMin, *hdPhiMax;
    TH1 *hdPhiMinGenMatch, *hdPhiMaxGenMatch;

    HistoContainer(const std::string& csName = "") : csName_(csName)
    {
        if(csName_.size() > 0) csName_ += "_";

        hMET       = bookHisto<TH1D>("MET",100,0, 1000);
        hNJets     = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
        hNBJets    = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
        hNVertices = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);

        topPt   = bookHisto<TH1D>("topPt",   100,  0, 1000);
        topMass = bookHisto<TH1D>("topMass", 100,  0, 500);
        topEta  = bookHisto<TH1D>("topEta",  100, -5, 5);
        topCandPt   = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
        topCandMass = bookHisto<TH1D>("topCandMass", 100,  0, 500);
        topCandEta  = bookHisto<TH1D>("topCandEta",  100, -5, 5);

        topPtGenMatch   = bookHisto<TH1D>("topPtGenMatch",   100,  0, 1000);
        topMassGenMatch = bookHisto<TH1D>("topMassGenMatch", 100,  0, 500);
        topEtaGenMatch  = bookHisto<TH1D>("topEtaGenMatch",  100, -5, 5);
        topCandPtGenMatch   = bookHisto<TH1D>("topCandPtGenMatch",   100,  0, 1000);
        topCandMassGenMatch = bookHisto<TH1D>("topCandMassGenMatch", 100,  0, 500);
        topCandEtaGenMatch  = bookHisto<TH1D>("topCandEtaGenMatch",  100, -5, 5);

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

        randomTopPt   = bookHisto<TH1D>("randomTopPt",   100,  0, 1000);
        randomTopMass = bookHisto<TH1D>("randomTopMass", 100,  0, 500);
        randomTopEta  = bookHisto<TH1D>("randomTopEta",  100, -5, 5);
        randomTopCandPt   = bookHisto<TH1D>("randomTopCandPt",   100,  0, 1000);
        randomTopCandMass = bookHisto<TH1D>("randomTopCandMass", 100,  0, 500);
        randomTopCandEta  = bookHisto<TH1D>("randomTopCandEta",  100, -5, 5);
        randomTopMassByPt = bookHisto<TH2D>("randomTopMassByPt", 100,  0, 500, 100, 0, 1000);
        randomTopCandMassByPt = bookHisto<TH2D>("randomTopCandMassByPt", 100,  0, 500, 100, 0, 1000);

        fakerateMET = bookHisto<TH1D>("fakerateMET", 100,0, 1000);
        fakerateNj  = bookHisto<TH1D>("fakerateNj", 21,-0.5, 20.5);
        fakerateNb  = bookHisto<TH1D>("fakerateNb", 21,-0.5, 20.5);

        fakerateMET2 = bookHisto<TH1D>("fakerateMET2", 100,0, 1000);
        fakerateNj2  = bookHisto<TH1D>("fakerateNj2", 21,-0.5, 20.5);
        fakerateNb2  = bookHisto<TH1D>("fakerateNb2", 21,-0.5, 20.5);

        massTemplateTop = bookHisto<TH1D>("massTemplateTop", 100,  0, 500);
        massTemplateNotTop = bookHisto<TH1D>("massTemplateBG", 100,  0, 500);

        topCandMassByPt = bookHisto<TH2D>("topCandMassByPt", 100,  0, 500, 100, 0, 1000);
        bestTopCandSumMassByPt = bookHisto<TH2D>("bestTopCandSumMassByPt", 100,  0, 500, 100, 0, 1000);
        bestTopCandSumMassRecoMatchByPt = bookHisto<TH2D>("bestTopCandSumMassRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateTopByPt = bookHisto<TH2D>("massTemplateTopByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateNotTopByPt = bookHisto<TH2D>("massTemplateBGByPt", 100,  0, 500, 100, 0, 1000);

        massTemplateGen0MatchByPt = bookHisto<TH2D>("massTemplateGen0MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen1MatchByPt = bookHisto<TH2D>("massTemplateGen1MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen2MatchByPt = bookHisto<TH2D>("massTemplateGen2MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen3MatchByPt = bookHisto<TH2D>("massTemplateGen3MatchByPt", 100,  0, 500, 100, 0, 1000);

        massTemplateGen0MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen0MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen1MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen1MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen2MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen2MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen3MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen3MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);

        hdPhiMin = bookHisto<TH1D>("dPhiMin", 100, 0, 3.1415);
        hdPhiMax = bookHisto<TH1D>("dPhiMax", 100, 0, 3.1415);
        hdPhiMinGenMatch = bookHisto<TH1D>("dPhiMinGenMatch", 100, 0, 3.1415);
        hdPhiMaxGenMatch = bookHisto<TH1D>("dPhiMaxGenMatch", 100, 0, 3.1415);
    }

    void fill(const NTupleReader& tr, const double& eWeight, TRandom* trand)
    {
        const double& met    = tr.getVar<double>("met");
        const double& metphi = tr.getVar<double>("metphi");

        const double& ht                   = tr.getVar<double>("HTTopTag");
        const int&    vtxSize              = tr.getVar<int>("vtxSize");
        const int&    cntCSVS              = tr.getVar<int>("cntCSVSTopTag");
        const TopTaggerResults* ttr        = tr.getVar<TopTaggerResults*>("ttrMVA");

        const int& cntNJetsPt30Eta24 = tr.getVar<int>("cntNJetsPt30Eta24TopTag");

        const std::vector<TLorentzVector>& vTops        = tr.getVec<TLorentzVector>("vTopsNewMVA");

        const TLorentzVector& lepton = tr.getVar<TLorentzVector>("lepton");

        hMET->Fill(met, eWeight);
        hNJets->Fill(cntNJetsPt30Eta24, eWeight);
        hNBJets->Fill(cntCSVS, eWeight);
        hNVertices->Fill(vtxSize,eWeight);

        const std::vector<TLorentzVector>& genTops = tr.getVec<TLorentzVector>("genTops");
        const std::vector<TLorentzVector>& genTopsRecoMatch = tr.getVec<TLorentzVector>("vTopsGenMatchTriNewMVA");

        //plots for gen efficiency
        for(const TLorentzVector& genTop : genTops)
        {
            genTopPt->Fill(genTop.Pt(), eWeight);
            genTopMass->Fill(genTop.M(), eWeight);
            genTopEta->Fill(genTop.Eta(), eWeight);
        }

        for(const TLorentzVector& genTop : genTopsRecoMatch)
        {
            genTopMatchPt->Fill(genTop.Pt(), eWeight);
            genTopMatchMass->Fill(genTop.M(), eWeight);
            genTopMatchEta->Fill(genTop.Eta(), eWeight);
        }

        //fakerate histograms
        //const auto& vTopsNCandNewMVA = tr.getVec<int>("vTopsNCandNewMVA");
        //const auto& vTopsMatchNewMVA = tr.getVec<int>("vTopsMatchNewMVABool");
        //for(unsigned int i = 0; i < vTopsNCandNewMVA.size(); ++i)
        //{
        //    if(vTopsNCandNewMVA[i] == 3 && !vTopsMatchNewMVA[i])
        //    {
        //        fakerateMET->Fill(met, eWeight);
        //        fakerateNj->Fill(cntNJetsPt30Eta24, eWeight);
        //        fakerateNb->Fill(cntCSVS, eWeight);
        //        break;
        //    }
        //}

        //SF plots
        if(tr.getVar<double>("bestTopMass") > 0.0)
        {
            const TLorentzVector& bestCandLV = tr.getVar<TLorentzVector>("bestTopMassLV");
            bestTopCandPt->Fill(bestCandLV.Pt(), eWeight);
            bestTopCandMass->Fill(bestCandLV.M(), eWeight);
            bestTopCandEta->Fill(bestCandLV.Eta(), eWeight);

            if(tr.getVar<bool>("bestTopMassTopTag"))
            {
                bestTopPt->Fill(bestCandLV.Pt(), eWeight);
                bestTopMass->Fill(bestCandLV.M(), eWeight);
                bestTopEta->Fill(bestCandLV.Eta(), eWeight);
            }

            if(tr.getVar<bool>("bestTopMassGenMatch"))
            {
                bestTopGenPt->Fill(bestCandLV.Pt(), eWeight);
                bestTopGenMass->Fill(bestCandLV.M(), eWeight);
                bestTopGenEta->Fill(bestCandLV.Eta(), eWeight);
            }
            else
            {
                bestTopNotGenPt->Fill(bestCandLV.Pt(), eWeight);
                bestTopNotGenMass->Fill(bestCandLV.M(), eWeight);
                bestTopNotGenEta->Fill(bestCandLV.Eta(), eWeight);
            }
        }

        for(auto& top : ttr->getTops())
        {
            massTemplateTop->Fill(top->p().M(), eWeight);
            massTemplateTopByPt->Fill(top->p().M(), top->p().Pt(), eWeight);

            topPt->Fill(top->p().Pt(), eWeight);
            topMass->Fill(top->p().M(), eWeight);
            topEta->Fill(top->p().Eta(), eWeight);

            if(top->getBestGenTopMatch() != nullptr)
            {
                topPtGenMatch->Fill(top->p().Pt(), eWeight);
                topMassGenMatch->Fill(top->p().M(), eWeight);
                topEtaGenMatch->Fill(top->p().Eta(), eWeight);
            }
        }

        for(auto& top : ttr->getTops())
        {
            if(top->getNConstituents() == 3 && top->getBestGenTopMatch() == nullptr)
            {
                fakerateMET->Fill(met, eWeight);
                fakerateNj->Fill(cntNJetsPt30Eta24, eWeight);
                fakerateNb->Fill(cntCSVS, eWeight);
                break;
            }
        }

        for(auto& top : ttr->getTops())
        {
            if(top->getNConstituents() == 3)
            {
                fakerateMET2->Fill(met, eWeight);
                fakerateNj2->Fill(cntNJetsPt30Eta24, eWeight);
                fakerateNb2->Fill(cntCSVS, eWeight);
                break;
            }
        }

        //Find best candiate

        //Find b jets
        std::vector<const Constituent*> bjets;
        for(const Constituent& constituent : ttr->getConstituents())
        {
            if(constituent.getBTagDisc() > 0.8484)
            {
                bjets.push_back(&constituent);
            }
        }

        //met TLorentz vector
        TLorentzVector MET;
        MET.SetPtEtaPhiM(met, 0.0, metphi, 0);

        double bestSumPtVal = 99999.999;
        const TopObject* bestCand = nullptr;
        std::vector<int> randCandIndicies;
        int iCand = 0;
        for(auto& topCand : ttr->getTopCandidates())
        {
            //delta R and Nb requirements
            bool passLepCand = lepton.DeltaR(topCand.p()) > 2;
            int nBConstituents = topCand.getNBConstituents(0.8484);
            if(passLepCand && nBConstituents == 1) randCandIndicies.push_back(iCand);

            //dPhi cut tests
            double dPhiMin = 999.99;
            double dPhiMax = -999.99;
            for(const auto* bjet : bjets)
            {
                double dPhi = fabs(ROOT::Math::VectorUtil::DeltaR(lepton + bjet->p() + MET, topCand.p()));
                for(const auto* constituent : topCand.getConstituents())
                {
                    if(bjet == constituent) continue;
                }
                dPhiMin = std::min(dPhiMin, dPhi);
                dPhiMax = std::max(dPhiMax, dPhi);
            }

            hdPhiMin->Fill(dPhiMin, eWeight);
            hdPhiMax->Fill(dPhiMax, eWeight);

            if(topCand.getBestGenTopMatch() != nullptr)
            {
                hdPhiMinGenMatch->Fill(dPhiMin, eWeight);
                hdPhiMaxGenMatch->Fill(dPhiMax, eWeight);
            }

            bool recoMatch = false;
            for(auto& top : ttr->getTops())
            {
                if(top == &topCand)
                {
                    recoMatch = true;
                    break;
                }
            }

            if(passLepCand && nBConstituents == 1)
            {
                int nGenMatch = 0;
                for(const auto& genMatch : topCand.getGenTopMatches())
                {
                    nGenMatch = std::max(nGenMatch, static_cast<int>(genMatch.second.size()));
                }
                switch(nGenMatch)
                {
                case 0:
                    massTemplateGen0MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen0MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 1:
                    massTemplateGen1MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen1MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 2:
                    massTemplateGen2MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen2MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 3:
                    massTemplateGen3MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen3MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                }

                topCandPt->Fill(topCand.p().Pt(), eWeight);
                topCandMass->Fill(topCand.p().M(), eWeight);
                topCandEta->Fill(topCand.p().Eta(), eWeight);

                if(topCand.getBestGenTopMatch() != nullptr)
                {
                    topCandPtGenMatch->Fill(topCand.p().Pt(), eWeight);
                    topCandMassGenMatch->Fill(topCand.p().M(), eWeight);
                    topCandEtaGenMatch->Fill(topCand.p().Eta(), eWeight);
                }

                topCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);

                if(topCand.getDiscriminator() > std::min(0.97, 0.8 + 0.0005*topCand.p().Pt()))
                {
                    //massTemplateTop->Fill(topCand.p().M(), eWeight);
                    //massTemplateTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                }
                else
                {
                    massTemplateNotTop->Fill(topCand.p().M(), eWeight);
                    massTemplateNotTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                }
            }

            ++iCand;
        }

        if(bestCand)
        {
            bestTopCandSumPt->Fill(bestCand->p().Pt(), eWeight);
            bestTopCandSumMass->Fill(bestCand->p().M(), eWeight);
            bestTopCandSumEta->Fill(bestCand->p().Eta(), eWeight);
            bestTopCandSumMassByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);

            for(const auto& topPtr : ttr->getTops())
            {
                if(topPtr == bestCand)
                {
                    bestTopCandSumMassRecoMatchByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);
                    break;
                }
            }
        }

        if(randCandIndicies.size() > 0)
        {
            int nCand = trand->Integer(randCandIndicies.size());

            const TopObject& topCand = ttr->getTopCandidates()[randCandIndicies[nCand]];

            randomTopCandPt->Fill(topCand.p().Pt(), eWeight);
            randomTopCandMass->Fill(topCand.p().M(), eWeight);
            randomTopCandEta->Fill(topCand.p().Eta(), eWeight);
            randomTopCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);;

            for(const auto& topPtr : ttr->getTops())
            {
                if(topPtr == &topCand)
                {
                    randomTopPt->Fill(topCand.p().Pt(), eWeight);
                    randomTopMass->Fill(topCand.p().M(), eWeight);
                    randomTopEta->Fill(topCand.p().Eta(), eWeight);
                    randomTopMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                }
            }
        }

    }

    void save(TFile *f)
    {
        f->cd();

        for(TH1* hist : histos) hist->Write();
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

    HistoContainer histsQCD("QCD"), histsTTbar("ttbar");

    TRandom* trand = new TRandom3();

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
                const bool&   passLeptonVeto       = tr.getVar<bool>("passLeptVetoTopTag");
                const bool&   passBJets            = tr.getVar<bool>("passBJetsTopTag");
                const bool&   passnJets            = tr.getVar<bool>("passnJetsTopTag");
                const bool&   passdPhis            = tr.getVar<bool>("passdPhisTopTag");
                const double& ht                   = tr.getVar<double>("HTTopTag");

                const bool& passMuTrigger     = tr.getVar<bool>("passMuTrigger");
                const bool& passElecTrigger   = tr.getVar<bool>("passElecTrigger");
                const bool& passMETMHTTrigger = tr.getVar<bool>("passMETMHTTrigger");
                const bool& passSearchTrigger = tr.getVar<bool>("passSearchTrigger");
                const bool& passHighHtTrigger = tr.getVar<bool>("passHighHtTrigger");
                const bool& passPhotonTrigger = tr.getVar<bool>("passPhotonTrigger");

                const std::vector<TLorentzVector>& cutMuVec = tr.getVec<TLorentzVector>("cutMuVec");
                const std::vector<double>& cutMuMTlepVec = tr.getVec<double>("cutMuMTlepVec");
                const std::vector<TLorentzVector>& cutElecVec = tr.getVec<TLorentzVector>("cutElecVec");
                const std::vector<double>& cutElecMTlepVec = tr.getVec<double>("cutElecMTlepVec");

                const double isData = !tr.checkBranch("genDecayLVec");

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

                    eWeight *= puWF * bTagWF;
                }

                const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>(jetVecLabel);
                const std::vector<double>& recoJetsBtag     = tr.getVec<double>("recoJetsBtag_0");
                int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

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

                // calculate passBLep
                bool passBLep = false;
                for(int i = 0; i < jetsLVec.size(); i++)
                {
                    //Is this a b-tagged jet (loose wp?)?
                    if(recoJetsBtag[i] < 0.8) continue;

                    passBLep =  passBLep || jetsLVec[i].DeltaR(lepton) < 1;
                }

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

                //semileptonic ttbar enriched control sample
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passSingleLep20
                    && passBJets
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && passBLep
                    //&& mTLep < 100
                    && (ht > 250)
                    && (met > 250)
                    )
                {
                    histsTTbar.fill(tr, eWeight, trand);
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

        TFile *f = new TFile(filename.c_str(),"RECREATE");
        if(f->IsZombie())
        {
            std::cout << "Cannot create " << filename << std::endl;
            throw "File is zombie";
        }

        histsQCD.save(f);
        histsTTbar.save(f);

        f->Write();
        f->Close();
    }
}
