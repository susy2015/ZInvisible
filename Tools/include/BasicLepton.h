#ifndef BASICLEPTON_H 
#define BASICLEPTON_H 

#include "TypeDefinitions.h"
#include "PhotonTools.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"

#include "TopTagger.h"
#include "TTModule.h"
#include "TopTaggerUtilities.h"
#include "TopTaggerResults.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"

#include "TopTagger/TopTagger/interface/TopObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <vector>
#include <iostream>
#include <string>
#include <set>

namespace plotterFunctions
{
    class BasicLepton
    {
    private:
        TRandom3 *tr3;
        void basicLepton(NTupleReader& tr)
        {
            const auto& muonsLVec                           = tr.getVec<TLorentzVector>("muonsLVec");
            const auto& muonsRelIso                         = tr.getVec<data_t>("muonsRelIso");
            const auto& muonsMiniIso                        = tr.getVec<data_t>("muonsMiniIso");
            const auto& muonsCharge                         = tr.getVec<data_t>("muonsCharge");
            const std::vector<data_t>& muonspfActivity      = tr.getVec<data_t>("muonspfActivity");
            const auto& muonsFlagIDVec                      = tr.getVec<int>("muonsFlagLoose");

            const auto& elesLVec                            = tr.getVec<TLorentzVector>("elesLVec");
            const auto& elesMiniIso                         = tr.getVec<data_t>("elesMiniIso");
            const auto& elesCharge                          = tr.getVec<data_t>("elesCharge");
            const auto& elesisEB                            = tr.getVec<unsigned int>("elesisEB");
            const std::vector<data_t>& elespfActivity       = tr.getVec<data_t>("elespfActivity");
            const auto& elesFlagIDVec                       = tr.getVec<int>("elesFlagVeto");

            //muons
            auto* cutMuVec                  = new std::vector<TLorentzVector>();
            auto* cutMuVecRecoOnly          = new std::vector<TLorentzVector>();
            auto* cutMuCharge               = new std::vector<data_t>();
            auto* cutMuActivity             = new std::vector<data_t>();
            
            //electrons
            auto* cutElecVec                = new std::vector<TLorentzVector>();
            auto* cutElecVecRecoOnly        = new std::vector<TLorentzVector>();
            auto* cutElecCharge             = new std::vector<data_t>();
            auto* cutElecActivity           = new std::vector<data_t>();

            //std::vector<TLorentzVector> cutMuVecRecoOnly;
            //std::vector<TLorentzVector> cutElecVecRecoOnly;
            
            //std::vector<TLorentzVector>* Zrecopt = new std::vector<TLorentzVector>();

            //muon selections
            int cutMuSummedCharge = 0;
            int nTriggerMuons = 0;

            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], 0.0, 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr)) // emulates muons with pt but no iso requirements (should this be 0.0 or -1, compare to electrons).
                {
                    cutMuVecRecoOnly->push_back(muonsLVec[i]);
                }
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                {
                    if(AnaFunctions::passMuon( muonsLVec[i], muonsRelIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                    {
                        if(nTriggerMuons == 0 && muonsLVec[i].Pt() > 17)  nTriggerMuons++;
                        else if(muonsLVec[i].Pt() > 8)  nTriggerMuons++;
                    }
                    cutMuVec->push_back(muonsLVec[i]);
                    //std::cout<<"cutMuVec PT "<<muonsLVec[i].Pt()<<std::endl; 
                    cutMuCharge->push_back(muonsCharge[i]);
                    cutMuActivity->push_back(muonspfActivity[i]);
                    if(muonsCharge[i] > 0) cutMuSummedCharge++;
                    else                   cutMuSummedCharge--;
                }
            }
            //std::cout<<"New Muon Selection "<<(*cutMuVec).size()<<std::endl;
            //electron selection
            int cutElecSummedCharge = 0;
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                if(AnaFunctions::passElectron(elesLVec[i], 0.0, -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr)) // emulates electrons with pt but no iso requirements.
                {
                    cutElecVecRecoOnly->push_back(elesLVec[i]);
                }

                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec->push_back(elesLVec[i]);
                    cutElecCharge->push_back(elesCharge[i]);
                    cutElecActivity->push_back(elespfActivity[i]);
                    if(elesCharge[i] > 0) cutElecSummedCharge++;
                    else                  cutElecSummedCharge--;
                }
            }

            //muons
            tr.registerDerivedVec("cutMuVec",             cutMuVec);
            tr.registerDerivedVec("cutMuVecRecoOnly",     cutMuVecRecoOnly);
            tr.registerDerivedVec("cutMuActivity",        cutMuActivity);
            tr.registerDerivedVec("cutMuCharge",          cutMuCharge);
            tr.registerDerivedVar("cutMuSummedCharge",    cutMuSummedCharge);
            tr.registerDerivedVar("nTriggerMuons",        nTriggerMuons);
            
            //electrons
            tr.registerDerivedVec("cutElecVec",           cutElecVec);
            tr.registerDerivedVec("cutElecVecRecoOnly",   cutElecVecRecoOnly);
            tr.registerDerivedVec("cutElecActivity",      cutElecActivity);
            tr.registerDerivedVec("cutElecCharge",        cutElecCharge);
            tr.registerDerivedVar("cutElecSummedCharge",  cutElecSummedCharge);
            

        }

    public:
        BasicLepton()
        {
            tr3 = new TRandom3();
        }

        void operator()(NTupleReader& tr)
        {
            basicLepton(tr);
        }
    };

}
#endif
