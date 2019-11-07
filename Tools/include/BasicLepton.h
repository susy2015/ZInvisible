#ifndef BASICLEPTON_H 
#define BASICLEPTON_H 

#include "TypeDefinitions.h"
#include "ZInvisible/Tools/PhotonTools.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "ZInvisible/Tools/ScaleFactors.h"
#include "ZInvisible/Tools/ScaleFactorsttBar.h"
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
        void basicLepton(NTupleReader& tr)
        {
            const auto& muonsLVec        = tr.getVec<TLorentzVector>("MuonTLV");
            const auto& muonsMiniIso     = tr.getVec<data_t>("Muon_miniPFRelIso_all");
            const auto& muonsCharge      = tr.getVec<int>("Muon_charge");
            const auto& muonsJetIndex    = tr.getVec<int>("Muon_jetIdx");
            const auto& muonsFlagIDVec   = tr.getVec<bool_t>("Muon_mediumId");
            const auto& elesLVec         = tr.getVec<TLorentzVector>("ElectronTLV");
            const auto& elesMiniIso      = tr.getVec<data_t>("Electron_miniPFRelIso_all");
            const auto& elesCharge       = tr.getVec<int>("Electron_charge");
            const auto& elesJetIndex     = tr.getVec<int>("Electron_jetIdx");
            const auto& elesFlagIDVec    = tr.getVec<int>("Electron_cutBasedNoIso");
            
            // the scale factors only exist in MC, not in Data
            bool useMuSF   = tr.checkBranch("Muon_MediumSF");
            bool useElecSF = tr.checkBranch("Electron_MediumSF");
            std::vector<data_t> muonsScaleFactor;
            std::vector<data_t> muonsScaleFactorError;
            std::vector<data_t> elesScaleFactor;
            std::vector<data_t> elesScaleFactorError;
            if (useMuSF)    
            {
                muonsScaleFactor        = tr.getVec<data_t>("Muon_MediumSF");    
                muonsScaleFactorError   = tr.getVec<data_t>("Muon_MediumSFErr"); 
            }   
            if (useElecSF)  
            {
                elesScaleFactor         = tr.getVec<data_t>("Electron_MediumSF"); 
                elesScaleFactorError    = tr.getVec<data_t>("Electron_MediumSFErr"); 
            }   

            //muons
            auto& cutMuVec            = tr.createDerivedVec<TLorentzVector>("cutMuVec");
            auto& cutMuVecRecoOnly    = tr.createDerivedVec<TLorentzVector>("cutMuVecRecoOnly");
            auto& cutMuSF             = tr.createDerivedVec<data_t>("cutMuSF");
            auto& cutMuSF_Up          = tr.createDerivedVec<data_t>("cutMuSF_Up");
            auto& cutMuSF_Down        = tr.createDerivedVec<data_t>("cutMuSF_Down");
            auto& cutMuCharge         = tr.createDerivedVec<int>("cutMuCharge");
            auto& cutMuJetIndex       = tr.createDerivedVec<int>("cutMuJetIndex");
            //electrons
            auto& cutElecVec          = tr.createDerivedVec<TLorentzVector>("cutElecVec");
            auto& cutElecVecRecoOnly  = tr.createDerivedVec<TLorentzVector>("cutElecVecRecoOnly");
            auto& cutElecSF           = tr.createDerivedVec<data_t>("cutElecSF");
            auto& cutElecSF_Up        = tr.createDerivedVec<data_t>("cutElecSF_Up");
            auto& cutElecSF_Down      = tr.createDerivedVec<data_t>("cutElecSF_Down");
            auto& cutElecCharge       = tr.createDerivedVec<int>("cutElecCharge");
            auto& cutElecJetIndex     = tr.createDerivedVec<int>("cutElecJetIndex");

            auto& cutMuSummedCharge     = tr.createDerivedVar<int>("cutMuSummedCharge");
            auto& cutElecSummedCharge   = tr.createDerivedVar<int>("cutElecSummedCharge");
            auto& nTriggerMuons         = tr.createDerivedVar<int>("nTriggerMuons");


            //muon selections
            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], 0.0, 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr)) // emulates muons with pt but no iso requirements (should this be 0.0 or -1, compare to electrons).
                {
                    cutMuVecRecoOnly.push_back(muonsLVec[i]);
                }
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                {
                    if(nTriggerMuons == 0 && muonsLVec[i].Pt() > 17)  nTriggerMuons++;
                    else if(muonsLVec[i].Pt() > 8)  nTriggerMuons++;
                    
                    cutMuVec.push_back(muonsLVec[i]);
                    cutMuCharge.push_back(muonsCharge[i]);
                    cutMuJetIndex.push_back(muonsJetIndex[i]);
                    if (useMuSF)
                    {
                        cutMuSF.push_back(muonsScaleFactor[i]); 
                        cutMuSF_Up.push_back(muonsScaleFactor[i] + muonsScaleFactorError[i]); 
                        cutMuSF_Down.push_back(muonsScaleFactor[i] - muonsScaleFactorError[i]); 
                    } 
                    else
                    {
                        cutMuSF.push_back(1.0); 
                        cutMuSF_Up.push_back(1.0); 
                        cutMuSF_Down.push_back(1.0); 
                    } 

                    if(muonsCharge[i] > 0) cutMuSummedCharge++;
                    else                   cutMuSummedCharge--;
                }
            }
            
            //electron selection
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                // Electron_cutBased    Int_t   cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
                // Electron_cutBasedNoIso: Removed isolation requirement from eGamma ID;  Int_t  (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
                bool passElectonID = (elesFlagIDVec[i] >= 3);
                if(AnaFunctions::passElectron(elesLVec[i], 0.0, -1, passElectonID, AnaConsts::elesMiniIsoArr)) // emulates electrons with pt but no iso requirements.
                {
                    cutElecVecRecoOnly.push_back(elesLVec[i]);
                }

                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, passElectonID, AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec.push_back(elesLVec[i]);
                    cutElecCharge.push_back(elesCharge[i]);
                    cutElecJetIndex.push_back(elesJetIndex[i]);
                    if (useElecSF)
                    {
                        cutElecSF.push_back(elesScaleFactor[i]);
                        cutElecSF_Up.push_back(elesScaleFactor[i] + elesScaleFactorError[i]);
                        cutElecSF_Down.push_back(elesScaleFactor[i] - elesScaleFactorError[i]);
                    }
                    else
                    {
                        cutElecSF.push_back(1.0);
                        cutElecSF_Up.push_back(1.0);
                        cutElecSF_Down.push_back(1.0);
                    }
                    if(elesCharge[i] > 0) cutElecSummedCharge++;
                    else                  cutElecSummedCharge--;
                }
            }

        }

    public:
        BasicLepton()
        {
        }

        void operator()(NTupleReader& tr)
        {
            basicLepton(tr);
        }
    };

}
#endif
