#ifndef GAMMA_H
#define GAMMA_H

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

#include "TopTagger/TopTagger/include/TopObject.h"

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                          //
// Photon Objects in CMSSW8028_2016 ntuples                                                                                 //
//                                                                                                                          //
// Parton TLorentzVector:                  genPartonLVec       form genParticles   pt > 10                                  //
// Generated Photon TLorentzVector:        gammaLVecGen        from genParticles   pt > 10                                  //
// Accepted Photon Variable (Loose):       loosePhotonID       from photonCands    passAcc                                  // 
// Accepted Photon Variable (Medium):      mediumPhotonID      from photonCands    passAcc                                  // 
// Accepted Photon Variable (Tight):       tightPhotonID       from photonCands    passAcc                                  // 
// Reconstructed Photon TLorentzVector:    gammaLVec           from photonCands    no cuts                                  //
// Reconstructed Photon Variable:          genMatched          from photonCands    passAcc                                  // 
// Full ID Isolated Photon Variable:       fullID              from photonCands    passAcc and passID and passIso           //
// Loose ID Isolated Photon Variable:      extraLooseID        from photonCands    passAcc passIDLoose and passIsoLoose     //
//                                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace plotterFunctions
{
    class Gamma {

    private:

      void generateGamma(NTupleReader& tr) {

        const auto& gammaLVec            = tr.getVec<TLorentzVector>("gammaLVec");     // reco photon
        const auto& gammaLVecGen         = tr.getVec<TLorentzVector>("gammaLVecGen");  // gen photon
        const auto& genPartonLVec        = tr.getVec<TLorentzVector>("genPartonLVec"); // gen parton 
        const auto& looseID              = tr.getVec<unsigned int>("loosePhotonID");
        const auto& mediumID             = tr.getVec<unsigned int>("mediumPhotonID");
        const auto& tightID              = tr.getVec<unsigned int>("tightPhotonID");
        //const auto& fullID               = tr.getVec<unsigned int>("fullID"); // not in CMSSW8028_2016 right now
        const auto& extraLooseID         = tr.getVec<unsigned int>("extraLooseID");
        const auto& genMatched           = tr.getVec<data_t>("genMatched");
        const auto& sigmaIetaIeta        = tr.getVec<data_t>("sigmaIetaIeta");
        const auto& pfNeutralIsoRhoCorr  = tr.getVec<data_t>("pfNeutralIsoRhoCorr");
        const auto& pfGammaIsoRhoCorr    = tr.getVec<data_t>("pfGammaIsoRhoCorr");
        const auto& pfChargedIsoRhoCorr  = tr.getVec<data_t>("pfChargedIsoRhoCorr");
        const auto& hadTowOverEM         = tr.getVec<data_t>("hadTowOverEM");
        const auto& MT2                  = tr.getVar<data_t>("best_had_brJet_MT2");
        const auto& met                  = tr.getVar<data_t>("met");
        const auto& nJets                = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
        const auto& ht                   = tr.getVar<data_t>("HT");
        const auto& nbJets               = tr.getVar<int>("cntCSVS");
        const auto& ntops                = tr.getVar<int>("nTopCandSortedCnt");

        //variables to be used in the analysis code
        double photonPtCut = 200.0;
        double photonMet = -999.9;
        auto* gammaLVecGenPtCut  = new std::vector<TLorentzVector>(); 
        auto* gammaLVecGenAcc    = new std::vector<TLorentzVector>(); 
        auto* gammaLVecRecoAcc   = new std::vector<TLorentzVector>();
        auto* gammaLVecIsoAcc    = new std::vector<TLorentzVector>(); 
        auto* promptPhotons      = new std::vector<TLorentzVector>(); 
        auto* fakePhotons        = new std::vector<TLorentzVector>();
        auto* fragmentationQCD   = new std::vector<TLorentzVector>();
        auto* loosePhotons       = new std::vector<TLorentzVector>();
        auto* mediumPhotons      = new std::vector<TLorentzVector>();
        auto* tightPhotons       = new std::vector<TLorentzVector>();
        auto* directPhotons      = new std::vector<TLorentzVector>();
        auto* totalPhotons       = new std::vector<TLorentzVector>();

        //Pass cuts; use some variables from ntuples
        
        //Select gen photons passing pt and eta cuts 
        for(int i = 0; i < gammaLVecGen.size(); ++i) {
          if (PhotonFunctions::passPhoton_Pt(gammaLVecGen[i]))    gammaLVecGenPtCut->push_back(gammaLVecGen[i]);
          if (PhotonFunctions::passPhoton_PtEta(gammaLVecGen[i])) gammaLVecGenAcc->push_back(gammaLVecGen[i]);
        }

        //Select reco photons within the ECAL acceptance region and Pt > 200 GeV 
        for(int i = 0; i < gammaLVec.size(); ++i) {
          if (   
                 PhotonFunctions::passPhoton_PtEta(gammaLVec[i])
              && PhotonFunctions::passPhoton_ECAL(gammaLVec[i])
              && bool(genMatched[i])
              ) 
          {
            gammaLVecRecoAcc->push_back(gammaLVec[i]);
            //Select iso photons passing passAcc, passIDLoose and passIsoLoose
            if(bool(extraLooseID[i]))
            {
              gammaLVecIsoAcc->push_back(gammaLVec[i]);
            }
          }
        }
        

        photonMet = met;
        //Get TLorentz vector for Loose, Medium and Tight ID photon selection
        for(int i = 0; i < gammaLVecRecoAcc->size(); i++){
          if ((*gammaLVecRecoAcc)[i].Pt() > photonPtCut){
            totalPhotons->push_back((*gammaLVecRecoAcc)[i]);
            
            if(looseID[i]) loosePhotons->push_back((*gammaLVecRecoAcc)[i]);
            if(mediumID[i]) mediumPhotons->push_back((*gammaLVecRecoAcc)[i]);
            if(tightID[i]) tightPhotons->push_back((*gammaLVecRecoAcc)[i]);

            //add loose photon pt to ptmiss
            if(looseID[i]) photonMet += (*gammaLVecRecoAcc)[i].Pt();
          } 
        }

        //Gen-Matching Photons (Pt > photonPtCut in GeV)
        if(   tr.checkBranch("gammaLVecGen")  && &gammaLVecGen != nullptr
           && tr.checkBranch("genPartonLVec") && &genPartonLVec != nullptr)
        {
          for(int i = 0; i < gammaLVecRecoAcc->size(); i++)
          {
            if((*gammaLVecRecoAcc)[i].Pt() > photonPtCut && looseID[i])
            {
              if(PhotonFunctions::isGenMatched_Method2((*gammaLVecRecoAcc)[i],gammaLVecGen))
              {
                promptPhotons->push_back((*gammaLVecRecoAcc)[i]);
                if(PhotonFunctions::isDirectPhoton((*gammaLVecRecoAcc)[i],genPartonLVec)) directPhotons->push_back((*gammaLVecRecoAcc)[i]);
                if(PhotonFunctions::isFragmentationPhoton((*gammaLVecRecoAcc)[i],genPartonLVec)) fragmentationQCD->push_back((*gammaLVecRecoAcc)[i]);
              }
              else fakePhotons->push_back((*gammaLVecRecoAcc)[i]);
            }
          }
        }

        tr.registerDerivedVar("photonMet", photonMet);
        tr.registerDerivedVar("passNphoton",totalPhotons->size() >= 1);
        tr.registerDerivedVar("passNloose",loosePhotons->size() >= 1);
        tr.registerDerivedVar("passNmedium",mediumPhotons->size() >= 1);
        tr.registerDerivedVar("passNtight",tightPhotons->size() >= 1);
        tr.registerDerivedVar("passFakes", fakePhotons->size() >= 1);
        tr.registerDerivedVar("passPrompt", promptPhotons->size() >= 1);
        tr.registerDerivedVar("passDirect", directPhotons->size() >= 1);
        tr.registerDerivedVar("passFragmentation", fragmentationQCD->size() >= 1);
        tr.registerDerivedVec("gammaLVecGenAcc",gammaLVecGenAcc);
        tr.registerDerivedVec("gammaLVecRecoAcc",gammaLVecRecoAcc);
        tr.registerDerivedVec("gammaLVecIsoAcc",gammaLVecIsoAcc);
        tr.registerDerivedVec("cutPhotons",loosePhotons);
        tr.registerDerivedVec("totalPhotons",totalPhotons);
        tr.registerDerivedVec("promptPhotons",promptPhotons);
        tr.registerDerivedVec("fakePhotons",fakePhotons);
        tr.registerDerivedVec("fragmentationQCD",fragmentationQCD);
        tr.registerDerivedVec("directPhotons",directPhotons);
        tr.registerDerivedVar("nPhotonNoID",totalPhotons->size());
        tr.registerDerivedVar("nPhoton",loosePhotons->size());
        tr.registerDerivedVar("nFakes",fakePhotons->size());
        tr.registerDerivedVar("nPrompt", promptPhotons->size());
      }

    public:

      Gamma(){}

      ~Gamma(){}

      void operator()(NTupleReader& tr)
      {
        generateGamma(tr);
      }
    };
}

#endif
