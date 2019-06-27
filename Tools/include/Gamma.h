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
        //std::cout << "Running Gamma.h" << std::endl;

        const auto& gammaLVec              = tr.getVec<TLorentzVector>("PhotonTLV");     // reco photon
        //const auto& gammaLVecGen         = tr.getVec<TLorentzVector>("gammaLVecGen");  // gen photon
        //const auto& genPartonLVec        = tr.getVec<TLorentzVector>("genPartonLVec"); // gen parton 
        //const auto& Photon_cutBased        = tr.getVec<int>("Photon_cutBased");
        const auto& Photon_jetIdx          = tr.getVec<int>("Photon_jetIdx");
        const auto& Photon_Stop0l          = tr.getVec<unsigned char>("Photon_Stop0l");
        //const auto& loosePhotonID        = tr.getVec<unsigned char>("Photon_mvaID_WP80");
        //const auto& tightPhotonID        = tr.getVec<unsigned char>("Photon_mvaID_WP90");
        //const auto& genMatched           = tr.getVec<data_t>("genMatched");
        const auto& met                    = tr.getVar<data_t>("MET_pt");
        const auto& metphi                 = tr.getVar<data_t>("MET_phi");
        
        
        // the scale factors only exist in MC, not in Data
        //const auto& Photon_LooseSF         = tr.getVec<data_t>("Photon_LooseSF");
        bool usePhotonSF = tr.checkBranch("Photon_LooseSF");
        std::vector<data_t> Photon_LooseSF;
        if (usePhotonSF) { Photon_LooseSF = tr.getVec<data_t>("Photon_LooseSF"); }   


        // toggle debugging print statements
        bool debug = false;

        //variables to be used in the analysis code
        //float photonMet = -999.9;
        //float photonPtCut = 200.0;
        float metWithPhoton = -999.9;
        float metphiWithPhoton = -999.9;
        float cutPhotonPt = -999.9;
        float cutPhotonEta = -999.9;
        float photonSF = 1.0;
        bool passPhotonSelection = false;
        
        // if you use new, you need to register it or destroy it yourself to clear memory
        auto* gammaLVecGen              = new std::vector<TLorentzVector>(); 
        auto* gammaLVecGenEta           = new std::vector<TLorentzVector>(); 
        auto* gammaLVecGenEtaPt         = new std::vector<TLorentzVector>(); 
        auto* gammaLVecGenEtaPtMatched  = new std::vector<TLorentzVector>(); 
        auto* gammaLVecReco             = new std::vector<TLorentzVector>();
        auto* gammaLVecRecoEta          = new std::vector<TLorentzVector>();
        auto* gammaLVecRecoEtaPt        = new std::vector<TLorentzVector>();
        auto* gammaLVecRecoEtaPtMatched = new std::vector<TLorentzVector>();
        auto* gammaLVecRecoIso          = new std::vector<TLorentzVector>(); 
        auto* gammaLVecPassLooseID      = new std::vector<TLorentzVector>();
        //auto* gammaLVecPassMediumID     = new std::vector<TLorentzVector>();
        //auto* gammaLVecPassTightID      = new std::vector<TLorentzVector>();
        auto* gammaJetIndexPassLooseID  = new std::vector<int>();
        auto* gammaSFPassLooseID        = new std::vector<float>();
        //auto* gammaJetIndexPassMediumID = new std::vector<int>();
        //auto* gammaJetIndexPassTightID  = new std::vector<int>();
        
        // don't use new if it will not be registered or destroyed
        TLorentzVector metWithPhotonLVec;
        
        //auto* promptPhotons             = new std::vector<TLorentzVector>(); 
        //auto* fakePhotons               = new std::vector<TLorentzVector>();
        //auto* fragmentationQCD          = new std::vector<TLorentzVector>();
        //auto* loosePhotons              = new std::vector<TLorentzVector>();
        //auto* mediumPhotons             = new std::vector<TLorentzVector>();
        //auto* tightPhotons              = new std::vector<TLorentzVector>();
        //auto* directPhotons             = new std::vector<TLorentzVector>();
        //auto* totalPhotons              = new std::vector<TLorentzVector>();

        // Determine Gen Photons 

        //Pass cuts; use some variables from ntuples
        
        //Select gen photons
        for(int i = 0; i < gammaLVecGen->size(); ++i)
        {
          // ECAL eta cuts
          if (PhotonFunctions::passPhotonECAL(gammaLVecGen->at(i)))
          {
            gammaLVecGenEta->push_back(gammaLVecGen->at(i));
          }
          // passing pt and eta cuts
          if (PhotonFunctions::passPhotonEtaPt(gammaLVecGen->at(i)))
          {
            gammaLVecGenEtaPt->push_back(gammaLVecGen->at(i));
            // passing ECAL barrel/endcap eta cuts and reco match
            if (PhotonFunctions::isRecoMatched(gammaLVecGen->at(i), gammaLVec)) 
            {
              gammaLVecGenEtaPtMatched->push_back(gammaLVecGen->at(i));
            }
          }
        }

        
        // check photon vector lengths
        bool passTest1 = (gammaLVec.size() == Photon_Stop0l.size());
        if (debug || !passTest1) // print debugging statements
        {
          printf("gammaLVecGen | gammaLVec == Photon_Stop0l: %d | %d == %d --- %s\n",
            int(gammaLVecGen->size()), int(gammaLVec.size()), int(Photon_Stop0l.size()),
            passTest1 ? "passTest1" : "failTest1");
        }
        if (!passTest1)
        {
          // we should probably throw an exception here
          printf(" - ERROR in include/Gamma.h: vectors of photon variables do not have the same size.\n");
          printf(" - Set debug=true in include/Gamma.h for more information.\n");
          // throw exception
          throw 20;
          // try
          // {
          //   throw 20;
          // }
          // catch (int e)
          // {
          //   std::cout << "Exception: photon ntuple vectors do not have the same length." << std::endl;
          // }
        }

        //Select reco photons within the ECAL acceptance region and Pt > 200 GeV 
        for(int i = 0; i < gammaLVec.size(); ++i)
        {
          gammaLVecReco->push_back(gammaLVec[i]);
          if (PhotonFunctions::passPhotonECAL(gammaLVec[i])) 
          {
            gammaLVecRecoEta->push_back(gammaLVec[i]);
            // pt and eta cuts
            if (PhotonFunctions::passPhotonEtaPt(gammaLVec[i])) 
            {
              gammaLVecRecoEtaPt->push_back(gammaLVec[i]);
              // Photon ID: Photon_cutBased from Photon_cutBasedBitmap  Int_t   cut-based ID bitmap, 2^(0:loose, 1:medium, 2:tight)
              // use parentheses so that & is before ==
              bool passLoosePhotonID = bool(Photon_Stop0l[i]);
              //bool passLoosePhotonID  = bool(Photon_cutBased[i] & 0x1);
              //bool passMediumPhotonID = bool(Photon_cutBased[i] & 0x2);
              //bool passTightPhotonID  = bool(Photon_cutBased[i] & 0x4);
              if(passLoosePhotonID)  
              {
                  gammaLVecPassLooseID->push_back(gammaLVec[i]);
                  gammaJetIndexPassLooseID->push_back(Photon_jetIdx[i]);
                  if (usePhotonSF)
                  {
                      gammaSFPassLooseID->push_back(Photon_LooseSF[i]);
                  }
                  else
                  {
                      gammaSFPassLooseID->push_back(1.0);
                  }
              }
              //if(passMediumPhotonID) 
              //{
              //    gammaLVecPassMediumID->push_back(gammaLVec[i]);
              //    gammaJetIndexPassMediumID->push_back(Photon_jetIdx[i]);
              //}
              //if(passTightPhotonID)  
              //{
              //    gammaLVecPassTightID->push_back(gammaLVec[i]);
              //    gammaJetIndexPassTightID->push_back(Photon_jetIdx[i]);
              //}
              if(passLoosePhotonID)
              {
                gammaLVecRecoIso->push_back(gammaLVec[i]);
                // gen match
                //if (bool(genMatched[i]))
                if (PhotonFunctions::isGenMatched_Method1(gammaLVec[i], *gammaLVecGen))
                {
                  gammaLVecRecoEtaPtMatched->push_back(gammaLVec[i]);
                }
              }
            }
          }
        }

        // set default met LVec using met and metphi
        // Pt, Eta, Phi, M
        metWithPhotonLVec.SetPtEtaPhiM(met, 0.0, metphi, 0.0);
        metWithPhoton     = metWithPhotonLVec.Pt();
        metphiWithPhoton  = metWithPhotonLVec.Phi();
        // pass photon selection and add to MET
        if (gammaLVecPassLooseID->size() == 1)
        {
            cutPhotonPt  = (*gammaLVecPassLooseID)[0].Pt();
            cutPhotonEta = (*gammaLVecPassLooseID)[0].Eta();
            // Add LVecs of MET and Photon
            metWithPhotonLVec  += (*gammaLVecPassLooseID)[0];
            metWithPhoton       = metWithPhotonLVec.Pt();
            metphiWithPhoton    = metWithPhotonLVec.Phi();
            photonSF            = (*gammaSFPassLooseID)[0];
            passPhotonSelection = true;
        }

// --- Beginning of section not used (as of October 19, 2018)        
//
//        //photonMet = met;
//        //Get TLorentz vector for Loose, Medium and Tight ID photon selection
//        for(int i = 0; i < gammaLVecRecoEta->size(); i++){
//          if ((*gammaLVecRecoEta)[i].Pt() > photonPtCut){
//            totalPhotons->push_back((*gammaLVecRecoEta)[i]);
//            
//            if(loosePhotonID[i]) loosePhotons->push_back((*gammaLVecRecoEta)[i]);
//            if(tightPhotonID[i]) tightPhotons->push_back((*gammaLVecRecoEta)[i]);
//
//            //add loose photon pt to ptmiss
//            //if(loosePhotonID[i]) photonMet += (*gammaLVecRecoEta)[i].Pt();
//          } 
//        }
//
//        //Gen-Matching Photons (Pt > photonPtCut in GeV)
//        if(   tr.checkBranch("gammaLVecGen")  && &gammaLVecGen != nullptr
//           && tr.checkBranch("genPartonLVec") && &genPartonLVec != nullptr)
//        {
//          for(int i = 0; i < gammaLVecRecoEta->size(); i++)
//          {
//            if((*gammaLVecRecoEta)[i].Pt() > photonPtCut && loosePhotonID[i])
//            {
//              if(PhotonFunctions::isGenMatched_Method2((*gammaLVecRecoEta)[i],gammaLVecGen))
//              {
//                promptPhotons->push_back((*gammaLVecRecoEta)[i]);
//                if(PhotonFunctions::isDirectPhoton((*gammaLVecRecoEta)[i],genPartonLVec)) directPhotons->push_back((*gammaLVecRecoEta)[i]);
//                if(PhotonFunctions::isFragmentationPhoton((*gammaLVecRecoEta)[i],genPartonLVec)) fragmentationQCD->push_back((*gammaLVecRecoEta)[i]);
//              }
//              else fakePhotons->push_back((*gammaLVecRecoEta)[i]);
//            }
//          }
//        }
//
// --- End of section not used (as of October 19, 2018)        

        // Register derived variables
        tr.registerDerivedVar("cutPhotonPt", cutPhotonPt);
        tr.registerDerivedVar("cutPhotonEta", cutPhotonEta);
        tr.registerDerivedVar("metWithPhoton", metWithPhoton);
        tr.registerDerivedVar("metphiWithPhoton", metphiWithPhoton);
        tr.registerDerivedVar("photonSF", photonSF);
        tr.registerDerivedVar("passPhotonSelection", passPhotonSelection);
        tr.registerDerivedVec("gammaLVecGen", gammaLVecGen);
        tr.registerDerivedVec("gammaLVecGenEta", gammaLVecGenEta);
        tr.registerDerivedVec("gammaLVecGenEtaPt", gammaLVecGenEtaPt);
        tr.registerDerivedVec("gammaLVecGenEtaPtMatched", gammaLVecGenEtaPtMatched);
        tr.registerDerivedVec("gammaLVecReco", gammaLVecReco);
        tr.registerDerivedVec("gammaLVecRecoEta", gammaLVecRecoEta);
        tr.registerDerivedVec("gammaLVecRecoEtaPt", gammaLVecRecoEtaPt);
        tr.registerDerivedVec("gammaLVecRecoEtaPtMatched", gammaLVecRecoEtaPtMatched);
        tr.registerDerivedVec("gammaLVecRecoIso", gammaLVecRecoIso);
        tr.registerDerivedVec("gammaLVecPassLooseID", gammaLVecPassLooseID);
        //tr.registerDerivedVec("gammaLVecPassMediumID", gammaLVecPassMediumID);
        //tr.registerDerivedVec("gammaLVecPassTightID", gammaLVecPassTightID);
        tr.registerDerivedVec("gammaJetIndexPassLooseID", gammaJetIndexPassLooseID);
        //tr.registerDerivedVec("gammaJetIndexPassMediumID", gammaJetIndexPassMediumID);
        //tr.registerDerivedVec("gammaJetIndexPassTightID", gammaJetIndexPassTightID);
        tr.registerDerivedVec("gammaSFPassLooseID", gammaSFPassLooseID);
        
        //tr.registerDerivedVar("photonMet", photonMet);
        //tr.registerDerivedVec("cutPhotons", loosePhotons);
        //tr.registerDerivedVec("totalPhotons", totalPhotons);
        //tr.registerDerivedVec("promptPhotons", promptPhotons);
        //tr.registerDerivedVec("fakePhotons", fakePhotons);
        //tr.registerDerivedVec("fragmentationQCD", fragmentationQCD);
        //tr.registerDerivedVec("directPhotons", directPhotons);
        //tr.registerDerivedVar("nPhotonNoID", totalPhotons->size());
        //tr.registerDerivedVar("nPhoton", loosePhotons->size());
        //tr.registerDerivedVar("nFakes", fakePhotons->size());
        //tr.registerDerivedVar("nPrompt", promptPhotons->size());
        //tr.registerDerivedVar("passNphoton", totalPhotons->size() >= 1);
        //tr.registerDerivedVar("passNloose", loosePhotons->size() >= 1);
        //tr.registerDerivedVar("passNmedium", mediumPhotons->size() >= 1);
        //tr.registerDerivedVar("passNtight", tightPhotons->size() >= 1);
        //tr.registerDerivedVar("passFakes", fakePhotons->size() >= 1);
        //tr.registerDerivedVar("passPrompt", promptPhotons->size() >= 1);
        //tr.registerDerivedVar("passDirect", directPhotons->size() >= 1);
        //tr.registerDerivedVar("passFragmentation", fragmentationQCD->size() >= 1);
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
