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

        const auto& gammaLVec            = tr.getVec<TLorentzVector>("PhotonTLV");     // reco photon
        //const auto& gammaLVecGen         = tr.getVec<TLorentzVector>("gammaLVecGen");  // gen photon
        //const auto& genPartonLVec        = tr.getVec<TLorentzVector>("genPartonLVec"); // gen parton 
        const auto& loosePhotonID        = tr.getVec<unsigned char>("Photon_mvaID_WP80");
        const auto& tightPhotonID        = tr.getVec<unsigned char>("Photon_mvaID_WP90");
        //const auto& genMatched           = tr.getVec<data_t>("genMatched");
        const auto& Photon_genPartFlav   = tr.getVec<char>("Photon_genPartFlav");
        const auto& Photon_genPartIdx    = tr.getVec<int>("Photon_genPartIdx");
        const auto& met                  = tr.getVar<data_t>("MET_pt");
        const auto& metphi               = tr.getVar<data_t>("MET_phi");


        // toggle debugging print statements
        bool debug = false;

        //variables to be used in the analysis code
        //float photonMet = -999.9;
        //float photonPtCut = 200.0;
        float metWithPhoton = -999.9;
        float metphiWithPhoton = -999.9;
        bool passPhotonSelection = false;
        
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
        auto* gammaLVecPassMediumID     = new std::vector<TLorentzVector>();
        auto* gammaLVecPassTightID      = new std::vector<TLorentzVector>();
        auto* metLVec                   = new TLorentzVector();
        auto* metWithPhotonLVec         = new TLorentzVector();
        
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

        //Select reco photons; only eta cuts for now
        for(int i = 0; i < gammaLVec.size(); ++i)
        {
            gammaLVecReco->push_back(gammaLVec[i]);
          // passing ECAL barrel/endcap eta cuts
          // this needs to be done prior to any other cuts (pt, gen matched, etc)
          // this cut should match passAcc which is done in StopTupleMaker/SkimsAUX/plugins/PhotonIDisoProducer.cc
          if (PhotonFunctions::passPhotonECAL(gammaLVec[i])) 
          {
            gammaLVecRecoEta->push_back(gammaLVec[i]);
          }
        }
        
        // check vector lengths: gammaLVecRecoEta should have the same length as photon ntuple values for which passAcc=true
        bool passTest1 = true;
        bool passTest2 = true;
        if (gammaLVecReco->size()    != gammaLVecRecoEta->size()) passTest1 = false;
        if (gammaLVecRecoEta->size() != loosePhotonID.size())     passTest2 = false;
        if (gammaLVecRecoEta->size() != tightPhotonID.size())     passTest2 = false;
        if (debug || !passTest2) // print debugging statements
        {
          printf("gammaLVecGen gammaLVecReco gammaLVecRecoEta loosePhotonID tightPhotonID: %d | %d == %d == %d %d --- %s, %s\n", \
            int(gammaLVecGen->size()), int(gammaLVecReco->size()), int(gammaLVecRecoEta->size()), int(loosePhotonID.size()), int(tightPhotonID.size()), \
            passTest1 ? "passTest1" : "failTest1", passTest2 ? "passTest2" : "failTest2");
        }
        if (!passTest2)
        {
          // we should probably throw an exception here
          printf(" - ERROR in include/Gamma.h: TLorentzVector gammaLVecRecoEta for reco photons does not have the same length as one or more photon ntuple vectors.\n");
          printf(" - Set debug=true in include/Gamma.h for more information.\n");
          // throw exception
          try
          {
            throw 20;
          }
          catch (int e)
          {
            std::cout << "Exception: TLorentzVector photonLVecRecoEta for reco photons does not have the same length as one or more photon ntuple vectors." << std::endl;
          }
        }

        //Select reco photons within the ECAL acceptance region and Pt > 200 GeV 
        for(int i = 0; i < gammaLVecRecoEta->size(); ++i)
        {
          // pt and eta cuts
          if (PhotonFunctions::passPhotonEtaPt((*gammaLVecRecoEta)[i])) 
          {
            gammaLVecRecoEtaPt->push_back((*gammaLVecRecoEta)[i]);
            //Select iso photons passing passAcc, passIDLoose and passIsoLoose
            if(bool(loosePhotonID[i]))  gammaLVecPassLooseID  -> push_back((*gammaLVecRecoEta)[i]);
            if(bool(tightPhotonID[i]))  gammaLVecPassTightID  -> push_back((*gammaLVecRecoEta)[i]);
            if(bool(loosePhotonID[i]))
            {
              gammaLVecRecoIso->push_back((*gammaLVecRecoEta)[i]);
              // gen match
              //if (bool(genMatched[i]))
              if (PhotonFunctions::isGenMatched_Method1((*gammaLVecRecoEta)[i], *gammaLVecGen))
              {
                gammaLVecRecoEtaPtMatched->push_back((*gammaLVecRecoEta)[i]);
              }
            }
          }
        }

        // set met LVec
        // Pt, Eta, Phi, E
        //metLVec->SetPtEtaPhiE(met, 0.0, metphi, met);
        // Pt, Eta, Phi, M
        metLVec->SetPtEtaPhiM(met, 0.0, metphi, 0.0);
        metWithPhotonLVec = metLVec;
        metWithPhoton     = metLVec->Pt();
        metphiWithPhoton  = metLVec->Phi();
        // pass photon selection and add to MET
        if (gammaLVecRecoIso->size() == 1)
        {
            // Add LVecs of MET and Photon
            *metWithPhotonLVec += (*gammaLVecRecoIso)[0];
            metWithPhoton       = metWithPhotonLVec->Pt();
            metphiWithPhoton    = metWithPhotonLVec->Phi();
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
        tr.registerDerivedVar("metWithPhoton", metWithPhoton);
        tr.registerDerivedVar("metphiWithPhoton", metphiWithPhoton);
        tr.registerDerivedVar("passPhotonSelection", passPhotonSelection);
        
        tr.registerDerivedVec("gammaLVecPassLooseID", gammaLVecPassLooseID);
        tr.registerDerivedVec("gammaLVecPassMediumID", gammaLVecPassMediumID);
        tr.registerDerivedVec("gammaLVecPassTightID", gammaLVecPassTightID);
        tr.registerDerivedVec("gammaLVecGen", gammaLVecGen);
        tr.registerDerivedVec("gammaLVecGenEta", gammaLVecGenEta);
        tr.registerDerivedVec("gammaLVecGenEtaPt", gammaLVecGenEtaPt);
        tr.registerDerivedVec("gammaLVecGenEtaPtMatched", gammaLVecGenEtaPtMatched);
        tr.registerDerivedVec("gammaLVecReco", gammaLVecReco);
        tr.registerDerivedVec("gammaLVecRecoEta", gammaLVecRecoEta);
        tr.registerDerivedVec("gammaLVecRecoEtaPt", gammaLVecRecoEtaPt);
        tr.registerDerivedVec("gammaLVecRecoEtaPtMatched", gammaLVecRecoEtaPtMatched);
        tr.registerDerivedVec("gammaLVecRecoIso", gammaLVecRecoIso);
        
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
