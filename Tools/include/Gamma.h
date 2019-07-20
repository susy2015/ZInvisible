#ifndef GAMMA_H
#define GAMMA_H

#include "TypeDefinitions.h"
#include "ZInvisible/Tools/PhotonTools.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"
#include "ZInvisible/Tools/ScaleFactors.h"
#include "ZInvisible/Tools/ScaleFactorsttBar.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TTModule.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
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
// Parton TLorentzVector:                  GenPartonTLV        from genParticles   pt > 10                                  //
// Generated Photon TLorentzVector:        GenPhotonTLV        from genParticles   pt > 10                                  //
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
        std::vector<TLorentzVector> GenPartTLV;
        std::vector<int> GenPart_pdgId;
        std::vector<int> GenPart_status;
        std::vector<int> GenPart_statusFlags;
        std::vector<int> Photon_genPartIdx;
        std::vector<unsigned char> Photon_genPartFlav;
        // get gen variables if they exist (MC only)
        bool isData = ! tr.checkBranch("GenPart_pt");
        if (! isData)
        {
            GenPartTLV            = tr.getVec<TLorentzVector>("GenPartTLV"); // gen particles
            GenPart_pdgId         = tr.getVec<int>("GenPart_pdgId");
            GenPart_status        = tr.getVec<int>("GenPart_status");
            GenPart_statusFlags   = tr.getVec<int>("GenPart_statusFlags");
            Photon_genPartIdx     = tr.getVec<int>("Photon_genPartIdx");
            Photon_genPartFlav    = tr.getVec<unsigned char>("Photon_genPartFlav");
        }
        //const auto& Photon_cutBased        = tr.getVec<int>("Photon_cutBased");
        const auto& gammaLVec              = tr.getVec<TLorentzVector>("PhotonTLV");  // reco photon
        const auto& Photon_jetIdx          = tr.getVec<int>("Photon_jetIdx");
        const auto& Photon_Stop0l          = tr.getVec<unsigned char>("Photon_Stop0l");
        const auto& met                    = tr.getVar<data_t>("MET_pt");
        const auto& metphi                 = tr.getVar<data_t>("MET_phi");
        
        // the scale factors only exist in MC, not in Data
        bool usePhotonSF = tr.checkBranch("Photon_LooseSF");
        std::vector<data_t> Photon_LooseSF;
        if (usePhotonSF) { Photon_LooseSF = tr.getVec<data_t>("Photon_LooseSF"); }   


        // toggle debugging print statements
        bool debug = false;

        //variables to be used in the analysis code
        float metWithPhoton = -999.9;
        float metphiWithPhoton = -999.9;
        float cutPhotonPt = -999.9;
        float cutPhotonEta = -999.9;
        float photonSF = 1.0;
        bool passPhotonSelection = false;
        
        // if you use new, you need to register it or destroy it yourself to clear memory
        auto* GenPartonTLV              = new std::vector<TLorentzVector>(); 
        auto* GenPhotonTLV              = new std::vector<TLorentzVector>(); 
        auto* GenPhotonTLVEta           = new std::vector<TLorentzVector>(); 
        auto* GenPhotonTLVEtaPt         = new std::vector<TLorentzVector>(); 
        auto* GenPhotonTLVEtaPtMatched  = new std::vector<TLorentzVector>(); 
        auto* RecoPhotonTLV             = new std::vector<TLorentzVector>();
        auto* RecoPhotonTLVEta          = new std::vector<TLorentzVector>();
        auto* RecoPhotonTLVEtaPt        = new std::vector<TLorentzVector>();
        auto* RecoPhotonTLVEtaPtMatched = new std::vector<TLorentzVector>();
        auto* RecoPhotonTLVIso          = new std::vector<TLorentzVector>(); 
        auto* gammaLVecPassLooseID      = new std::vector<TLorentzVector>();
        //auto* gammaLVecPassMediumID     = new std::vector<TLorentzVector>();
        //auto* gammaLVecPassTightID      = new std::vector<TLorentzVector>();
        auto* gammaGenParticleIndexPassLooseID      = new std::vector<int>();
        auto* gammaGenParticleFlavorPassLooseID     = new std::vector<unsigned char>();
        auto* gammaJetIndexPassLooseID              = new std::vector<int>();
        auto* gammaSFPassLooseID                    = new std::vector<float>();
        //auto* gammaJetIndexPassMediumID = new std::vector<int>();
        //auto* gammaJetIndexPassTightID  = new std::vector<int>();
        auto* promptPhotons     = new std::vector<TLorentzVector>();
        auto* directPhotons     = new std::vector<TLorentzVector>();
        auto* fakePhotons       = new std::vector<TLorentzVector>();
        auto* fragmentedPhotons = new std::vector<TLorentzVector>();
        
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

        //NanoAOD Gen Particles Ref: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#GenPart
        //Particle Status Codes Ref: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
        //Particle ID Numbering Ref: http://pdg.lbl.gov/2018/reviews/rpp2018-rev-monte-carlo-numbering.pdf

        // Determine GenPhotons and GenPartons from GenPart (gen particles) 
        if (! isData)
        {
            for (int i = 0; i < GenPartTLV.size(); ++i)
            {
                int pdgId       = GenPart_pdgId[i];
                int status      = GenPart_status[i];
                int statusFlags = GenPart_statusFlags[i];
                // Particle IDs
                // quarks: +/- (1 to 6)
                // gluons: + (9 and 21)
                // outgoing particles of the hardest subprocess: status == 23
                // fromHardProcess: stautsFlags == 8
                if ( ( (abs(pdgId) > 1 && abs(pdgId) < 7) || pdgId == 9 || pdgId == 21 ) && status == 23 && statusFlags == 8)
                {
                    GenPartonTLV->push_back(GenPartTLV[i]);
                }
                // Particle IDs
                // photons: +22
                // stable: status == 1
                // isPrompt: statusFlags == 0
                if (pdgId == 22 && status == 1 && statusFlags == 0)
                {
                    GenPhotonTLV->push_back(GenPartTLV[i]);
                }
            }
            //Apply cuts to Gen Photons
            for(int i = 0; i < GenPhotonTLV->size(); ++i)
            {
                // ECAL eta cuts
                if (PhotonFunctions::passPhotonECAL(GenPhotonTLV->at(i)))
                {
                    GenPhotonTLVEta->push_back(GenPhotonTLV->at(i));
                }
                // passing pt and eta cuts
                if (PhotonFunctions::passPhotonEtaPt(GenPhotonTLV->at(i)))
                {
                    GenPhotonTLVEtaPt->push_back(GenPhotonTLV->at(i));
                    // passing ECAL barrel/endcap eta cuts and reco match
                    if (PhotonFunctions::isRecoMatched(GenPhotonTLV->at(i), gammaLVec)) 
                    {
                        GenPhotonTLVEtaPtMatched->push_back(GenPhotonTLV->at(i));
                    }
                }
            }
        }


        
        // check photon vector lengths
        bool passTest1 = (gammaLVec.size() == Photon_Stop0l.size());
        if (debug || !passTest1) // print debugging statements
        {
            printf("gammaLVec == Photon_Stop0l: %d == %d --- %s\n", int(gammaLVec.size()), int(Photon_Stop0l.size()),
            passTest1 ? "passTest1" : "failTest1");
        }
        if (!passTest1)
        {
            printf(" - ERROR in include/Gamma.h: vectors of photon variables do not have the same size.\n");
            printf(" - Set debug=true in include/Gamma.h for more information.\n");
            // throw exception
            throw 20;
        }

        //Select reco photons within the ECAL acceptance region and Pt > 200 GeV 
        for(int i = 0; i < gammaLVec.size(); ++i)
        {
            RecoPhotonTLV->push_back(gammaLVec[i]);
            if (PhotonFunctions::passPhotonECAL(gammaLVec[i])) 
            {
                RecoPhotonTLVEta->push_back(gammaLVec[i]);
                // pt and eta cuts
                if (PhotonFunctions::passPhotonEtaPt(gammaLVec[i])) 
                {
                    RecoPhotonTLVEtaPt->push_back(gammaLVec[i]);
                    // Photon ID: Photon_cutBased from Photon_cutBasedBitmap  Int_t   cut-based ID bitmap, 2^(0:loose, 1:medium, 2:tight)
                    // use parentheses so that & is before ==
                    bool passLoosePhotonID = bool(Photon_Stop0l[i]);
                    //bool passLoosePhotonID  = bool(Photon_cutBased[i] & 0x1);
                    //bool passMediumPhotonID = bool(Photon_cutBased[i] & 0x2);
                    //bool passTightPhotonID  = bool(Photon_cutBased[i] & 0x4);
                    if(passLoosePhotonID)  
                    {
                        RecoPhotonTLVIso->push_back(gammaLVec[i]);
                        gammaLVecPassLooseID->push_back(gammaLVec[i]);
                        gammaJetIndexPassLooseID->push_back(Photon_jetIdx[i]);
                        // MC Only
                        if (! isData)
                        {
                            gammaGenParticleIndexPassLooseID->push_back(Photon_genPartIdx[i]);
                            gammaGenParticleFlavorPassLooseID->push_back(Photon_genPartFlav[i]);
                            if (PhotonFunctions::isGenMatched_Method1(gammaLVec[i], *GenPhotonTLV))
                            {
                                RecoPhotonTLVEtaPtMatched->push_back(gammaLVec[i]);
                                promptPhotons->push_back(gammaLVec[i]);
                                if(PhotonFunctions::isDirectPhoton(gammaLVec[i],        *GenPartonTLV))  directPhotons->push_back(gammaLVec[i]);
                                if(PhotonFunctions::isFragmentationPhoton(gammaLVec[i], *GenPartonTLV))  fragmentedPhotons->push_back(gammaLVec[i]);
                            }
                            else
                            {
                                fakePhotons->push_back(gammaLVec[i]);
                            }
                        }
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

// - Beginning of section not used (as of October 19, 2018)        
//
//      //photonMet = met;
//      //Get TLorentz vector for Loose, Medium and Tight ID photon selection
//      for(int i = 0; i < RecoPhotonTLVEta->size(); i++){
//        if ((*RecoPhotonTLVEta)[i].Pt() > photonPtCut){
//          totalPhotons->push_back((*RecoPhotonTLVEta)[i]);
//          
//          if(loosePhotonID[i]) loosePhotons->push_back((*RecoPhotonTLVEta)[i]);
//          if(tightPhotonID[i]) tightPhotons->push_back((*RecoPhotonTLVEta)[i]);
//
//          //add loose photon pt to ptmiss
//          //if(loosePhotonID[i]) photonMet += (*RecoPhotonTLVEta)[i].Pt();
//        } 
//      }
//
//      //Gen-Matching Photons (Pt > photonPtCut in GeV)
//      if(   tr.checkBranch("GenPhotonTLV")  && &GenPhotonTLV != nullptr
//         && tr.checkBranch("GenPartonTLV") && &GenPartonTLV != nullptr)
//      {
//        for(int i = 0; i < RecoPhotonTLVEta->size(); i++)
//        {
//          if((*RecoPhotonTLVEta)[i].Pt() > photonPtCut && loosePhotonID[i])
//          {
//            if(PhotonFunctions::isGenMatched_Method1((*RecoPhotonTLVEta)[i],GenPhotonTLV))
//            {
//              promptPhotons->push_back((*RecoPhotonTLVEta)[i]);
//              if(PhotonFunctions::isDirectPhoton((*RecoPhotonTLVEta)[i],GenPartonTLV)) directPhotons->push_back((*RecoPhotonTLVEta)[i]);
//              if(PhotonFunctions::isFragmentationPhoton((*RecoPhotonTLVEta)[i],GenPartonTLV)) fragmentationQCD->push_back((*RecoPhotonTLVEta)[i]);
//            }
//            else fakePhotons->push_back((*RecoPhotonTLVEta)[i]);
//          }
//        }
//      }
//
// - End of section not used (as of October 19, 2018)        

        // Register derived variables
        tr.registerDerivedVar("cutPhotonPt", cutPhotonPt);
        tr.registerDerivedVar("cutPhotonEta", cutPhotonEta);
        tr.registerDerivedVar("metWithPhoton", metWithPhoton);
        tr.registerDerivedVar("metphiWithPhoton", metphiWithPhoton);
        tr.registerDerivedVar("photonSF", photonSF);
        tr.registerDerivedVar("passPhotonSelection", passPhotonSelection);
        tr.registerDerivedVec("GenPartonTLV", GenPartonTLV);
        tr.registerDerivedVec("GenPhotonTLV", GenPhotonTLV);
        tr.registerDerivedVec("GenPhotonTLVEta", GenPhotonTLVEta);
        tr.registerDerivedVec("GenPhotonTLVEtaPt", GenPhotonTLVEtaPt);
        tr.registerDerivedVec("GenPhotonTLVEtaPtMatched", GenPhotonTLVEtaPtMatched);
        tr.registerDerivedVec("RecoPhotonTLV", RecoPhotonTLV);
        tr.registerDerivedVec("RecoPhotonTLVEta", RecoPhotonTLVEta);
        tr.registerDerivedVec("RecoPhotonTLVEtaPt", RecoPhotonTLVEtaPt);
        tr.registerDerivedVec("RecoPhotonTLVEtaPtMatched", RecoPhotonTLVEtaPtMatched);
        tr.registerDerivedVec("RecoPhotonTLVIso", RecoPhotonTLVIso);
        tr.registerDerivedVec("gammaLVecPassLooseID", gammaLVecPassLooseID);
        //tr.registerDerivedVec("gammaLVecPassMediumID", gammaLVecPassMediumID);
        //tr.registerDerivedVec("gammaLVecPassTightID", gammaLVecPassTightID);
        tr.registerDerivedVec("gammaGenParticleIndexPassLooseID",   gammaGenParticleIndexPassLooseID);
        tr.registerDerivedVec("gammaGenParticleFlavorPassLooseID",  gammaGenParticleFlavorPassLooseID);
        tr.registerDerivedVec("gammaJetIndexPassLooseID",           gammaJetIndexPassLooseID);
        //tr.registerDerivedVec("gammaJetIndexPassMediumID", gammaJetIndexPassMediumID);
        //tr.registerDerivedVec("gammaJetIndexPassTightID", gammaJetIndexPassTightID);
        tr.registerDerivedVec("gammaSFPassLooseID", gammaSFPassLooseID);
        tr.registerDerivedVec("promptPhotons", promptPhotons);
        tr.registerDerivedVec("directPhotons", directPhotons);
        tr.registerDerivedVec("fakePhotons", fakePhotons);
        tr.registerDerivedVec("fragmentedPhotons", fragmentedPhotons);
        
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
