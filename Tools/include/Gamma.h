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

namespace plotterFunctions
{
    class Gamma {

    private:
        std::string year_;
        bool verbose = false;

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
        const auto& gammaLVec              = tr.getVec<TLorentzVector>("PhotonTLV");  // reco photon
        const auto& Photon_jetIdx          = tr.getVec<int>("Photon_jetIdx");
        const auto& Photon_Stop0l          = tr.getVec<unsigned char>("Photon_Stop0l");
        const auto& met                    = tr.getVar<data_t>("MET_pt");
        const auto& metphi                 = tr.getVar<data_t>("MET_phi");
        
        // the scale factors only exist in MC, not in Data
        bool usePhotonSF = tr.checkBranch("Photon_LooseSF");
        std::vector<data_t> Photon_LooseSF;
        if (usePhotonSF) { Photon_LooseSF = tr.getVec<data_t>("Photon_LooseSF"); }   
        
        // --- Photon ID --- //
        // 2016:      Use Photon_cutBased       : Int_t cut-based Spring16-V2p2 ID (0:fail, 1: :loose, 2:medium, 3:tight)
        // 2017,2018: Use Photon_cutBasedBitmap : Int_t cut-based ID bitmap, 2^(0:loose, 1: medium, 2:tight); should be 2017 V2
        std::vector<int>  Photon_ID;
        std::vector<bool> Photon_PassLooseID;
        std::vector<bool> Photon_PassMediumID;
        std::vector<bool> Photon_PassTightID;
        // 2016
        if (year_.compare("2016") == 0)
        {
            Photon_ID = tr.getVec<int>("Photon_cutBased");
            for (const auto& id : Photon_ID)
            {
                Photon_PassLooseID.push_back(  bool(id > 0) );
                Photon_PassMediumID.push_back( bool(id > 1) );
                Photon_PassTightID.push_back(  bool(id > 2) );
            }
        }
        // 2017, 2018
        else
        {
            Photon_ID = tr.getVec<int>("Photon_cutBasedBitmap");
            for (const auto& id : Photon_ID)
            {
                Photon_PassLooseID.push_back(  bool(id & 1) );
                Photon_PassMediumID.push_back( bool(id & 2) );
                Photon_PassTightID.push_back(  bool(id & 4) );
            }
        }
        for (int i = 0; i < Photon_ID.size(); ++i)
        {
            std::cout << "id = " << Photon_ID[i] << "; (L, M, T) = (" << Photon_PassLooseID[i] << ", " << Photon_PassMediumID[i] << ", " << Photon_PassTightID[i] << ")" << std::endl;
        }
        
        float metWithPhoton = -999.9;
        float metphiWithPhoton = -999.9;
        float cutPhotonPt = -999.9;
        float cutPhotonEta = -999.9;
        float photonSF = 1.0;
        bool passPhotonSelection            = false;
        bool passPhotonSelectionDirect      = false;
        bool passPhotonSelectionFragmented  = false;
        bool passPhotonSelectionFake        = false;
        
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
        auto* PromptPhotons     = new std::vector<TLorentzVector>();
        auto* DirectPhotons     = new std::vector<TLorentzVector>();
        auto* FragmentedPhotons = new std::vector<TLorentzVector>();
        auto* FakePhotons       = new std::vector<TLorentzVector>();
        
        // don't use new if it will not be registered or destroyed
        TLorentzVector metWithPhotonLVec;

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
                // stautsFlags is bitwise; do not apply it
                if ( ( (abs(pdgId) > 1 && abs(pdgId) < 7) || pdgId == 9 || pdgId == 21 ) && status == 23)
                {
                    //printf("Found GenParton: pdgId = %d, status = %d, statusFlags = %d\n", pdgId, status, statusFlags);
                    GenPartonTLV->push_back(GenPartTLV[i]);
                }
                // Particle IDs
                // photons: +22
                // stable: status == 1
                // stautsFlags is bitwise; do not apply it
                if (pdgId == 22 && status == 1)
                {
                    //printf("Found GenPhoton: pdgId = %d, status = %d, statusFlags = %d\n", pdgId, status, statusFlags);
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


        
        // toggle debugging print statements
        bool debug = false;
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
                        if (verbose) printf("Found LoosePhoton; ");
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
                                if (verbose) printf("Found PromptPhoton; ");
                                RecoPhotonTLVEtaPtMatched->push_back(gammaLVec[i]);
                                PromptPhotons->push_back(gammaLVec[i]);
                                if (PhotonFunctions::isFragmentationPhoton(gammaLVec[i], *GenPartonTLV))
                                {
                                    if (verbose) printf("Found FragmentedPhoton\n");
                                    FragmentedPhotons->push_back(gammaLVec[i]);
                                }
                                // direct photon if not fragmented
                                else
                                {
                                    if (verbose) printf("Found DirectPhoton\n");
                                    DirectPhotons->push_back(gammaLVec[i]);
                                }
                            }
                            // fake photon if not prompt
                            else
                            {
                                if (verbose) printf("Found FakePhoton\n");
                                FakePhotons->push_back(gammaLVec[i]);
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
            // MC Only
            if (! isData)
            {
                if      (DirectPhotons->size() == 1)     passPhotonSelectionDirect       = true;
                else if (FragmentedPhotons->size() == 1) passPhotonSelectionFragmented   = true;
                else if (FakePhotons->size() == 1)       passPhotonSelectionFake         = true;
            }                                               
        }
        
        // Register derived variables
        tr.registerDerivedVar("cutPhotonPt", cutPhotonPt);
        tr.registerDerivedVar("cutPhotonEta", cutPhotonEta);
        tr.registerDerivedVar("metWithPhoton", metWithPhoton);
        tr.registerDerivedVar("metphiWithPhoton", metphiWithPhoton);
        tr.registerDerivedVar("photonSF", photonSF);
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
        tr.registerDerivedVec("PromptPhotons", PromptPhotons);
        tr.registerDerivedVec("DirectPhotons", DirectPhotons);
        tr.registerDerivedVec("FragmentedPhotons", FragmentedPhotons);
        tr.registerDerivedVec("FakePhotons", FakePhotons);
        tr.registerDerivedVar("passPhotonSelection", passPhotonSelection);
        tr.registerDerivedVar("passPhotonSelectionDirect", passPhotonSelectionDirect);
        tr.registerDerivedVar("passPhotonSelectionFragmented", passPhotonSelectionFragmented);
        tr.registerDerivedVar("passPhotonSelectionFake", passPhotonSelectionFake);
    }

    public:

        Gamma(std::string year = "") : year_(year)
        {
        }

        ~Gamma(){}

        void operator()(NTupleReader& tr)
        {
            generateGamma(tr);
        }
    };
}

#endif
