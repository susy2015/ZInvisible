#ifndef GAMMA_H
#define GAMMA_H

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

namespace plotterFunctions
{
    class Gamma {

    private:

      void generateGamma(NTupleReader& tr) {

        const std::vector<TLorentzVector>& gammaLVec  = tr.getVec<TLorentzVector>("gammaLVec");        // reco
        const std::vector<TLorentzVector>& gammaLVecGen  = tr.getVec<TLorentzVector>("gammaLVecGen");  // gen
        const std::vector<TLorentzVector>& genPartonLVec = tr.getVec<TLorentzVector>("genPartonLVec"); // gen parton 
        const std::vector<bool>& looseID = tr.getVec<bool>("loosePhotonID");
        const std::vector<bool>& mediumID = tr.getVec<bool>("mediumPhotonID");
        const std::vector<bool>& tightID = tr.getVec<bool>("tightPhotonID");
        const std::vector<int>& extraLooseID = tr.getVec<int>("extraLooseID");
        const std::vector<double>& sigmaIetaIeta = tr.getVec<double>("sigmaIetaIeta");
        const std::vector<double>& pfNeutralIsoRhoCorr = tr.getVec<double>("pfNeutralIsoRhoCorr");
        const std::vector<double>& pfGammaIsoRhoCorr = tr.getVec<double>("pfGammaIsoRhoCorr");
        const std::vector<double>& pfChargedIsoRhoCorr = tr.getVec<double>("pfChargedIsoRhoCorr");
        const std::vector<double>& hadTowOverEM = tr.getVec<double>("hadTowOverEM");
        const double& MT2 = tr.getVar<double>("best_had_brJet_MT2");
        const double& met = tr.getVar<double>("met");
        const int& nJets = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
        const double& ht = tr.getVar<double>("HT");
        const int& nbJets = tr.getVar<int>("cntCSVS");
        const int& ntops = tr.getVar<int>("nTopCandSortedCnt");

        //variables to be used in the analysis code
        double photonPtCut = 100.0;
        double photonMet = -999.9;
        std::vector<TLorentzVector> *gammaLVecGenAcc    = new std::vector<TLorentzVector>(); 
        std::vector<TLorentzVector> *promptPhotons      = new std::vector<TLorentzVector>(); 
        std::vector<TLorentzVector> *fakePhotons        = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *fragmentationQCD   = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *loosePhotons       = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *mediumPhotons      = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *tightPhotons       = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *directPhotons      = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> *totalPhotons       = new std::vector<TLorentzVector>();
        std::vector<TLorentzVector> gammaLVecRecoAcc;
        //std::vector<TLorentzVector> *tempVec            = new std::vector<TLorentzVector>();
        //std::vector<TLorentzVector> gammaLVecRecoAcc, tempVec;

        //Pass cuts that were not applied in the ntuple
        
        //Select gen photons within the ECAL acceptance region and Pt > 100 GeV 
        for(int i = 0; i < gammaLVecGen.size(); ++i) {
          if (PhotonFunctions::passPhoton(gammaLVecGen[i])) gammaLVecGenAcc->push_back(gammaLVecGen[i]);
        }

        //Select reco photons within the ECAL acceptance region and Pt > 100 GeV 
        for(int i = 0; i < gammaLVec.size(); ++i) {
          if (PhotonFunctions::passPhoton(gammaLVec[i])) gammaLVecRecoAcc.push_back(gammaLVec[i]);
        }

        photonMet = met;
        //Get TLorentz vector for Loose, Medium and Tight ID photon selection
        for(int i = 0; i < gammaLVecRecoAcc.size(); i++){
          if (gammaLVecRecoAcc[i].Pt() > photonPtCut){
            totalPhotons->push_back(gammaLVecRecoAcc[i]);
            
            if(looseID[i]) loosePhotons->push_back(gammaLVecRecoAcc[i]);
            if(mediumID[i]) mediumPhotons->push_back(gammaLVecRecoAcc[i]);
            if(tightID[i]) tightPhotons->push_back(gammaLVecRecoAcc[i]);

            //add loose photon pt to ptmiss
            if(looseID[i]) photonMet += gammaLVecRecoAcc[i].Pt();
          } 
        }

        //Gen-Matching Photons (Pt > photonPtCut in GeV)
        if(tr.checkBranch("gammaLVecGen") && tr.checkBranch("genPartonLVec") && &gammaLVecGen != nullptr && &genPartonLVec != nullptr){
          for(int i = 0; i < gammaLVecRecoAcc.size(); i++){
            if(gammaLVecRecoAcc[i].Pt() > photonPtCut && looseID[i]){
              if(PhotonFunctions::isGenMatched_Method2(gammaLVecRecoAcc[i],gammaLVecGen)){
                promptPhotons->push_back(gammaLVecRecoAcc[i]);
                if(PhotonFunctions::isDirectPhoton(gammaLVecRecoAcc[i],genPartonLVec)) directPhotons->push_back(gammaLVecRecoAcc[i]);
                if(PhotonFunctions::isFragmentationPhoton(gammaLVecRecoAcc[i],genPartonLVec)) fragmentationQCD->push_back(gammaLVecRecoAcc[i]);
              }
              else fakePhotons->push_back(gammaLVecRecoAcc[i]);
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
