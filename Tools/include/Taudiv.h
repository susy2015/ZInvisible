#ifndef TAUDIV_H
#define TAUDIV_H

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

namespace plotterFunctions
{
    class Taudiv {
    private:
        std::shared_ptr<TopTagger> ttPtr_mine;
        void generateTaudiv(NTupleReader& tr) {
          const auto& tau1                  = tr.getVec<data_t>("tau1");
          const auto& tau2                  = tr.getVec<data_t>("tau2");
          const auto& tau3                  = tr.getVec<data_t>("tau3");
          const auto& puppitau1             = tr.getVec<data_t>("puppitau1");
          const auto& puppitau2             = tr.getVec<data_t>("puppitau2");
          const auto& puppitau3             = tr.getVec<data_t>("puppitau3");
          const auto& softDropMass          = tr.getVec<data_t>("softDropMass");
          const auto& puppisoftDropMass     = tr.getVec<data_t>("puppisoftDropMass");
          const auto& jetsLVec              = tr.getVec<TLorentzVector>("jetsLVec");
          const auto& ak8JetsLVec           = tr.getVec<TLorentzVector>("ak8JetsLVec");
          const auto& puppiJetsLVec         = tr.getVec<TLorentzVector>("puppiJetsLVec");
          const auto& genDecayLVec          = tr.getVec<TLorentzVector>("genDecayLVec");
          const auto& genDecayPdgIdVec      = tr.getVec<int>("genDecayPdgIdVec");
          const auto& genDecayIdxVec        = tr.getVec<int>("genDecayIdxVec");
          const auto& genDecayMomIdxVec     = tr.getVec<int>("genDecayMomIdxVec");
          const auto& nTopCandSortedCnt     = tr.getVar<int>("nTopCandSortedCntZinv");

          auto* puppiLVecLoose_top          = new std::vector<TLorentzVector>();
          auto* puppiLVectight_top          = new std::vector<TLorentzVector>();
          auto* puppiLVecLoose_w            = new std::vector<TLorentzVector>();
          auto* puppiLVectight_w            = new std::vector<TLorentzVector>();
          auto* puppitau2Dtau1              = new std::vector<data_t>();
          auto* puppitau3Dtau2              = new std::vector<data_t>();
          auto* puppitau2Dtau1_SDM          = new std::vector<data_t>();
          auto* puppitau3Dtau2_SDM          = new std::vector<data_t>();
          auto* hadWLVec                    = new std::vector<TLorentzVector>();

          //std::shared_ptr<TopTagger> ttPtr;
          //const TopTaggerResults& ttr = ttPtr->getResults();
          int monoJet=0;
          int diJet=0;
          int triJet=0;
          //TopTagger tt;
          //tt.setCfgFile("TopTagger.cfg");
          //const TopTaggerResults& ttr = ttPtr_mine.getResults();
          const TopTaggerResults& ttr =ttPtr_mine->getResults();
          std::vector<TopObject*> Ntop = ttr.getTops();
          for(int i=0; i<nTopCandSortedCnt; i++){
              if(Ntop[i]->getNConstituents() == 1) monoJet++;
              else if(Ntop[i]->getNConstituents() == 2) diJet++;
              else if(Ntop[i]->getNConstituents() == 3) triJet++;
              //std::cout<<monoJet<<std::endl;
          }
          //std::cout<<"Ntop: " << nTopCandSortedCnt<<std::endl;
          //std::cout<<"Monojet: " << monoJet<<std::endl;
          //std::cout<<"Dijet: " << diJet<<std::endl;
          //std::cout<<"Trijet: " << triJet<<std::endl;
          
          const int& nJetsAk8 = ak8JetsLVec.size(); 
          const int& nJetsPuppi = puppiJetsLVec.size();
          tr.registerDerivedVar("nJetsAk8", nJetsAk8);
          tr.registerDerivedVar("nJetsPuppi", nJetsPuppi);
          tr.registerDerivedVar("typeMono",monoJet);
          tr.registerDerivedVar("typeDi",diJet);
          tr.registerDerivedVar("typeTri",triJet);
          
          //for(unsigned int i=1; i< ak8JetsLVec.size(); i++){
          //std::cout<<"AK8 size "<<njetsAk8 << std::endl;
          //std::cout<<"AK8 pt "<<ak8JetsLVec[i].Pt() << std::endl;
          //}
           
          if(puppitau2.size()!=0 && puppitau1.size()!=0 && puppitau2.size()==puppitau1.size()){
              for(int iJet = 0; iJet < nJetsPuppi; ++iJet){
                  puppitau2Dtau1->push_back(puppitau2[iJet]/(puppitau1[iJet]));
              }
          }
          else { 
              puppitau2Dtau1->push_back( -1);
          }

          if(puppitau2.size()!=0 && puppitau3.size()!=0 && puppitau2.size()==puppitau3.size()){
              for(int iJet = 0; iJet < nJetsPuppi; ++iJet){
                  puppitau3Dtau2->push_back(puppitau3[iJet]/(puppitau2[iJet]));
              }
          }
          else {
              puppitau3Dtau2->push_back( -1);
          }
          tr.registerDerivedVec("puppitau2Dtau1", puppitau2Dtau1);
          tr.registerDerivedVec("puppitau3Dtau2", puppitau3Dtau2);
          
          ///WTagging
         /* 
          for(int tau = 0; tau < (*puppitau2Dtau1).size(); ++tau){
             if (puppisoftDropMass[tau]>65 && puppisoftDropMass[tau]<100){
                 // push back tau variables after mass cut
                 puppitau2Dtau1_SDM->push_back(puppitau2Dtau1->at(tau));

                 if ((*puppitau2Dtau1)[tau] >= 0 && (*puppitau2Dtau1)[tau] < 0.6){ // loose
                     puppiLVecLoose_w->push_back(puppiJetsLVec[tau]);  
                     //std::cout <<"PT_puupi"<< (*puppiLVectight_w).size()  << std::endl;    // (*puppiLVectight_w)[0].Pt() 

                     if ((*puppitau2Dtau1)[tau] < 0.45){ // tight
                         puppiLVectight_w->push_back(puppiJetsLVec[tau]); 
                     }
                 }
             } 
          }
          */
          //Top 1%
          /*
          for(int tau = 0; tau < (*puppitau3Dtau2).size(); ++tau){
              if (puppisoftDropMass[tau]>105 && puppisoftDropMass[tau]<210){
                  puppitau3Dtau2_SDM->push_back(puppitau3Dtau2->at(tau));
                  
                  if ((*puppitau3Dtau2)[tau] >= 0 && (*puppitau3Dtau2)[tau] < 0.65){
                      puppiLVecLoose_top->push_back(puppiJetsLVec[tau]);

                      if ((*puppitau3Dtau2)[tau] < 0.54){ 
                          puppiLVectight_top->push_back(puppiJetsLVec[tau]);
                      }
                  }
              }
          }
          */
          tr.registerDerivedVec("puppiLVectight_top", puppiLVectight_top);
          tr.registerDerivedVec("puppiLVecLoose_top", puppiLVecLoose_top);
          tr.registerDerivedVec("puppiLVectight_w", puppiLVectight_w);
          tr.registerDerivedVec("puppiLVecLoose_w", puppiLVecLoose_w);
          tr.registerDerivedVec("puppitau2Dtau1_SDM", puppitau2Dtau1_SDM);
          tr.registerDerivedVec("puppitau3Dtau2_SDM", puppitau3Dtau2_SDM);

          //(*hadWLVec) = genUtility::GetHadWLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
        }

    public:
        Taudiv(std::shared_ptr<TopTagger> ttPtr) { 
          //std::cout << "OMG! OMG! OMG! What's the STD?" << std::endl;
         ttPtr_mine = ttPtr;
        }
        ~Taudiv() {}
        void operator()(NTupleReader& tr)
        {
          generateTaudiv(tr);
        }       

    };
}

#endif
