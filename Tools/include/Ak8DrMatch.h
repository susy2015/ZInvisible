#ifndef AK8DRMATCH_H
#define AK8DRMATCH_H

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
     class Ak8DrMatch {
     private:
         void generateAk8DrMatch(NTupleReader& tr) {
             const auto& jetsLVec           = tr.getVec<TLorentzVector>("jetsLVec");
             const auto& puppiJetsLVec      = tr.getVec<TLorentzVector>("puppiJetsLVec");
             const auto& puppiLVectight_top = tr.getVec<TLorentzVector>("puppiLVectight_top");
             const auto& puppiLVecLoose_top = tr.getVec<TLorentzVector>("puppiLVecLoose_top");
             const auto& puppiLVectight_w   = tr.getVec<TLorentzVector>("puppiLVectight_w");
             const auto& puppiLVecLoose_w   = tr.getVec<TLorentzVector>("puppiLVecLoose_w");  

             int nJetsAK41_min = 0;
             int nJetsAK41_med = 0; 
             int nJetsAK41_lar = 0;
             int nJetsPuppi_T1_min = 0;
             int nJetsPuppi_T1_med = 0;
             int nJetsPuppi_T1_lar = 0;
             int nJetsPuppi_L1_min = 0;
             int nJetsPuppi_L1_med = 0;
             int nJetsPuppi_L1_lar = 0;
             int nJetsAK42_min = 0;
             int nJetsAK42_med = 0;
             int nJetsAK42_lar = 0;
 
             auto* ak81dRMin = new std::vector<data_t>();
             auto* ak82dRMin = new std::vector<data_t>(); 
             auto* puppi_top_L_1dRMin = new std::vector<data_t>();
             auto* puppi_top_L_2dRMin = new std::vector<data_t>();
             auto* puppi_top_T_1dRMin = new std::vector<data_t>();
             auto* puppi_top_T_2dRMin = new std::vector<data_t>(); 
             for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
             {
                 if(puppiJetsLVec.size() >= 1) ak81dRMin->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiJetsLVec[0]));
                 if(puppiJetsLVec.size() >= 2) ak82dRMin->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiJetsLVec[1]));
                 //if(puppiLVectight_top.size() >= 1) puppi_top_L_1dRMin-> push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVectight_top[0]));
                 //if(puppiLVectight_top.size() >= 2) puppi_top_L_2dRMin ->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVectight_top[1]));
                 //if(puppiLVecLoose_top.size() >= 1) puppi_top_T_1dRMin-> push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVecLoose_top[0]));
                 //if(puppiLVecLoose_top.size() >= 2) puppi_top_T_2dRMin ->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVecLoose_top[1]));
                 std::sort( ak81dRMin->begin(),ak81dRMin->end() );
                 std::sort( ak82dRMin->begin(),ak82dRMin->end() );
             }
             tr.registerDerivedVec("ak81dRMin", ak81dRMin);
             tr.registerDerivedVec("ak82dRMin", ak82dRMin);
             //tr.registerDerivedVec("puppi_top_L_1dRMin", puppi_top_L_1dRMin);
             //tr.registerDerivedVec("puppi_top_L_2dRMin", puppi_top_L_2dRMin);    
             //tr.registerDerivedVec("puppi_top_T_1dRMin", puppi_top_T_1dRMin);
             //tr.registerDerivedVec("puppi_top_T_2dRMin", puppi_top_T_2dRMin);     
 
             for(int iJet1 = 0; iJet1 < ak81dRMin->size(); ++iJet1)
             {
                 if(ak81dRMin->at(iJet1) <=0.2) nJetsAK41_min++;
                 if(ak81dRMin->at(iJet1) <=0.4 && ak81dRMin->at(iJet1) > 0.2) nJetsAK41_med++;
                 if(ak81dRMin->at(iJet1) <=0.8 && ak81dRMin->at(iJet1) > 0.4) nJetsAK41_lar++;
            
             } 
             for(int iJet1 = 0; iJet1 < puppi_top_L_1dRMin->size(); ++iJet1)
             {
                 if(puppi_top_L_1dRMin->at(iJet1) <=0.2) nJetsPuppi_L1_min++;
                 if(puppi_top_L_1dRMin->at(iJet1) <=0.4 && puppi_top_L_1dRMin->at(iJet1) > 0.2) nJetsPuppi_L1_med++;
                 if(puppi_top_L_1dRMin->at(iJet1) <=0.8 && puppi_top_L_1dRMin->at(iJet1) > 0.4) nJetsPuppi_L1_lar++;

             }
             for(int iJet1 = 0; iJet1 < puppi_top_T_1dRMin->size(); ++iJet1)
             {
                 if(puppi_top_T_1dRMin->at(iJet1) <=0.2) nJetsPuppi_T1_min++;
                 if(puppi_top_T_1dRMin->at(iJet1) <=0.4 && puppi_top_T_1dRMin->at(iJet1) > 0.2) nJetsPuppi_T1_med++;
                 if(puppi_top_T_1dRMin->at(iJet1) <=0.8 && puppi_top_T_1dRMin->at(iJet1) > 0.4) nJetsPuppi_T1_lar++;

             }
             for(int iJet2 = 0; iJet2 < ak82dRMin->size(); ++iJet2)
             {
                 if(ak82dRMin->at(iJet2) <=0.2) nJetsAK42_min++;
                 if(ak82dRMin->at(iJet2) <=0.4 && ak82dRMin->at(iJet2) > 0.2) nJetsAK42_med++;                         
                 if(ak82dRMin->at(iJet2) <=0.8 && ak82dRMin->at(iJet2) > 0.4) nJetsAK42_lar++;
   
             }          
             tr.registerDerivedVar("nJetsAK41_min",nJetsAK41_min);
             tr.registerDerivedVar("nJetsAK41_med",nJetsAK41_med);
             tr.registerDerivedVar("nJetsAK41_lar",nJetsAK41_lar);
             tr.registerDerivedVar("nJetsPuppi_L1_min",nJetsPuppi_L1_min);
             tr.registerDerivedVar("nJetsPuppi_L1_med",nJetsPuppi_L1_med);
             tr.registerDerivedVar("nJetsPuppi_L1_lar",nJetsPuppi_L1_lar);
             tr.registerDerivedVar("nJetsPuppi_T1_min",nJetsPuppi_T1_min);
             tr.registerDerivedVar("nJetsPuppi_T1_med",nJetsPuppi_T1_med);
             tr.registerDerivedVar("nJetsPuppi_T1_lar",nJetsPuppi_T1_lar);
             tr.registerDerivedVar("nJetsAK42_min",nJetsAK42_min);
             tr.registerDerivedVar("nJetsAK42_med",nJetsAK42_med);
             tr.registerDerivedVar("nJetsAK42_lar",nJetsAK42_lar);


             // Also start looking at subjet information
             const std::vector<TLorentzVector>& puppiSubJetsLVec  = tr.getVec<TLorentzVector>("puppiSubJetsLVec");
             const std::vector<data_t>& puppiSubJetsBdisc = tr.getVec<data_t>("puppiSubJetsBdisc");

             // For each tagged top/W, find the corresponding subjets
             /*
             std::vector< std::vector<TLorentzVector> > W_subjets;
             std::vector<data_t>* W_subjets_pt_reldiff = new std::vector<data_t>();
             for( TLorentzVector myW : puppiLVectight_w)
             {
                 std::vector<TLorentzVector> myW_subjets;
                 int i = 0;
                 for(TLorentzVector puppiSubJet : puppiSubJetsLVec)
                 {
                     data_t myDR = ROOT::Math::VectorUtil::DeltaR(myW, puppiSubJet);
                     if (myDR < 0.8)
                     {
                         myW_subjets.push_back(puppiSubJet);
                     }
                     ++i;
                 }
                 // If more than 2 matches, find the best combination of two subjets by checking diff in 4-vector
                 if (myW_subjets.size() > 2) {
                     data_t min_diff = 999999.;
                     int min_j=0, min_k=1;
                     for (int j=0 ; j<myW_subjets.size(); ++j)
                     {
                         for (int k=j+1; k<myW_subjets.size(); ++k)
                         {
                             TLorentzVector diff_LV = myW - myW_subjets[j] - myW_subjets[k];
                             data_t diff = abs(diff_LV.M());
                             if(diff < min_diff)
                             {
                                 min_diff = diff;
                                 min_j = j;
                                 min_k = k;
                             }
                         }
                     }
                     std::vector<TLorentzVector> mynewW_subjets = {myW_subjets[min_j], myW_subjets[min_k]};
                     W_subjets.push_back(mynewW_subjets);
                     W_subjets_pt_reldiff->push_back( ((myW_subjets[min_j]+myW_subjets[min_k]).Pt()-myW.Pt())/myW.Pt());
                 } else {
                     W_subjets.push_back(myW_subjets);
                     W_subjets_pt_reldiff->push_back( ((myW_subjets[0]+myW_subjets[1]).Pt()-myW.Pt())/myW.Pt());
                 }
             }
             tr.registerDerivedVec("W_subjets_pt_reldiff", W_subjets_pt_reldiff);
             */
             // For each tagged top/W, find the corresponding subjets
             std::vector< std::vector< TLorentzVector> > top_subjets;
             std::vector<data_t>* top_subjets_pt_reldiff = new std::vector<data_t>();
             /*
             for( TLorentzVector mytop : puppiLVectight_top)
             {
                 std::vector<TLorentzVector> mytop_subjets;
                 int i = 0;
                 for(TLorentzVector puppiSubJet : puppiSubJetsLVec)
                 {
                     data_t myDR = ROOT::Math::VectorUtil::DeltaR(mytop, puppiSubJet);
                     if (myDR < 0.8)
                     {
                         mytop_subjets.push_back(puppiSubJet);
                     }
                     ++i;
                 }
                 // If more than 2 matches, find the best combination of two subjets
                 if (mytop_subjets.size() > 2) {
                     data_t min_diff = 999999.;
                     int min_j=0, min_k=1;
                     for (int j=0 ; j<mytop_subjets.size(); ++j)
                     {
                         for (int k=j+1; k<mytop_subjets.size(); ++k)
                         {
                             TLorentzVector diff_LV = mytop - mytop_subjets[j] - mytop_subjets[k];
                             data_t diff = abs(diff_LV.M());
                             if(diff < min_diff)
                             {
                                 min_diff = diff;
                                 min_j = j;
                                 min_k = k;
                             }
                         }
                     }
                     std::vector<TLorentzVector> mynewtop_subjets = {mytop_subjets[min_j], mytop_subjets[min_k]};
                     top_subjets.push_back(mynewtop_subjets);
                     top_subjets_pt_reldiff->push_back( ((mytop_subjets[min_j]+mytop_subjets[min_k]).Pt()-mytop.Pt())/mytop.Pt());
                 } else {
                     top_subjets.push_back(mytop_subjets);
                     top_subjets_pt_reldiff->push_back( ((mytop_subjets[0]+mytop_subjets[1]).Pt()-mytop.Pt())/mytop.Pt());
                 }
             }
             tr.registerDerivedVec("top_subjets_pt_reldiff", top_subjets_pt_reldiff);
             */
             // Figure out gen matching..
             const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
             const std::vector<int>& genDecayIdxVec          = tr.getVec<int>("genDecayIdxVec");
             const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("genDecayMomIdxVec");
             const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

             auto* gentop_match                                 = new std::vector<bool>(); // helpful to make plots of matched and unmatched number of tops
             auto* dR_top_gentop                                = new std::vector<data_t>(); 
             auto* dR_AK4_topsubjet_genmatched                  = new std::vector<data_t>(); 
             auto* dR_AK4_top_genmatched                        = new std::vector<data_t>(); 
             auto* top_N_AK4_matched_genmatched                 = new std::vector<int>(); 
             auto* top_N_AK4_matched_notgenmatched              = new std::vector<int>(); 
             auto* top_N_AK4_notmatched_genmatched              = new std::vector<int>(); 
             auto* top_N_AK4_notmatched_notgenmatched           = new std::vector<int>(); 
             auto* top_N_AK4_matched_genmatched_0p6             = new std::vector<int>(); 
             auto* top_N_AK4_matched_notgenmatched_0p6          = new std::vector<int>(); 
             auto* top_N_AK4_notmatched_genmatched_0p6          = new std::vector<int>(); 
             auto* top_N_AK4_notmatched_notgenmatched_0p6       = new std::vector<int>(); 
             auto* top_N_AK4_matched_genmatchedother            = new std::vector<int>(); 
             auto* top_N_AK4_matched_notgenmatchedother         = new std::vector<int>(); 
             auto* top_N_AK4_matched_genmatchedother_0p6        = new std::vector<int>(); 
             auto* top_N_AK4_matched_notgenmatchedother_0p6     = new std::vector<int>(); 
             if(tr.checkBranch("genDecayPdgIdVec") && &genDecayLVec != nullptr)
             {
                 // For each tagged top, find the matching gen particles

                 // These are the hadronically decaying top quarks in the event:
                 std::vector<TLorentzVector> hadtopLVec = genUtility::GetHadTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
                 std::vector< std::vector<TLorentzVector> > hadtopdauLVec;
                 for(TLorentzVector hadtop : hadtopLVec)
                 {
                     hadtopdauLVec.push_back(genUtility::GetTopdauLVec(hadtop, genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec));
                 }

                 // check all tagged tops
                 /*
                 for(unsigned int imytop=0; imytop<puppiLVectight_top.size(); ++imytop) 
                 {
                     TLorentzVector mytop = puppiLVectight_top[imytop];
                     //std::cout << "Mytop info: " << mytop.Pt() << " " << mytop.Eta() << " " << mytop.Phi() << std::endl;
                     // For now find the closest hadtop in deltaR
                     TLorentzVector temp_gentop_match_LV;
                     double min_DR = 99.;
                     int matched_hadtop_index = -1;
                     for(unsigned int myhadtop_i=0; myhadtop_i<hadtopLVec.size(); ++myhadtop_i)
                     {
                         TLorentzVector myhadtop = hadtopLVec[myhadtop_i];
                         double DR_top = ROOT::Math::VectorUtil::DeltaR(mytop, myhadtop);
                         if (DR_top < min_DR) 
                         {
                             temp_gentop_match_LV = myhadtop;
                             min_DR = DR_top;
                             matched_hadtop_index = myhadtop_i;
                         }
                     }
                     dR_top_gentop->push_back(min_DR);
                     // DR should be small for it to actually be a match
                     if(min_DR < 0.4)
                     {
                         //std::cout << "Mytop info: " << mytop.Pt() << " " << mytop.Eta() << " " << mytop.Phi() << std::endl;
                         gentop_match->push_back(true);
                         // Now find the gen daughters for this gentop
                         std::vector<TLorentzVector> gentopdauLVec = hadtopdauLVec[matched_hadtop_index];

                         // Now we have the tagged top (mytop), the gen had top (temp_gentop_match_LV), and the gen daughters (gentopdauLVec)
                         // ready for some matching FUN!

                         // Removing AK4 jets based on DR matching with subjets of tagged top
                         std::vector<TLorentzVector> mysubjets = top_subjets[imytop];
                         std::vector<int> ak4_removed;
                         if(mysubjets.size() != 2)
                             std::cout << "Attention: found " << mysubjets.size() << " subjets instead of 2" << std::endl;
                         //std::cout << "Subjet 0: " << mysubjets[0].Pt() << " " << mysubjets[0].Eta() << " " << mysubjets[0].Phi() << std::endl;
                         //std::cout << "Subjet 0: " << mysubjets[1].Pt() << " " << mysubjets[1].Eta() << " " << mysubjets[1].Phi() << std::endl;
                         
                         // some counters
                         int N_AK4_matched_genmatched = 0;
                         int N_AK4_matched_notgenmatched = 0;
                         int N_AK4_notmatched_genmatched = 0;
                         int N_AK4_notmatched_notgenmatched = 0;
                         // some counters
                         int N_AK4_matched_genmatched_0p6 = 0;
                         int N_AK4_matched_notgenmatched_0p6 = 0;
                         int N_AK4_notmatched_genmatched_0p6 = 0;
                         int N_AK4_notmatched_notgenmatched_0p6 = 0;
                         // matched to another gentop
                         int N_AK4_matched_genmatchedother = 0;
                         int N_AK4_matched_notgenmatchedother = 0;
                         int N_AK4_matched_genmatchedother_0p6 = 0;
                         int N_AK4_matched_notgenmatchedother_0p6 = 0;

                         for (unsigned int j=0; j<jetsLVec.size(); ++j)
                         {
                             double DR1 = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mysubjets[0]);
                             double DR2 = DR1;
                             if(mysubjets.size()>1) {
                             double DR2 = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mysubjets[1]);
                             }
                             //std::cout << "DR1, DR2: " << DR1 << " " << DR2 << std::endl;
                             // Check if it matches a gen daughter
                             bool genmatch = false;
                             for (TLorentzVector gendau : gentopdauLVec)
                             {
                                 double DR_AK4_gen = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], gendau);
                                 //std::cout << "gen DR " << DR_AK4_gen << std::endl;
                                 if (DR_AK4_gen < 0.4)
                                 {
                                     // matches gendaughter
                                     genmatch = true;
                                     break;
                                 }
                             }
                             if(genmatch){
                                 dR_AK4_topsubjet_genmatched->push_back(std::min(DR1,DR2));
                                 dR_AK4_top_genmatched->push_back(ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mytop));
                             }
                             // should merge this with 'genmatch' finding...
                             bool genmatch_other = false;
                             for (unsigned int other=0; other<hadtopdauLVec.size() && !genmatch_other; ++other)
                             {
                                 if(other == matched_hadtop_index)
                                     continue;
                                 for (TLorentzVector gendau : hadtopdauLVec[other])
                                 {
                                     double DR_AK4_gen = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], gendau);
                                     //std::cout << "gen DR " << DR_AK4_gen << std::endl;
                                     if (DR_AK4_gen < 0.4)
                                     {
                                         // matches gendaughter
                                         genmatch_other = true;
                                         break;
                                     }
                                 }
                             }

                             bool subjetmatch = false;
                             bool subjetmatch_0p6 = false;
                             if (DR1 < 0.4 || DR2 < 0.4)
                             {
                                 //std::cout << "Found AK4 jet matching a subjet" << std::endl;
                                 // found a match
                                 subjetmatch = true;
                                 ak4_removed.push_back(j);
                             }
                             if (DR1 < 0.6 || DR2 < 0.6)
                             {
                                 //std::cout << "Found AK4 jet matching a subjet" << std::endl;
                                 // found a match
                                 subjetmatch_0p6 = true;
                             }


                             if(genmatch){
                                 if(subjetmatch)
                                     N_AK4_matched_genmatched++;
                                 else
                                     N_AK4_notmatched_genmatched++;
                                 if(subjetmatch_0p6)
                                     N_AK4_matched_genmatched_0p6++;
                                 else
                                     N_AK4_notmatched_genmatched_0p6++;
                             } else { // Not genmatched to any of the correct top daughters
                                 if(subjetmatch)
                                     N_AK4_matched_notgenmatched++;
                                 else
                                     N_AK4_notmatched_notgenmatched++;
                                 if(subjetmatch_0p6)
                                     N_AK4_matched_notgenmatched_0p6++;
                                 else
                                     N_AK4_notmatched_notgenmatched_0p6++;

                                 if(genmatch_other){
                                     if(subjetmatch)
                                         N_AK4_matched_genmatchedother++;
                                     if(subjetmatch_0p6)
                                         N_AK4_matched_genmatchedother_0p6++;
                                 } else {
                                     if(subjetmatch)
                                         N_AK4_matched_notgenmatchedother++;
                                     if(subjetmatch_0p6)
                                         N_AK4_matched_notgenmatchedother_0p6++;
                                 }

                             }




                         }
                         top_N_AK4_matched_genmatched->push_back(N_AK4_matched_genmatched);
                         top_N_AK4_matched_notgenmatched->push_back(N_AK4_matched_notgenmatched);
                         top_N_AK4_notmatched_genmatched->push_back(N_AK4_notmatched_genmatched);
                         top_N_AK4_notmatched_notgenmatched->push_back(N_AK4_notmatched_notgenmatched);
                         top_N_AK4_matched_genmatched_0p6->push_back(N_AK4_matched_genmatched_0p6);
                         top_N_AK4_matched_notgenmatched_0p6->push_back(N_AK4_matched_notgenmatched_0p6);
                         top_N_AK4_notmatched_genmatched_0p6->push_back(N_AK4_notmatched_genmatched_0p6);
                         top_N_AK4_notmatched_notgenmatched_0p6->push_back(N_AK4_notmatched_notgenmatched_0p6);

                         top_N_AK4_matched_genmatchedother->push_back(N_AK4_matched_genmatchedother);
                         top_N_AK4_matched_notgenmatchedother->push_back(N_AK4_matched_notgenmatchedother);
                         top_N_AK4_matched_genmatchedother_0p6->push_back(N_AK4_matched_genmatchedother_0p6);
                         top_N_AK4_matched_notgenmatchedother_0p6->push_back(N_AK4_matched_notgenmatchedother_0p6);
                         
                     } else // No match
                     { 
                         gentop_match->push_back(false);
                     }

                 }
             */
             }
             
             tr.registerDerivedVec("gentop_match", gentop_match);
             tr.registerDerivedVec("dR_top_gentop", dR_top_gentop);
             tr.registerDerivedVec("dR_AK4_topsubjet_genmatched", dR_AK4_topsubjet_genmatched);
             tr.registerDerivedVec("dR_AK4_top_genmatched", dR_AK4_top_genmatched);
             tr.registerDerivedVec("top_N_AK4_matched_genmatched", top_N_AK4_matched_genmatched);
             tr.registerDerivedVec("top_N_AK4_matched_notgenmatched", top_N_AK4_matched_notgenmatched);
             tr.registerDerivedVec("top_N_AK4_notmatched_genmatched", top_N_AK4_notmatched_genmatched);
             tr.registerDerivedVec("top_N_AK4_notmatched_notgenmatched", top_N_AK4_notmatched_notgenmatched);
             tr.registerDerivedVec("top_N_AK4_matched_genmatched_0p6", top_N_AK4_matched_genmatched_0p6);
             tr.registerDerivedVec("top_N_AK4_matched_notgenmatched_0p6", top_N_AK4_matched_notgenmatched_0p6);
             tr.registerDerivedVec("top_N_AK4_notmatched_genmatched_0p6", top_N_AK4_notmatched_genmatched_0p6);
             tr.registerDerivedVec("top_N_AK4_notmatched_notgenmatched_0p6", top_N_AK4_notmatched_notgenmatched_0p6);

             tr.registerDerivedVec("top_N_AK4_matched_genmatchedother", top_N_AK4_matched_genmatchedother);
             tr.registerDerivedVec("top_N_AK4_matched_notgenmatchedother", top_N_AK4_matched_notgenmatchedother);
             tr.registerDerivedVec("top_N_AK4_matched_genmatchedother_0p6", top_N_AK4_matched_genmatchedother_0p6);
             tr.registerDerivedVec("top_N_AK4_matched_notgenmatchedother_0p6", top_N_AK4_matched_notgenmatchedother_0p6);

         }

     public:
         Ak8DrMatch() {
         }
         ~Ak8DrMatch() {}
         void operator()(NTupleReader& tr)
         {
             generateAk8DrMatch(tr);
         }
     };
}

#endif
