#ifndef TRIGGERINFO_H 
#define TRIGGERINFO_H 

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
    class TriggerInfo
    {
    private:
        int indexMuTrigger;
        int indexElecTrigger;
        int indexMETMHTTrigger;
        int indexPhotonTrigger;
        bool miniTuple_, noMC_;

        double GetMuonTriggerEff(const double& muEta) 
        {
            if (-2.6 <= muEta && muEta < -2.2) return 0.7861842;
            else if(-2.2 <= muEta && muEta < -1.8) return 0.8233438;
            else if(-1.8 <= muEta && muEta < -1.4) return 0.8151685;
            else if(-1.4 <= muEta && muEta < -1.0) return 0.8991723;
            else if(-1.0 <= muEta && muEta < -0.6) return 0.9125786;
            else if(-0.6 <= muEta && muEta < -0.2) return 0.8880085;
            else if(-0.2 <= muEta && muEta <  0.2) return 0.9334851;
            else if( 0.2 <= muEta && muEta <  0.6) return 0.8857523;
            else if( 0.6 <= muEta && muEta <  1.0) return 0.9052119;
            else if( 1.0 <= muEta && muEta <  1.4) return 0.9004312;
            else if( 1.4 <= muEta && muEta <  1.8) return 0.8384009;
            else if( 1.8 <= muEta && muEta <  2.2) return 0.8218332;
            else if( 2.2 <= muEta && muEta <  2.6) return 0.7781818;
            else                                   return 0.000;
        }

        double GetTriggerEffWeight(const double& met, const double& ht) 
        {
            if (ht<1000)
            {
                if (met<25) return 0.001542561;
                else if (met<50) return 0.003222389;
                else if (met<75) return 0.00987073;
                else if (met<100) return 0.03865682;
                else if (met<125) return 0.1387231;
                else if (met<150) return 0.3564816;
                else if (met<175) return 0.6276442;
                else if (met<200) return 0.8154821;
                else if (met<275) return 0.9340538;
                else if (met<400) return 0.9858562; 
                else if (met<600) return 0.9931507;
                else if (met<1000) return 1.00;
                else return 1.00;
            } 
            else 
            {
                if (met<25) return  0.02067183;
                else if (met<50) return 0.02504944;
                else if (met<75) return 0.04486466;
                else if (met<100) return 0.07434402;
                else if (met<125) return 0.1518288;
                else if (met<150) return 0.2802669;
                else if (met<175) return 0.4642409;
                else if (met<200) return 0.6596434;
                else if (met<275) return 0.8510453;
                else if (met<400) return 0.9563492;
                else if (met<600) return 0.9874214;
                else if (met<1000) return 0.9736842; 
                else return 0.9736842;
            }
        }
        double GetTriggerEffStatUncHi(const double& met, const double& ht) 
        {
            if (ht<1000)
            {
                if (met<25) return 0.0001251554;
                else if (met<50) return 0.0001310897;
                else if (met<75) return 0.0002597269;
                else if (met<100) return 0.0006525702;
                else if (met<125) return 0.001545856;
                else if (met<150) return 0.002821274;
                else if (met<200) return 0.003691577;
                else if (met<275) return 0.003877182;
                else if (met<400) return 0.002294442; 
                else if (met<600) return 0.002045071;
                else if (met<1000) return 0.003725375;
                else return 0.00;
            } 
            else 
            {
                if (met<25) return 0.004283915;
                else if (met<50) return 0.003169914;
                else if (met<75) return 0.004349597;
                else if (met<100) return 0.006241982;
                else if (met<125) return 0.01001983;
                else if (met<150) return 0.01455422;
                else if (met<175) return 0.0183275;
                else if (met<200) return 0.01960093;
                else if (met<275) return 0.01062354;
                else if (met<400) return 0.007445741;
                else if (met<600) return 0.006010458;
                else if (met<1000) return 0.01697945; 
                else return 0.01697945;
            }
        }
        double GetTriggerEffStatUncLo(const double& met, const double& ht) 
        {
            if (ht<1000)
            {
                if (met<25) return 0.0001160878;
                else if (met<50) return 0.000126075;
                else if (met<75) return 0.0002532144;
                else if (met<100) return 0.000642253;
                else if (met<125) return 0.001531628;
                else if (met<150) return 0.002811409;
                else if (met<175) return 0.003706407;
                else if (met<200) return 0.003940439;
                else if (met<275) return 0.00236968;
                else if (met<400) return 0.002358961;
                else if (met<600) return 0.006617554;
                else if (met<1000) return 0.1422293;
                else return 0.1422293;
            }
           
            else 
            {
                if (met<25) return 0.003609465;
                else if (met<50) return 0.002838673;
                else if (met<75) return 0.003996443;
                else if (met<100) return 0.005811049;
                else if (met<125) return 0.009521872;
                else if (met<150) return 0.01412113;
                else if (met<175) return 0.01823465;
                else if (met<200) return 0.02013986;
                else if (met<275) return 0.01126014;
                else if (met<400) return 0.008759573;
                else if (met<600) return 0.009833846;
                else if (met<1000) return 0.03365661; 
                else return 0.03365661;
            }
        }
        double GetTriggerEffSystUncHi(const double& met, const double& ht) 
        {
            return 0.0;
            /* if (met<100) return 0.0272; */
            /* else if (met<150) return 0.0872; */
            /* else if (met<175) return 0.1505; */
            /* else if (met<200) return 0.0423; */
            /* else if (met<275) return 0.0112; */
            /* else if (met<400) return 0.0001;  */
            /* else return 0.0; */
        }
        double GetTriggerEffSystUncLo(const double& met, const double& ht) 
        {
            return 0.0;
            /* if (met<100) return 0.0120; */
            /* else if (met<150) return 0.0872; */
            /* else if (met<175) return 0.1505; */
            /* else if (met<200) return 0.0792; */
            /* else if (met<275) return 0.0112; */
            /* else if (met<400) return 0.0001;  */
            /* else return 0.0018; */
        }

        void triggerInfo(NTupleReader& tr)
        {
            const std::vector<std::string>& triggerNames = tr.getVec<std::string>("TriggerNames");
            const std::vector<int>& passTrigger          = tr.getVec<int>("PassTrigger");

            bool passMuTrigger = false;
            bool passElecTrigger = false;
            bool passMETMHTTrigger = false;
            bool passPhotonTrigger = false;

            const std::string muTrigName = "HLT_Mu50_v";//"HLT_Mu45_eta2p1_v";
            const std::string elecTrigName = "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v";
            const std::string metmhtTrigName = "HLT_PFMET110_PFMHT110_IDTight_v";
            const std::string photonTrigName = "HLT_Photon175_v";

            // Find the index of our triggers if we don't know them already
            if(indexMuTrigger == -1 || indexElecTrigger == -1 || indexMETMHTTrigger == -1 || indexPhotonTrigger == -1)
            {
                for(int i = 0; i < triggerNames.size(); ++i)
                {
                  if(triggerNames[i].find(muTrigName) != std::string::npos)
                    {
                      indexMuTrigger = i;
                    }
                  else if(triggerNames[i].find(elecTrigName) != std::string::npos)
                    {
                      indexElecTrigger = i;
                    }
                  else if(triggerNames[i].find(metmhtTrigName) != std::string::npos)
                    {
                      indexMETMHTTrigger = i;
                    }
                  else if(triggerNames[i].find(photonTrigName) != std::string::npos)
                    {
                      indexPhotonTrigger = i;
                    }

                }
            }
            if(indexMuTrigger != -1 && indexElecTrigger != -1 && indexPhotonTrigger != -1)
            {
                // Check if the event passes the trigger, and double check that we are looking at the right trigger
                if(triggerNames[indexMuTrigger].find(muTrigName) != std::string::npos && passTrigger[indexMuTrigger])
                    passMuTrigger = true;
                if(triggerNames[indexElecTrigger].find(elecTrigName) != std::string::npos && passTrigger[indexElecTrigger])
                    passElecTrigger = true;
                if(triggerNames[indexMETMHTTrigger].find(metmhtTrigName) != std::string::npos && passTrigger[indexMETMHTTrigger])
                    passMETMHTTrigger = true;
                if(triggerNames[indexPhotonTrigger].find(photonTrigName) != std::string::npos && passTrigger[indexPhotonTrigger])
                  passPhotonTrigger = true;
            }
            else
            {
                std::cout << "Could not find trigger in the list of trigger names" << std::endl;
            }

            bool passSearchTrigger = false;
            for(int it = 0; it < triggerNames.size(); ++it)
            {
                if( triggerNames[it].find("HLT_PFMET170_NoiseCleaned_v")             != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET170_JetIdCleaned_v")             != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET170_HBHECleaned_v")              != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v")         != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != std::string::npos ||
                    triggerNames[it].find("HLT_Photon175_v") != std::string::npos
                    )
                {
                    if( passTrigger[it] ) 
                    {
                        passSearchTrigger = true;
                        break;
                    }
                }
            }

            tr.registerDerivedVar("passPhotonTrigger",passPhotonTrigger);
            tr.registerDerivedVar("passMuTrigger",passMuTrigger);
            tr.registerDerivedVar("passElecTrigger",passElecTrigger);
            tr.registerDerivedVar("passMETMHTTrigger",passMETMHTTrigger);
            tr.registerDerivedVar("passSearchTrigger",passSearchTrigger);
        }

        void triggerInfoMC(NTupleReader& tr)
        {
            const double& met                            = tr.getVar<double>("cleanMetPt");
            const double& ht                             = tr.getVar<double>("HTZinv");
            const std::vector<TLorentzVector>& cutMuVec  = tr.getVec<TLorentzVector>("cutMuVec");

            // MC trigger efficiencies
            double triggerEff = GetTriggerEffWeight(met,ht);
            double triggerEffStatUncUp = GetTriggerEffStatUncHi(met,ht);
            double triggerEffSystUncUp = GetTriggerEffSystUncHi(met,ht);
            double triggerEffUncUp     = TMath::Sqrt(triggerEffStatUncUp*triggerEffStatUncUp + triggerEffSystUncUp*triggerEffSystUncUp);
            double triggerEffStatUncDown = GetTriggerEffStatUncLo(met,ht);
            double triggerEffSystUncDown = GetTriggerEffSystUncLo(met,ht);
            double triggerEffUncDown     = TMath::Sqrt(triggerEffStatUncDown*triggerEffStatUncDown + triggerEffSystUncDown*triggerEffSystUncDown);

            //Calculate muon trigger weights
            double muTrigWgt = 0.0;
            if(cutMuVec.size() >= 2 && cutMuVec[0].Pt() > 50 && cutMuVec[1].Pt() > 50)
            {
                double muEff1 = GetMuonTriggerEff(cutMuVec[0].Eta());
                double muEff2 = GetMuonTriggerEff(cutMuVec[1].Eta());

                muTrigWgt = 1 - (1 - muEff1)*(1 - muEff2);
            }
            else if(cutMuVec.size() >= 1 && cutMuVec[0].Pt() > 50)
            {
                //For events with only 1 muon (emu events in particular or events with a subleading muon below 45 GeV) just use the single muon eff
                muTrigWgt = GetMuonTriggerEff(cutMuVec[0].Eta());
            }

            tr.registerDerivedVar("TriggerEffMC",triggerEff);
            tr.registerDerivedVar("TriggerEffUpMC",triggerEff+triggerEffUncUp);
            tr.registerDerivedVar("TriggerEffDownMC",triggerEff-triggerEffUncDown);

            tr.registerDerivedVar("muTrigWgt", muTrigWgt);
        }

    public:
        TriggerInfo(bool miniTuple = false, bool noMC = false)
        {
            indexMuTrigger = -1;
            indexElecTrigger = -1;
            indexMETMHTTrigger = -1;
            miniTuple_ = miniTuple;
            noMC_ = noMC;
        }

        void operator()(NTupleReader& tr)
        {
            if(!miniTuple_) triggerInfo(tr);
            if(!noMC_)       triggerInfoMC(tr);
        }

    };
}

#endif
