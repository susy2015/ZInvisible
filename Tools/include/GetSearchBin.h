#ifndef GETSEARCHBIN_H
#define GETSEARCHBIN_H

#include "TypeDefinitions.h"
#include "PhotonTools.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/SB2018.h"
#include "SusyAnaTools/Tools/SusyUtility.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"
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
    class GetSearchBin
    {
    private:

        void getSearchBin(NTupleReader& tr)
        {
            std::vector<std::string> JetPtCuts = {"_jetpt30"};
            // begin loop over jet pt cuts 
            for (const auto& suffix : JetPtCuts) 
            {
                // Note: Only HT, s_met, nJets, and dPhi are calculated with different jet pt cuts
                const auto& nMergedTops         = tr.getVar<int>("nMergedTops");
                const auto& nResolvedTops       = tr.getVar<int>("nResolvedTops");
                const auto& nWs                 = tr.getVar<int>("nWs");
                const auto& nBottoms            = tr.getVar<int>("nBottoms");
                const auto& nSoftBottoms        = tr.getVar<int>("nSoftBottoms");
                const auto& nJets               = tr.getVar<int>("nJets" + suffix);
                const auto& ht                  = tr.getVar<data_t>("HT" + suffix);
                const auto& met                 = tr.getVar<data_t>("MET_pt");
                const auto& ptb                 = tr.getVar<data_t>("ptb");
                const auto& mtb                 = tr.getVar<data_t>("mtb");
                const auto& ISRJetPt            = tr.getVar<data_t>("ISRJetPt");
                
                float mtb_cut = 175.0;
                
                //-----------------------------------------//
                //--- Updated Search Bins (August 2019) ---//
                //-----------------------------------------//
                // int SBv4_lowdm(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
                // int SBv4_highdm(float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
                int nSearchBinLowDM  = SBv4_lowdm(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
                int nSearchBinHighDM = SBv4_highdm(mtb, nJets, nBottoms, nMergedTops, nWs, nResolvedTops, ht, met);
                
                //------------------------------------------------//
                //--- Updated Validation Bins (September 2019) ---//
                //------------------------------------------------//
                // int SBv3_lowdm_validation(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
                // int SBv3_lowdm_validation_high_MET(int nb, int nSV, float ISRpt, float met)
                // int SBv3_highdm_validation(float mtb, int njets, int ntop, int nw, int nres, int nb, float met)
                int nValidationBinLowDM        = SBv3_lowdm_validation(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
                int nValidationBinLowDMHighMET = SBv3_lowdm_validation_high_MET(nBottoms, nSoftBottoms, ISRJetPt, met);
                int nValidationBinHighDM       = SBv3_highdm_validation(mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met); 
                
                tr.registerDerivedVar("nSearchBinLowDM"             + suffix, nSearchBinLowDM);
                tr.registerDerivedVar("nSearchBinHighDM"            + suffix, nSearchBinHighDM);
                tr.registerDerivedVar("nValidationBinLowDM"         + suffix, nValidationBinLowDM);
                tr.registerDerivedVar("nValidationBinLowDMHighMET"  + suffix, nValidationBinLowDMHighMET);
                tr.registerDerivedVar("nValidationBinHighDM"        + suffix, nValidationBinHighDM);
            }
        }

    public:

        GetSearchBin()
        {
            // void printJson(const std::string& fileName, const std::string& key, const std::string& title)
            const std::string fileName = "dc_BkgPred_BinMaps_master.json";
            SusyUtility::printJson(fileName, "binNum", "Search Bins");
        }
        
        ~GetSearchBin(){}

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };
}

#endif
