#ifndef GETSEARCHBIN_H
#define GETSEARCHBIN_H

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/SB2018.h"

namespace plotterFunctions
{
    class GetSearchBin
    {
    private:

        void getSearchBin(NTupleReader& tr)
        {
            const auto& nBottoms            = tr.getVar<int>("nBottoms");
            const auto& nSoftBottoms        = tr.getVar<int>("nSoftBottoms");
            const auto& nMergedTops         = tr.getVar<int>("nMergedTops");
            const auto& nJets               = tr.getVar<int>("nJets");
            const auto& nWs                 = tr.getVar<int>("nWs");
            const auto& nResolvedTops       = tr.getVar<int>("nResolvedTops");
            const auto& met                 = tr.getVar<data_t>("MET_pt");
            const auto& ht                  = tr.getVar<data_t>("HT");
            const auto& ptb                 = tr.getVar<data_t>("ptb");
            const auto& mtb                 = tr.getVar<data_t>("mtb");
            const auto& ISRJetPt            = tr.getVar<data_t>("ISRJetPt");
            
            float mtb_cut = 175.0;
            
            //------------------------------------------//
            //--- Updated Search Bins (January 2019) ---//
            //------------------------------------------//
            //================================SUS-16-049 (team_A) search bin low dm================================
            // int SB_team_A_lowdm(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
            //================================search bin v2 high dm================================================
            // int SBv2_highdm(float mtb_cut, float mtb, int njets, int ntop, int nw, int nres, int nb, float met, float ht)
            int nSearchBinLowDM  = SBv2_lowdm(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
            int nSearchBinHighDM = SBv2_highdm(mtb_cut, mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met, ht);
            
            //-------------------------------------------//
            //--- Updated Validation Bins (June 2019) ---//
            //-------------------------------------------//
            //================================SUS-16-049 (team_A) low dm validation=================================================
            // int SBv2_lowdm_validation(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
            //================================low dm validation high MET=================================================
            // int SBv2_lowdm_validation_high_MET(int nb, int nSV, float ISRpt, float met)
            //================================SBv2 high dm validation=================================================
            // int SBv2_highdm_validation(float mtb, int njets, int ntop, int nw, int nres, int nb, float met)
            int nValidationBinLowDM        = SBv2_lowdm_validation(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
            int nValidationBinLowDMHighMET = SBv2_lowdm_validation_high_MET(nBottoms, nSoftBottoms, ISRJetPt, met);
            int nValidationBinHighDM       = SBv2_highdm_validation(mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met); 
            
            tr.registerDerivedVar("nSearchBinLowDM",              nSearchBinLowDM);
            tr.registerDerivedVar("nSearchBinHighDM",             nSearchBinHighDM);
            tr.registerDerivedVar("nValidationBinLowDM",          nValidationBinLowDM);
            tr.registerDerivedVar("nValidationBinLowDMHighMET",   nValidationBinLowDMHighMET);
            tr.registerDerivedVar("nValidationBinHighDM",         nValidationBinHighDM);
        }

    public:

        GetSearchBin(){}
        
        ~GetSearchBin(){}

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };
}

#endif
