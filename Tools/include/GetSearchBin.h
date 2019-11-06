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

			double mtb_cut = 175.0;
            
            int nSearchBinLowDM  = SBv4_lowdm(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
            int nSearchBinHighDM = SBv4_highdm(mtb, nJets, nBottoms, nMergedTops, nWs, nResolvedTops, ht, met);
			std::vector<int> *nSearchBinHighDMLoose = new std::vector<int>;
			for(int bin : SBv2_highdm_loose_bin(mtb_cut, mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met, ht))
				nSearchBinHighDMLoose->push_back(bin);
            
            int nValidationBinLowDM        = SBv3_lowdm_validation(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
            int nValidationBinLowDMHighMET = SBv3_lowdm_validation_high_MET(nBottoms, nSoftBottoms, ISRJetPt, met);
            int nValidationBinHighDM       = SBv3_highdm_validation(mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met); 
            
            tr.registerDerivedVar("nSearchBinLowDM",              nSearchBinLowDM);
            tr.registerDerivedVar("nSearchBinHighDM",             nSearchBinHighDM);
			tr.registerDerivedVec("nSearchBinHighDMLoose",        nSearchBinHighDMLoose);
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
