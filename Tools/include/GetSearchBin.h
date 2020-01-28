#ifndef GETSEARCHBIN_H
#define GETSEARCHBIN_H

#include "TypeDefinitions.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
//#include "SusyAnaTools/Tools/SusyUtility.h"
#include "SusyAnaTools/Tools/SB2018.h"
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <fstream>

namespace plotterFunctions
{
	class GetSearchBin
	{
		private:

			void getSearchBin(NTupleReader& tr)
			{
				const auto& nJets               = tr.getVar<int>("Stop0l_nJets");
				const auto& nBottoms            = tr.getVar<int>("Stop0l_nbtags");
				const auto& nSoftBottoms        = tr.getVar<int>("Stop0l_nSoftb");
				const auto& nMergedTops         = tr.getVar<int>("Stop0l_nTop");
				const auto& nResolvedTops       = tr.getVar<int>("Stop0l_nResolved");
				const auto& nWs                 = tr.getVar<int>("Stop0l_nW");
				//const auto& ht                  = tr.getVar<data_t>("HT");
				const auto& ptb                 = tr.getVar<data_t>("Stop0l_Ptb");
				const auto& mtb                 = tr.getVar<data_t>("Stop0l_Mtb");
				const auto& ISRJetPt            = tr.getVar<data_t>("Stop0l_ISRJetPt");
				const auto& met                 = tr.getVar<data_t>("MET_pt");

				//------------------------------------------------------//
				//--- Updated Validation Bins: SBv3 (September 2019) ---//
				//------------------------------------------------------//
				int nValidationBinLowDM        = SBv3_lowdm_validation(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
				int nValidationBinLowDMHighMET = SBv3_lowdm_validation_high_MET(nBottoms, nSoftBottoms, ISRJetPt, met);
				int nValidationBinHighDM       = SBv3_highdm_validation(mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met); 

				// validation bins
				tr.registerDerivedVar("nValidationBinLowDM", nValidationBinLowDM);
				tr.registerDerivedVar("nValidationBinLowDMHighMET", nValidationBinLowDMHighMET);
				tr.registerDerivedVar("nValidationBinHighDM", nValidationBinHighDM);
			}

		public:

			GetSearchBin() {}
			~GetSearchBin(){}

			void operator()(NTupleReader& tr)
			{
				getSearchBin(tr);
			}
	};
}

#endif
