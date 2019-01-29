#ifndef GETSEARCHBIN_H
#define GETSEARCHBIN_H

#include "TypeDefinitions.h"
#include "PhotonTools.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/SB2018.h"
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
    class GetSearchBin
    {
    private:
        SearchBins sbins;

        void getSearchBin(NTupleReader& tr)
        {
            const auto& nb                  = tr.getVar<int>("cntCSVSZinv");
            const auto& nt                  = tr.getVar<int>("nTopCandSortedCntZinv");
            const auto& njets               = tr.getVar<int>("cntNJetsPt20Eta24Zinv");
            const auto& nWs                 = tr.getVar<int>("nWs");
            const auto& nResolvedTops       = tr.getVar<int>("nResolvedTops");
            const auto& met                 = tr.getVar<data_t>("cleanMetPt");
            const auto& ht                  = tr.getVar<data_t>("HTZinv");
            const auto& bottompt_scalar_sum = tr.getVar<data_t>("ptbZinv");
            const auto& mtb                 = tr.getVar<data_t>("mtbZinv");
            const auto& softbLVec           = tr.getVec<TLorentzVector>("softbLVecZinv");
            const auto& ISRLVec             = tr.getVec<TLorentzVector>("vISRJetZinv");
            
            //------------------------------------------//
            //--- Updated Search Bins (January 2019) ---//
            //------------------------------------------//

            //================================SUS-16-049 (team_A) search bin low dm================================
            //int SB_team_A_lowdm(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
            //================================search bin v2 high dm================================================
            //int SBv2_highdm(float mtb_cut, float mtb, int njets, int ntop, int nw, int nres, int nb, float met, float ht)
            
            int nSV = softbLVec.size();
            float ISRpt = 0.0;
            if(ISRLVec.size() == 1) ISRpt = ISRLVec.at(0).Pt();
            float mtb_cut = 175.0;

            int nSearchBinLowDM = SB_team_A_lowdm(njets, nb, nSV, ISRpt, bottompt_scalar_sum, met);
            int nSearchBinHighDM = SBv2_highdm(mtb_cut, mtb, njets, nt, nWs, nResolvedTops, nb, met, ht);


            tr.registerDerivedVar("nSearchBinLowDM", nSearchBinLowDM);
            tr.registerDerivedVar("nSearchBinHighDM", nSearchBinHighDM);
        }

    public:

        GetSearchBin(std::string sb_era) : sbins(sb_era) {}

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };
}

#endif
