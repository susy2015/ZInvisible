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
            const auto& nBottoms            = tr.getVar<int>("nBottoms");
            const auto& nSoftBottoms        = tr.getVar<int>("nSoftBottoms");
            const auto& nMergedTops         = tr.getVar<int>("nMergedTops");
            const auto& nJets               = tr.getVar<int>("nJets");
            const auto& nWs                 = tr.getVar<int>("nWs");
            const auto& nResolvedTops       = tr.getVar<int>("nResolvedTops");
            //const auto& met                 = tr.getVar<data_t>("cleanMetPt");
            const auto& met                 = tr.getVar<data_t>("metWithPhoton");
            const auto& ht                  = tr.getVar<data_t>("HT");
            const auto& bottompt_scalar_sum = tr.getVar<data_t>("ptb");
            const auto& mtb                 = tr.getVar<data_t>("mtb");
            //onst auto& softbLVec           = tr.getVec<TLorentzVector>("softbLVec");
            const auto& ISRJet              = tr.getVar<TLorentzVector>("ISRJet");
            
            //------------------------------------------//
            //--- Updated Search Bins (January 2019) ---//
            //------------------------------------------//

            //================================SUS-16-049 (team_A) search bin low dm================================
            //int SB_team_A_lowdm(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
            //================================search bin v2 high dm================================================
            //int SBv2_highdm(float mtb_cut, float mtb, int njets, int ntop, int nw, int nres, int nb, float met, float ht)
            
            //int nSV = softbLVec.size();
            float ISRpt = 0.0;
            ISRpt = ISRJet.Pt();
            //if(ISRJet.size() == 1) ISRpt = ISRJet.at(0).Pt();
            float mtb_cut = 175.0;

            int nSearchBinLowDM = SBv2_lowdm(nJets, nBottoms, nSoftBottoms, ISRpt, bottompt_scalar_sum, met);
                                                                    //n_top_merged, n_top_resolved
            int nSearchBinHighDM = SBv2_highdm(mtb_cut, mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met, ht);


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
