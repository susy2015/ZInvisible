#ifndef SYSTEMATICPREP_H
#define SYSTEMATICPREP_H

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
    class SystematicPrep
    {
    private:

        void systematicPrep(NTupleReader& tr)
        {
            const auto& jetsLVec            = tr.getVec<TLorentzVector>("prodJetsNoLep_jetsLVec");
            const auto& recoJetsJecUnc      = tr.getVec<data_t>("recoJetsJecUncLepCleaned");
            const auto& metMagUp            = tr.getVec<data_t>("metMagUp");
            const auto& metMagDown          = tr.getVec<data_t>("metMagDown");
            const auto& metPhiUp            = tr.getVec<data_t>("metPhiUp");
            const auto& metPhiDown          = tr.getVec<data_t>("metPhiDown");
            const auto& met                 = tr.getVar<data_t>("met");
            const auto& metphi              = tr.getVar<data_t>("metphi");

            auto* jetLVecUp                 = new std::vector<TLorentzVector>;
            auto* jetLVecDn                 = new std::vector<TLorentzVector>;
            auto* dPtMet                    = new std::vector<data_t>;
            auto* dPhiMet                   = new std::vector<data_t>;

            data_t metUp = 0.0, metDn = 99990.0;

            for(int iMet = 0; iMet < metMagUp.size(); ++iMet)
            {
                metUp = std::max(metUp, metMagUp[iMet]);
                metDn = std::min(metDn, metMagDown[iMet]);
                
                dPtMet->push_back((metMagUp[iMet] - met)/met);
                dPtMet->push_back((metMagDown[iMet] - met)/met);
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiUp[iMet] - metphi));
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiDown[iMet] - metphi));
            }

            for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
            {
                jetLVecUp->push_back(jetsLVec[iJet] * (1 + recoJetsJecUnc[iJet]));
                jetLVecDn->push_back(jetsLVec[iJet] * (1 - recoJetsJecUnc[iJet]));
            }

            tr.registerDerivedVar("metMEUUp", metUp);
            tr.registerDerivedVar("metMEUDn", metDn);

            tr.registerDerivedVec("dPtMet", dPtMet);
            tr.registerDerivedVec("dPhiMet", dPhiMet);

            tr.registerDerivedVec("jetLVecUp", jetLVecUp);
            tr.registerDerivedVec("jetLVecDn", jetLVecDn);
        }

    public:
        SystematicPrep()
        {

        }

        void operator()(NTupleReader& tr)
        {
            systematicPrep(tr);
        }

    };
}

#endif
