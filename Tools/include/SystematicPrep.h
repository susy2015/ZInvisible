#ifndef SYSTEMATICPREP_H
#define SYSTEMATICPREP_H

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
            const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            const std::vector<double>& recoJetsJecUnc    = tr.getVec<double>("recoJetsJecUncLepCleaned");

            const std::vector<double>& metMagUp   = tr.getVec<double>("metMagUp");
            const std::vector<double>& metMagDown = tr.getVec<double>("metMagDown");
            const std::vector<double>& metPhiUp   = tr.getVec<double>("metPhiUp");
            const std::vector<double>& metPhiDown = tr.getVec<double>("metPhiDown");

            const double& met    = tr.getVar<double>("met");
            const double& metphi = tr.getVar<double>("metphi");

            std::vector<TLorentzVector> *jetLVecUp = new std::vector<TLorentzVector>;
            std::vector<TLorentzVector> *jetLVecDn = new std::vector<TLorentzVector>;

            std::vector<double> *dPtMet = new std::vector<double>;
            std::vector<double> *dPhiMet = new std::vector<double>;

            double metUp = 0.0, metDn = 99990.0;

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
