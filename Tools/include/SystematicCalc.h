#ifndef SYSTEMATICCALC_H
#define SYSTEMATICCALC_H

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
    class SystematicCalc
    {
    private:
        SearchBins sbins;

        void systematicCalc(NTupleReader& tr)
        {
            const int& cntCSVSJEUUp = tr.getVar<int>("cntCSVSZinvJEUUp");
            const int& nTopCandSortedCntJEUUp = tr.getVar<int>("nTopCandSortedCntZinvJEUUp");
            const double& MT2JEUUp = tr.getVar<double>("best_had_brJet_MT2ZinvJEUUp");

            const int& cntCSVSJEUDn = tr.getVar<int>("cntCSVSZinvJEUDn");
            const int& nTopCandSortedCntJEUDn = tr.getVar<int>("nTopCandSortedCntZinvJEUDn");
            const double& MT2JEUDn = tr.getVar<double>("best_had_brJet_MT2ZinvJEUDn");

            const double& cleanMet = tr.getVar<double>("cleanMetPt");

            const int& cntCSVSMEUUp = tr.getVar<int>("cntCSVSZinvMEUUp");
            const int& nTopCandSortedCntMEUUp = tr.getVar<int>("nTopCandSortedCntZinvMEUUp");
            const double& MT2MEUUp = tr.getVar<double>("best_had_brJet_MT2ZinvMEUUp");
            const double& cleanMetMEUUp = tr.getVar<double>("metMEUUp");

            const int& cntCSVSMEUDn = tr.getVar<int>("cntCSVSZinvMEUDn");
            const int& nTopCandSortedCntMEUDn = tr.getVar<int>("nTopCandSortedCntZinvMEUDn");
            const double& MT2MEUDn = tr.getVar<double>("best_had_brJet_MT2ZinvMEUDn");
            const double& cleanMetMEUDn = tr.getVar<double>("metMEUDn");

            const double& HTUp           = tr.getVar<double>("HTZinvJEUUp");
            const double& HTDn           = tr.getVar<double>("HTZinvJEUDn");
            const double& HTMEUUp           = tr.getVar<double>("HTZinvMEUUp");
            const double& HTMEUDn           = tr.getVar<double>("HTZinvMEUDn");

            int nSearchBinJEUUp = sbins.find_Binning_Index(cntCSVSJEUUp, nTopCandSortedCntJEUUp, MT2JEUUp, cleanMet, HTUp);
            int nSearchBinJEUDn = sbins.find_Binning_Index(cntCSVSJEUDn, nTopCandSortedCntJEUDn, MT2JEUDn, cleanMet, HTDn);

            int nSearchBinMEUUp = sbins.find_Binning_Index(cntCSVSMEUUp, nTopCandSortedCntMEUUp, MT2MEUUp, cleanMetMEUUp, HTMEUUp);
            int nSearchBinMEUDn = sbins.find_Binning_Index(cntCSVSMEUDn, nTopCandSortedCntMEUDn, MT2MEUDn, cleanMetMEUDn, HTMEUDn);
            
            tr.registerDerivedVar("nSearchBinJEUUp", nSearchBinJEUUp);
            tr.registerDerivedVar("nSearchBinJEUDn", nSearchBinJEUDn);

            tr.registerDerivedVar("nSearchBinMEUUp", nSearchBinMEUUp);
            tr.registerDerivedVar("nSearchBinMEUDn", nSearchBinMEUDn);
        }

    public:
        SystematicCalc(std::string sb_era) : sbins(sb_era)
        {
            
        }

        void operator()(NTupleReader& tr)
        {
            systematicCalc(tr);
        }

    };
}

#endif
