#ifndef SYSTEMATICCALC_H
#define SYSTEMATICCALC_H

#include "TypeDefinitions.h"
#include "PhotonTools.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
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
    class SystematicCalc
    {
    private:
        SearchBins sbins;

        void systematicCalc(NTupleReader& tr)
        {
            const auto& cntCSVSJEUUp                = tr.getVar<int>("cntCSVSZinvJEUUp");
            const auto& nTopCandSortedCntJEUUp      = tr.getVar<int>("nTopCandSortedCntZinvJEUUp");
            const auto& MT2JEUUp                    = tr.getVar<data_t>("best_had_brJet_MT2ZinvJEUUp");
            const auto& cntCSVSJEUDn                = tr.getVar<int>("cntCSVSZinvJEUDn");
            const auto& nTopCandSortedCntJEUDn      = tr.getVar<int>("nTopCandSortedCntZinvJEUDn");
            const auto& MT2JEUDn                    = tr.getVar<data_t>("best_had_brJet_MT2ZinvJEUDn");
            const auto& cleanMet                    = tr.getVar<data_t>("cleanMetPt");
            const auto& cntCSVSMEUUp                = tr.getVar<int>("cntCSVSZinvMEUUp");
            const auto& nTopCandSortedCntMEUUp      = tr.getVar<int>("nTopCandSortedCntZinvMEUUp");
            const auto& MT2MEUUp                    = tr.getVar<data_t>("best_had_brJet_MT2ZinvMEUUp");
            const auto& cleanMetMEUUp               = tr.getVar<data_t>("metMEUUp");
            const auto& cntCSVSMEUDn                = tr.getVar<int>("cntCSVSZinvMEUDn");
            const auto& nTopCandSortedCntMEUDn      = tr.getVar<int>("nTopCandSortedCntZinvMEUDn");
            const auto& MT2MEUDn                    = tr.getVar<data_t>("best_had_brJet_MT2ZinvMEUDn");
            const auto& cleanMetMEUDn               = tr.getVar<data_t>("metMEUDn");
            const auto& HTUp                        = tr.getVar<data_t>("HTZinvJEUUp");
            const auto& HTDn                        = tr.getVar<data_t>("HTZinvJEUDn");
            const auto& HTMEUUp                     = tr.getVar<data_t>("HTZinvMEUUp");
            const auto& HTMEUDn                     = tr.getVar<data_t>("HTZinvMEUDn");

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
