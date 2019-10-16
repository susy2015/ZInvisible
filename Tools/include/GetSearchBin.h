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
#include <fstream>

// Repo for json.hpp: https://github.com/nlohmann/json/tree/master
#include "../../json/single_include/nlohmann/json.hpp"
using json = nlohmann::json;

namespace plotterFunctions
{
    class GetSearchBin
    {
    private:
        
        // json file containing unit mapping
        json json_;

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
                
                //------------------------------------------------//
                //--- Updated Search Bins: SBv4 (October 2019) ---//
                //------------------------------------------------//
                // int SBv4_lowdm(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
                // int SBv4_highdm(float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
                int nSearchBinLowDM  = SBv4_lowdm(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
                int nSearchBinHighDM = SBv4_highdm(mtb, nJets, nBottoms, nMergedTops, nWs, nResolvedTops, ht, met);
                
                //------------------------------------------------------//
                //--- Updated Validation Bins: SBv3 (September 2019) ---//
                //------------------------------------------------------//
                // int SBv3_lowdm_validation(int njets, int nb, int nSV, float ISRpt, float bottompt_scalar_sum, float met)
                // int SBv3_lowdm_validation_high_MET(int nb, int nSV, float ISRpt, float met)
                // int SBv3_highdm_validation(float mtb, int njets, int ntop, int nw, int nres, int nb, float met)
                int nValidationBinLowDM        = SBv3_lowdm_validation(nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
                int nValidationBinLowDMHighMET = SBv3_lowdm_validation_high_MET(nBottoms, nSoftBottoms, ISRJetPt, met);
                int nValidationBinHighDM       = SBv3_highdm_validation(mtb, nJets, nMergedTops, nWs, nResolvedTops, nBottoms, met); 

                //----------------------------------------//
                //--- Updated Unit Bins (October 2019) ---//
                //----------------------------------------//
                //int getUnitNumLowDM(const std::string& key, int njets, int nb, int nsv, float ISRpt, float ptb, float met)
                //int getUnitNumHighDM(const std::string& key, float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
                int nSearchbinUnitLowDM  = getUnitNumLowDM("binNum", nJets, nBottoms, nSoftBottoms, ISRJetPt, ptb, met);
                int nSearchbinUnitHighDM = getUnitNumHighDM("binNum", mtb, nJets, nBottoms, nMergedTops, nWs, nResolvedTops, ht, met);
                
                // search bins
                tr.registerDerivedVar("nSearchBinLowDM"             + suffix, nSearchBinLowDM);
                tr.registerDerivedVar("nSearchBinHighDM"            + suffix, nSearchBinHighDM);
                // validation bins
                tr.registerDerivedVar("nValidationBinLowDM"         + suffix, nValidationBinLowDM);
                tr.registerDerivedVar("nValidationBinLowDMHighMET"  + suffix, nValidationBinLowDMHighMET);
                tr.registerDerivedVar("nValidationBinHighDM"        + suffix, nValidationBinHighDM);
                // unit bins
                tr.registerDerivedVar("nSearchbinUnitLowDM"         + suffix, nSearchbinUnitLowDM);
                tr.registerDerivedVar("nSearchbinUnitHighDM"        + suffix, nSearchbinUnitHighDM);
            }
        }

    public:

        GetSearchBin()
        {
            bool print = false;
            const std::string fileName = "dc_BkgPred_BinMaps_master.json";
            loadJson(fileName);
            // print json file for testing
            if (print)
            {
                // void printJson(const std::string& fileName, const std::string& key, const std::string& title)
                SusyUtility::printJson(fileName, "binNum",      "Search Bins");
                SusyUtility::printJson(fileName, "unitCRNum",   "Control Region Units");
                SusyUtility::printJson(fileName, "unitSRNum",   "Search Region Units");
            }
        }
        
        ~GetSearchBin(){}
        
        // load json file
        void loadJson(const std::string& fileName)
        {
            // read json file
            std::ifstream i(fileName);
            i >> json_;
        }
        
        // return true if event passes unit selection, otherwise return false 
        bool passUnitLowDM(const std::string& unit, int njets, int nb, int nsv, float ISRpt, float ptb, float met)
        {
            bool pass = false;
            std::vector<std::string> cuts;
            const char delim = '_';
            std::string start = "bin_lm_";
            std::string met_name = "MET_pt";
            int start_len = start.length();
            int met_name_len = met_name.length();
            if (unit.find(start) == 0)
            {
                int met_pos = unit.find(met_name); 
                int final_len = met_pos - start_len - 1;
                std::string parsedUnit = unit.substr(start_len, final_len);
                std::string met_cut = unit.substr(met_pos);
                SusyUtility::splitString(parsedUnit, delim, cuts); 
                cuts.push_back(met_cut);
                //printf("unit: %s, %s, %s\n", unit.c_str(), parsedUnit.c_str(), met_cut.c_str());
                printf("%s: ", unit.c_str());
                for (const auto& c : cuts)
                {
                    printf("%s, ", c.c_str());
                }
                printf("\n");
            }
            return pass;
        }
        
        // return true if event passes unit selection, otherwise return false 
        bool passUnitHighDM(const std::string& unit, float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
        {
            return false;
        }
        
        // return unit number: can be used for search bins, CR units and SR units
        int getUnitNumLowDM(const std::string& key, int njets, int nb, int nsv, float ISRpt, float ptb, float met)
        {
            printf("njets = %d, nb = %d, nsv = %d, ISRpt = %f, pt = %f, met = %f\n", njets, nb, nsv, ISRpt, ptb, met);
            for (const auto& element : json_[key].items())
            {
                bool pass = passUnitLowDM(element.key(), njets, nb, nsv, ISRpt, ptb, met);
                if (pass)
                {
                    int bin = std::stoi(std::string(element.value()));
                    return bin;
                }
                //std::cout << element.key() << " : " << element.value() << std::endl;
            }
            return -1;
        }
        
        // return unit number: can be used for search bins, CR units and SR units
        int getUnitNumHighDM(const std::string& key, float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
        {
            for (const auto& element : json_[key].items())
            {
                bool pass = passUnitHighDM(element.key(), mtb, njets, nb, ntop, nw, nres, ht, met);
                if (pass)
                {
                    int bin = std::stoi(std::string(element.value()));
                    return bin;
                }
                //std::cout << element.key() << " : " << element.value() << std::endl;
            }
            return -1;
        }

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };
}

#endif
