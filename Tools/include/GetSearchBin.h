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
        std::string met_name_  = "MET_pt";
        std::map<std::string, std::string> prefixMap = {
            {"binNum",      "bin_"},
            {"unitCRNum",   "bin_lepcr_"},
            {"unitSRNum",   "bin_"},
        };

        void getSearchBin(NTupleReader& tr)
        {
            std::vector<std::string> JetPtCuts = {"_jetpt30"};
            // begin loop over jet pt cuts 
            for (const auto& suffix : JetPtCuts) 
            {
                // Note: Only HT, s_met, nJets, and dPhi are calculated with different jet pt cuts
                const auto& SAT_Pass_lowDM      = tr.getVar<bool>("SAT_Pass_lowDM" + suffix);
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
                
                // compare search bins (hui) vs. search bin units (matt) 
                //if (SAT_Pass_lowDM && ! (nSearchBinLowDM < 0 && nSearchbinUnitLowDM < 0))
                //{
                //    printf("nSB_hui = %d; nSB_matt = %d\n", nSearchBinLowDM, nSearchbinUnitLowDM);
                //}
                if (SAT_Pass_lowDM && (nSearchBinLowDM != nSearchbinUnitLowDM))
                {
                    printf("nSB_hui = %d; nSB_matt = %d\n", nSearchBinLowDM, nSearchbinUnitLowDM);
                }
                
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

        // See this link for cut definitions
        // https://github.com/mkilpatr/EstToolsSUSY/blob/SBv4/SUSYNano19/SRParameters_dc.hh#L122
        // 11 variables: 11 pass fuctions
        // also 1 function for using total number of top/W 
        // 12 pass functions in total
        bool pass_njets(const std::string& cut, int value)
        {
            if      (cut.compare("nj2to5")  == 0)   return bool(value >= 2 && value <= 5);
            else if (cut.compare("nj6")     == 0)   return bool(value >= 6);
            else if (cut.compare("nj7")     == 0)   return bool(value >= 7);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_nb(const std::string& cut, int value)
        {
            if      (cut.compare("nb0")     == 0)   return bool(value == 0);
            else if (cut.compare("nb1")     == 0)   return bool(value == 1);
            else if (cut.compare("nbgeq1")  == 0)   return bool(value >= 1);
            else if (cut.compare("nb2")     == 0)   return bool(value >= 2);
            else if (cut.compare("nbeq2")   == 0)   return bool(value == 2);
            else if (cut.compare("nb3")     == 0)   return bool(value >= 3);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_nsv(const std::string& cut, int value)
        {
            if      (cut.compare("nivf0")   == 0)   return bool(value == 0);
            else if (cut.compare("nivf1")   == 0)   return bool(value >= 1);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        // note: nrtntnw* and nrt* both start with nrt; be careful about this case
        //       check for nrtntnw* before nrt*
        bool pass_nTotalTopW(const std::string& cut, int value)
        {
            if      (cut.compare("nrtntnwgeq2")     == 0)   return bool(value >= 2);
            else if (cut.compare("nrtntnwgeq3")     == 0)   return bool(value >= 3);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_ntop(const std::string& cut, int value)
        {
            if      (cut.compare("nt0")     == 0)   return bool(value == 0);
            else if (cut.compare("nt1")     == 0)   return bool(value == 1);
            else if (cut.compare("nt2")     == 0)   return bool(value == 2);
            else if (cut.compare("ntgeq1")  == 0)   return bool(value >= 1);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_nw(const std::string& cut, int value)
        {
            if      (cut.compare("nw0")     == 0)   return bool(value == 0);
            else if (cut.compare("nw1")     == 0)   return bool(value == 1);
            else if (cut.compare("nw2")     == 0)   return bool(value == 2);
            else if (cut.compare("nwgeq1")  == 0)   return bool(value >= 1);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_nres(const std::string& cut, int value)
        {
            if      (cut.compare("nrt0")     == 0)   return bool(value == 0);
            else if (cut.compare("nrt1")     == 0)   return bool(value == 1);
            else if (cut.compare("nrt2")     == 0)   return bool(value == 2);
            else if (cut.compare("nrtgeq1")  == 0)   return bool(value >= 1);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_ISRpt(const std::string& cut, float value)
        {
            if      (cut.compare("lowptisr")     == 0)   return bool(value >= 300 && value < 500);
            else if (cut.compare("medptisr")     == 0)   return bool(value >= 300);
            else if (cut.compare("highptisr")    == 0)   return bool(value >= 500);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_mtb(const std::string& cut, float value)
        {
            if      (cut.compare("lowmtb")     == 0)   return bool(value <  175);
            else if (cut.compare("highmtb")    == 0)   return bool(value >= 175);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_ptb(const std::string& cut, float value)
        {
            if      (cut.compare("lowptb")       == 0)   return bool(value <  40);
            else if (cut.compare("medptb")       == 0)   return bool(value >= 40 && value < 70);
            else if (cut.compare("highptb")      == 0)   return bool(value >= 70);
            else if (cut.compare("lowptb12")     == 0)   return bool(value <  80);
            else if (cut.compare("medptb12")     == 0)   return bool(value >= 80 && value < 140);
            else if (cut.compare("highptb12")    == 0)   return bool(value >= 140);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_ht(const std::string& cut, float value)
        {
            if      (cut.compare("htlt1000")      == 0)   return bool(value <  1000);
            else if (cut.compare("htgt1000")      == 0)   return bool(value >= 1000);
            else if (cut.compare("ht1000to1500")  == 0)   return bool(value >= 1000 && value < 1500);
            else if (cut.compare("htgt1500")      == 0)   return bool(value >= 1500);
            else if (cut.compare("htlt1300")      == 0)   return bool(value <  1300);
            else if (cut.compare("htgt1300")      == 0)   return bool(value >= 1300);
            else if (cut.compare("ht1000to1300")  == 0)   return bool(value >= 1000 && value < 1300);
            else if (cut.compare("ht1300to1500")  == 0)   return bool(value >= 1300 && value < 1500);
            else    std::cout << "ERROR in " << __func__ << ": No string match found for " << cut << std::endl;
            return false;
        }
        bool pass_met(const std::string& cut, float value)
        {
            // example: MET_pt450to550, MET_pt550to650, MET_pt650to750, MET_pt750toinf 
            //std::string met_name  = "MET_pt";
            std::string separator = "to";
            // check that string begins with met name
            if (cut.find(met_name_) != 0)
            {
                std::cout << "ERROR in " << __func__ << ": The cut " << cut << " does not being with " << met_name_ << std::endl;
                return false;
            }
            int met_len = met_name_.length();
            int sep_len = separator.length();
            int sep_pos = cut.find(separator);
            int min_len = sep_pos - met_len;
            std::string min = cut.substr(met_len, min_len);
            std::string max = cut.substr(sep_pos + sep_len);
            //printf("%s: [%s, %s]\n", cut.c_str(), min.c_str(), max.c_str());
            // if max in inf, only apply min cut
            if (max.compare("inf") == 0)
            {
                float min_val = std::stoi(min);
                return bool(value >= min_val);
            }
            // otherwise, apply both min and max cuts
            else
            {
                float min_val = std::stoi(min);
                float max_val = std::stoi(max);
                return bool(value >= min_val && value < max_val);
            }
        }

        // return vector of strings of cuts from unit string 
        // split out beginning of unit name
        // note that MET_pt has '_' in name
        std::vector<std::string> getCutVec(const std::string& unit, const std::string& start)
        {
            std::vector<std::string> cuts;
            const char delim     = '_';
            //std::string met_name = "MET_pt";
            int start_len = start.length();
            int met_pos = unit.find(met_name_); 
            int final_len = met_pos - start_len - 1;
            std::string parsedUnit = unit.substr(start_len, final_len);
            std::string met_cut = unit.substr(met_pos);
            SusyUtility::splitString(parsedUnit, delim, cuts); 
            cuts.push_back(met_cut);
            return cuts;
        }
        
        // return true if event passes unit selection, otherwise return false 
        bool passUnitLowDM(const std::string& unit, const std::string& prefix, int njets, int nb, int nsv, float ISRpt, float ptb, float met)
        {
            //printf("%s: %s\n", __func__, unit.c_str());
            //printf("%s: ", unit.c_str());
            //std::string start = "bin_lm_";
            // check if unit is low dm
            if (unit.find(prefix) == 0)
            {
                std::vector<std::string> cuts = getCutVec(unit, prefix);
                //printf("%s: ", unit.c_str());
                for (const auto& c : cuts)
                {
                    // note; be careful about order as some cuts may begin with the same string
                    // optimization: if we do not pass a cut, return false
                    //printf("%s, ", c.c_str());
                    if (c.find("nj") == 0)
                    {
                        if (! pass_njets(c, njets))
                        {
                            return false; 
                        }
                    }
                    else if (c.find("nb") == 0)
                    {
                        if (! pass_nb(c, nb))
                        {
                            return false; 
                        }
                    }
                    else if (c.find("nivf") == 0)
                    {
                        if (! pass_nsv(c, nsv))
                        {
                            return false; 
                        }
                    }
                    else if (c.find("isr") != std::string::npos)
                    {
                        if (! pass_ISRpt(c, ISRpt))
                        {
                            return false; 
                        }
                    }
                    else if (c.find("ptb") != std::string::npos)
                    {
                        if (! pass_ptb(c, ptb))
                        {
                            return false; 
                        }
                    }
                    else if (c.find("MET_pt") == 0)
                    {
                        if (! pass_met(c, met))
                        {
                            return false; 
                        }
                    }
                    // skip lowmtb
                    // - lowmtb appears in some, but not all low dm search bins 
                    // - lowmtb does not need to be applied for low dm search bins as it is included in low dm baseline
                    else if (c.compare("lowmtb") == 0)
                    {
                        ;//null statement, similar to pass in python; no operation required
                    }
                    // if cut is not matched to any variable, print error and return false
                    else
                    {
                        std::cout << "ERROR in " << __func__ << ": No string match found for " << c << std::endl;
                        return false;
                    }
                }
                // if we reach the end then no cut is false; return true
                //printf("\n");
                return true;
            }
            // if unit is not low dm, return false
            else 
            {
                return false;
            }
        }
        
        // return true if event passes unit selection, otherwise return false 
        bool passUnitHighDM(const std::string& unit, float mtb, int njets, int nb, int ntop, int nw, int nres, float ht, float met)
        {
            return false;
        }
        
        // return unit number: can be used for search bins, CR units and SR units
        int getUnitNumLowDM(const std::string& key, int njets, int nb, int nsv, float ISRpt, float ptb, float met)
        {
            bool verbose = false;
            std::string prefix = prefixMap[key];
            prefix = prefix + "lm_";
            //printf("njets = %d, nb = %d, nsv = %d, ISRpt = %f, ptb = %f, met = %f\n", njets, nb, nsv, ISRpt, ptb, met);
            for (const auto& element : json_[key].items())
            {
                std::string unit = element.key();
                // only check units with low dm previx
                if (unit.find(prefix) == 0)
                {
                    bool pass = passUnitLowDM(unit, prefix, njets, nb, nsv, ISRpt, ptb, met);
                    //std::cout << unit << " : " << element.value() << std::endl;
                    printf("%s: pass = %s\n", unit.c_str(), pass ? "true" : "false");
                    if (pass)
                    {
                        int bin = std::stoi(std::string(element.value()));
                        if (verbose)
                        {
                            //printf("njets = %d, nb = %d, nsv = %d, ISRpt = %f, ptb = %f, met = %f\n", njets, nb, nsv, ISRpt, ptb, met);
                            //printf("event passes unit selection; %d : %s\n", bin, unit.c_str());
                            printf("pass selection for unit %d, %s; njets = %d, nb = %d, nsv = %d, ISRpt = %f, ptb = %f, met = %f\n", bin, unit.c_str(), njets, nb, nsv, ISRpt, ptb, met);
                        }
                        return bin;
                    }
                }
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
