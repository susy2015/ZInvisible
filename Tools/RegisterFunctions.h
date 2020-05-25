#ifndef REGISTERFUNCTIONS_H
#define REGISTERFUNCTIONS_H

#include <set>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include "SusyAnaTools/Tools/GetVectors.h"
#include "SusyAnaTools/Tools/RunTopTagger.h"

class NTupleReader;

namespace plotterFunctions
{
    class BasicLepton;
    class GetSearchBin;
	class CountWLep;
	class DownsizeBootstrap;
    class Gamma;
	class JetSort;
	class RJet;
	class TopPt;
}

class SystWeights;

class RegisterFunctions
{
private:

public:
    RegisterFunctions() {}
        
    virtual void registerFunctions(NTupleReader& tr) {};
    virtual void activateBranches(std::set<std::string>& activeBranches) {};
    
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
    GetVectors                     *getVectors;
    RunTopTagger                   *runTopTagger;
    plotterFunctions::Gamma        *gamma;
    plotterFunctions::BasicLepton  *basicLepton;
    plotterFunctions::GetSearchBin *getSearchBin;
	plotterFunctions::CountWLep    *countWLep;
	plotterFunctions::DownsizeBootstrap * downBoot;
	plotterFunctions::JetSort      *jetSort;
	plotterFunctions::RJet         *rJet;
	plotterFunctions::TopPt        *toppt;

public:
    RegisterFunctionsNTuple(bool isCondor, std::string year, std::map<std::string, std::string> var_name_map);
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
};

#endif
