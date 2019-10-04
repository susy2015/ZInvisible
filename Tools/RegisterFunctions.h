#ifndef REGISTERFUNCTIONS_H
#define REGISTERFUNCTIONS_H

#include <set>
#include <string>
#include <vector>
#include <functional>
#include "SusyAnaTools/Tools/GetVectors.h"
#include "SusyAnaTools/Tools/RunTopTagger.h"

class NTupleReader;
class BaselineVessel;

namespace plotterFunctions
{
    class BasicLepton;
    class GetSearchBin;
    class Gamma;
	class JetSort;
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
    BaselineVessel                 *myBLV;
    GetVectors                     *getVectors;
    RunTopTagger                   *runTopTagger;
    plotterFunctions::Gamma        *gamma;
    plotterFunctions::BasicLepton  *basicLepton;
    plotterFunctions::GetSearchBin *getSearchBin;
	plotterFunctions::JetSort      *jetSort;

public:
    RegisterFunctionsNTuple(bool isCondor, std::string year);
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
};

#endif
