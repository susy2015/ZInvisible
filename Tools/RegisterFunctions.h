#ifndef REGISTERFUNCTION_H
#define REGISTERFUNCTION_H

#include <set>
#include <string>

class NTupleReader;
class BaselineVessel;

namespace plotterFunctions
{
    class GenerateWeight;
    class LepInfo;
    class Fakebtagvectors;
    class GetSearchBin;
    class TriggerInfo;
    class PrepareMiniTupleVars;
}

class RegisterFunctions
{
private:

public:
    RegisterFunctions() {}
        
    virtual void registerFunctions(NTupleReader& tr) {};
    virtual void activateBranches(std::set<std::string>& activeBranches) {};
    virtual const std::set<std::string> getMiniTupleSet();
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
    BaselineVessel *myBLV;
    BaselineVessel *blvZinv;
    BaselineVessel *blvZinv1b;
    BaselineVessel *blvZinv2b;
    BaselineVessel *blvZinv3b;
    plotterFunctions::GenerateWeight *weights;
    plotterFunctions::LepInfo *lepInfo;
    plotterFunctions::Fakebtagvectors *fakebtagvectors;
    plotterFunctions::GetSearchBin *getSearchBin;
    plotterFunctions::TriggerInfo *triggerInfo;
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;

public:
    RegisterFunctionsNTuple();
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
};

class RegisterFunctionsCalcEff : public RegisterFunctions
{
private:
    BaselineVessel *myBLV;
    plotterFunctions::LepInfo *lepInfo;

public:
    RegisterFunctionsCalcEff();
    ~RegisterFunctionsCalcEff();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
};

#endif
