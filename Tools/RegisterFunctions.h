#ifndef REGISTERFUNCTION_H
#define REGISTERFUNCTION_H

#include <set>
#include <string>
#include <vector>
#include <functional>

class NTupleReader;
class BaselineVessel;
class PDFUncertainty;
class BTagCorrector;

namespace plotterFunctions
{
    class GenerateWeight;
    class LepInfo;
    class Fakebtagvectors;
    class GetSearchBin;
    class TriggerInfo;
    class PrepareMiniTupleVars;
    class NJetWeight;
    class SystematicPrep;
    class SystematicCalc;
    class MetSmear;
    class PrepareTopVars;
}

class SystWeights;

class RegisterFunctions
{
private:

public:
    RegisterFunctions() {}
        
    virtual void registerFunctions(NTupleReader& tr) {};
    virtual void activateBranches(std::set<std::string>& activeBranches) {};
    virtual const std::set<std::string> getMiniTupleSet();
    virtual const std::set<std::string> getMiniTupleSetData();
    virtual void remakeBTagCorrector(std::string sampleName) {};
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
    BaselineVessel *myBLV;
    BaselineVessel *blvZinv;
    BaselineVessel *blvZinv1b;
    BaselineVessel *blvZinv2b;
    BaselineVessel *blvZinv3b;
    BaselineVessel *blvZinvJEUUp;
    BaselineVessel *blvZinvJEUDn;
    BaselineVessel *blvZinvMEUUp;
    BaselineVessel *blvZinvMEUDn;
    plotterFunctions::GenerateWeight *weights;
    plotterFunctions::NJetWeight *njWeight;
    plotterFunctions::LepInfo *lepInfo;
    plotterFunctions::Fakebtagvectors *fakebtagvectors;
    plotterFunctions::GetSearchBin *getSearchBin;
    plotterFunctions::TriggerInfo *triggerInfo;
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;
    plotterFunctions::SystematicPrep *systematicPrep;
    plotterFunctions::SystematicCalc *systematicCalc;
    PDFUncertainty *myPDFUnc;
    BTagCorrector *bTagCorrector;

public:
    RegisterFunctionsNTuple(bool isCondor = false, std::string sbEra = "SB_37_2015");
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
    void remakeBTagCorrector(std::string sampleName);
};

class RegisterFunctionsMiniTuple : public RegisterFunctions
{
private:
    plotterFunctions::NJetWeight *njWeight;
    plotterFunctions::TriggerInfo *triggerInfo;
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;
    plotterFunctions::MetSmear *metSmear;

public:
    RegisterFunctionsMiniTuple();
    ~RegisterFunctionsMiniTuple();
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

class RegisterFunctionsSyst : public RegisterFunctions
{
private:
    std::vector<std::function<void(NTupleReader&)> > funcs_;
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;
    plotterFunctions::NJetWeight *njWeight;
    plotterFunctions::TriggerInfo *triggerInfo;
    SystWeights *systWeights;
public:
    RegisterFunctionsSyst();
    ~RegisterFunctionsSyst();
    void addFunction(std::function<void(NTupleReader&)> func);
    void registerFunctions(NTupleReader& tr);
};

class RegisterFunctions2Dplot : public RegisterFunctions
{
private:
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;
public:
    RegisterFunctions2Dplot();
    ~RegisterFunctions2Dplot();
    void registerFunctions(NTupleReader& tr);
};

class RegisterFunctionsTopStudy : public RegisterFunctions
{
private:
    plotterFunctions::PrepareTopVars *prepareTopVars;
public:
    RegisterFunctionsTopStudy();
    ~RegisterFunctionsTopStudy();
    void registerFunctions(NTupleReader& tr);
};

void drawSBregionDefCopy(const double ymin_Yields = 0.05, const double ymax_Yields = 500., const bool logscale = true);

#endif
