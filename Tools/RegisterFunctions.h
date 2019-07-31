#ifndef REGISTERFUNCTIONS_H
#define REGISTERFUNCTIONS_H

#include <set>
#include <string>
#include <vector>
#include <functional>
#include "SusyAnaTools/Tools/GetVectors.h"
#include "SusyAnaTools/Tools/CleanedJets.h"
#include "SusyAnaTools/Tools/RunTopTagger.h"

class NTupleReader;
class BaselineVessel;
class PDFUncertainty;
class BTagCorrector;
class ISRCorrector;
class Pileup_Sys; 

namespace plotterFunctions
{
    class BasicLepton;
    class PrepareTopVars;
    class Gamma;
}

class SystWeights;

class RegisterFunctions
{
private:

public:
    RegisterFunctions() {}
        
    virtual void registerFunctions(NTupleReader& tr) {};
    virtual void activateBranches(std::set<std::string>& activeBranches) {};
    virtual void remakeBTagCorrector(std::string sampleName) {};
    virtual void remakeISRreweight(std::string sampleName) {};
    
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
    BaselineVessel *myBLV;
    GetVectors                                  *getVectors;
    CleanedJets                                 *cleanedJets;
    RunTopTagger                                *runTopTagger;
    plotterFunctions::Gamma                     *gamma;
    PDFUncertainty                              *myPDFUnc;
    BTagCorrector                               *bTagCorrector;
    ISRCorrector                                *ISRcorrector;
    Pileup_Sys                                  *pileup;
    plotterFunctions::BasicLepton               *basicLepton;

public:
    RegisterFunctionsNTuple(bool isCondor, std::string year);
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
    void remakeBTagCorrector(std::string sampleName);
    void remakeISRreweight(std::string sampleName);
};

#endif
