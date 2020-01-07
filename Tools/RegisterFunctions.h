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
    class GenerateWeight;
    class GeneratePhotonEfficiency;
    class LepInfo;
    class BasicLepton;
    class GetSearchBin;
    class PrepareMiniTupleVars;
    class NJetWeight;
    class SystematicPrep;
    class SystematicCalc;
    class PrepareTopVars;
    class Gamma;
    class ShapeNJets;
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
    virtual void remakeISRreweight(std::string sampleName) {};
    
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
    bool doSystematics_; 
    BaselineVessel                              *myBLV;
    BaselineVessel                              *myBLV_jetpt30;
    BaselineVessel                              *myBLV_jetpt30_jesTotalUp;
    BaselineVessel                              *myBLV_jetpt30_jesTotalDown;
    BaselineVessel                              *myBLV_jetpt30_METUnClustUp;
    BaselineVessel                              *myBLV_jetpt30_METUnClustDown;
    BaselineVessel                              *blv_drLeptonCleaned_jetpt30;
    BaselineVessel                              *blv_drLeptonCleaned_jetpt30_jesTotalUp;
    BaselineVessel                              *blv_drLeptonCleaned_jetpt30_jesTotalDown;
    BaselineVessel                              *blv_drPhotonCleaned_jetpt30;
    BaselineVessel                              *blv_drPhotonCleaned_jetpt30_jesTotalUp;
    BaselineVessel                              *blv_drPhotonCleaned_jetpt30_jesTotalDown;
    GetVectors                                  *getVectors;
    CleanedJets                                 *cleanedJets;
    RunTopTagger                                *runTopTagger;
    RunTopTagger                                *runTopTagger_jesTotalUp;
    RunTopTagger                                *runTopTagger_jesTotalDown;
    RunTopTagger                                *runTopTagger_drLeptonCleaned;
    RunTopTagger                                *runTopTagger_drLeptonCleaned_jesTotalUp;
    RunTopTagger                                *runTopTagger_drLeptonCleaned_jesTotalDown;
    RunTopTagger                                *runTopTagger_drPhotonCleaned;
    RunTopTagger                                *runTopTagger_drPhotonCleaned_jesTotalUp;
    RunTopTagger                                *runTopTagger_drPhotonCleaned_jesTotalDown;
    plotterFunctions::Gamma                     *gamma;
    PDFUncertainty                              *myPDFUnc;
    BTagCorrector                               *bTagCorrector;
    ISRCorrector                                *ISRcorrector;
    Pileup_Sys                                  *pileup;
    plotterFunctions::GenerateWeight            *weights;
    plotterFunctions::GeneratePhotonEfficiency  *generatePhotonEfficiency;
    plotterFunctions::NJetWeight                *njWeight;
    plotterFunctions::BasicLepton               *basicLepton;
    plotterFunctions::LepInfo                   *lepInfo;
    plotterFunctions::GetSearchBin              *getSearchBin;
    plotterFunctions::GetSearchBin              *getSearchBin_jesTotalUp;
    plotterFunctions::GetSearchBin              *getSearchBin_jesTotalDown;
    plotterFunctions::GetSearchBin              *getSearchBin_METUnClustUp;
    plotterFunctions::GetSearchBin              *getSearchBin_METUnClustDown;
    plotterFunctions::GetSearchBin              *getSearchBin_drPhotonCleaned;
    plotterFunctions::GetSearchBin              *getSearchBin_drPhotonCleaned_jesTotalUp;
    plotterFunctions::GetSearchBin              *getSearchBin_drPhotonCleaned_jesTotalDown;
    plotterFunctions::PrepareMiniTupleVars      *prepareMiniTupleVars;
    plotterFunctions::SystematicPrep            *systematicPrep;
    plotterFunctions::SystematicCalc            *systematicCalc;
    plotterFunctions::ShapeNJets                *shapeNJets;

public:
    RegisterFunctionsNTuple(bool doSystematics, std::string sbEra, std::string year);
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
    void remakeBTagCorrector(std::string sampleName);
    void remakeISRreweight(std::string sampleName);
};

class RegisterFunctionsMiniTuple : public RegisterFunctions
{
private:
    plotterFunctions::NJetWeight *njWeight;
    plotterFunctions::PrepareMiniTupleVars *prepareMiniTupleVars;

public:
    RegisterFunctionsMiniTuple();
    ~RegisterFunctionsMiniTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
};

class RegisterFunctionsCalcEff : public RegisterFunctions
{
private:
    BaselineVessel                              *myBLV;
    GetVectors                                  *getVectors;
    CleanedJets                                 *cleanedJets;
    RunTopTagger                                *runTopTagger;
    plotterFunctions::Gamma                     *gamma;
    plotterFunctions::BasicLepton               *basicLepton;
    plotterFunctions::GeneratePhotonEfficiency  *generatePhotonEfficiency;

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
    BaselineVessel *myBLV;
    plotterFunctions::PrepareTopVars *prepareTopVars;

public:
    RegisterFunctionsTopStudy();
    ~RegisterFunctionsTopStudy();
    void registerFunctions(NTupleReader& tr);
};

void drawSBregionDefCopy(const double ymin_Yields = 0.05, const double ymax_Yields = 500., const bool logscale = true);

#endif
