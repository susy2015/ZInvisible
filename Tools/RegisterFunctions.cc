#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "Gamma.h"
#include "BasicLepton.h"
#include "baselineDef.h"
#include "PDFUncertainty.h"
#include "BTagCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"

void activateBranches(std::set<std::string>& activeBranches)
{
    for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
    for(auto& bn : AnaConsts::activatedBranchNames_DataOnly) activeBranches.insert(bn);
    activeBranches.insert("run");
}

////////////////////////////////

RegisterFunctionsNTuple::RegisterFunctionsNTuple(bool isCondor, std::string year) : RegisterFunctions()
{            
    // Important: create objects!!
    myBLV               = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "");
    getVectors                   = new GetVectors;
    cleanedJets                  = new CleanedJets;
    runTopTagger                 = new RunTopTagger;
    gamma                        = new plotterFunctions::Gamma(year);
    basicLepton                  = new plotterFunctions::BasicLepton;
    myPDFUnc = new PDFUncertainty();
    bTagCorrector = nullptr;
    ISRcorrector = nullptr;
    pileup = nullptr;
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(getVectors)                   delete getVectors;
    if(cleanedJets)                  delete cleanedJets;
    if(runTopTagger)                 delete runTopTagger;
    if(myBLV)                        delete myBLV;
    if(basicLepton)                  delete basicLepton;
    if(myPDFUnc)                     delete myPDFUnc;
    if(bTagCorrector)                delete bTagCorrector;
    if(ISRcorrector)                 delete ISRcorrector;
    if(pileup)                       delete pileup;
    if(gamma)                        delete gamma;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    // order matters
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    tr.registerFunction(*cleanedJets);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*myBLV);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}

void RegisterFunctionsNTuple::remakeBTagCorrector(std::string sampleName)
{
    if(sampleName.find("Data") == std::string::npos){
        if(bTagCorrector) bTagCorrector->resetEffs(sampleName);
    }
}


void RegisterFunctionsNTuple::remakeISRreweight(std::string sampleName)
{
    if(sampleName.find("Data") == std::string::npos){
        if(ISRcorrector) ISRcorrector->resetSample(sampleName);
    }
}
