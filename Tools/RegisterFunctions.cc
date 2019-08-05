#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "Gamma.h"
#include "GetSearchBin.h"
#include "BasicLepton.h"
#include "baselineDef.h"

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
    myBLV        = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "");
    getVectors   = new GetVectors;
    runTopTagger = new RunTopTagger;
    gamma        = new plotterFunctions::Gamma(year);
    basicLepton  = new plotterFunctions::BasicLepton;
    getSearchBin = new plotterFunctions::GetSearchBin;
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(getVectors)   delete getVectors;
    if(runTopTagger) delete runTopTagger;
    if(myBLV)        delete myBLV;
    if(basicLepton)  delete basicLepton;
    if(getSearchBin) delete getSearchBin;
    if(gamma)        delete gamma;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    // order matters
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*myBLV);
    tr.registerFunction(*getSearchBin);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}
