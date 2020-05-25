#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "Gamma.h"
#include "GetSearchBin.h"
#include "CountWLep.h"
#include "DownsizeBootstrap.h"
#include "BasicLepton.h"
#include "JetSort.h"
#include "RJet.h"
#include "TopPt.h"

void activateBranches(std::set<std::string>& activeBranches)
{
    for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
    for(auto& bn : AnaConsts::activatedBranchNames_DataOnly) activeBranches.insert(bn);
    activeBranches.insert("run");
}

////////////////////////////////

RegisterFunctionsNTuple::RegisterFunctionsNTuple(bool isCondor, std::string year, std::map<std::string, std::string> var_name_map) : RegisterFunctions()
{            
    // Important: create objects!!
    getVectors   = new GetVectors;
    runTopTagger = new RunTopTagger;
    gamma        = new plotterFunctions::Gamma(year);
    basicLepton  = new plotterFunctions::BasicLepton;
    getSearchBin = new plotterFunctions::GetSearchBin;
	countWLep    = new plotterFunctions::CountWLep;
	downBoot     = new plotterFunctions::DownsizeBootstrap;
	jetSort      = new plotterFunctions::JetSort;
	rJet         = new plotterFunctions::RJet(var_name_map);
	toppt        = new plotterFunctions::TopPt;
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(getVectors)   delete getVectors;
    if(runTopTagger) delete runTopTagger;
    if(basicLepton)  delete basicLepton;
    if(getSearchBin) delete getSearchBin;
	if(countWLep)    delete countWLep;
	if(downBoot)     delete downBoot;
    if(gamma)        delete gamma;
	if(jetSort)      delete jetSort;
	if(rJet)         delete rJet;
	if(toppt)        delete toppt;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    // order matters
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*getSearchBin);
	tr.registerFunction(*countWLep);
	tr.registerFunction(*downBoot);
	tr.registerFunction(*jetSort);
	tr.registerFunction(*rJet);
	tr.registerFunction(*toppt);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}
