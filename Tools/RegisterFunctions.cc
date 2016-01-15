#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "Systematic.h"
#include "PDFUncertainty.h"

const std::set<std::string> RegisterFunctions::getMiniTupleSet()
{
    return std::set<std::string>({"HTZinv","cleanMetPt","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts"});
}

void activateBranches(std::set<std::string>& activeBranches)
{
    for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
    for(auto& bn : AnaConsts::activatedBranchNames_DataOnly) activeBranches.insert(bn);
    activeBranches.insert("ht");
    activeBranches.insert("run");
    activeBranches.insert("lumi");
    activeBranches.insert("event");
    activeBranches.insert("mht");
    activeBranches.insert("mhtphi");
    activeBranches.insert("genDecayPdgIdVec");
    activeBranches.insert("genDecayLVec");
    activeBranches.insert("muonsLVec");
    activeBranches.insert("muonsRelIso");
    activeBranches.insert("muonsMiniIso");
    activeBranches.insert("W_emuVec");
    activeBranches.insert("muonsCharge");
    activeBranches.insert("muonsMtw");
    activeBranches.insert("met");
    activeBranches.insert("metphi");
    activeBranches.insert("jetsLVec");
    activeBranches.insert("elesLVec");
    activeBranches.insert("elesRelIso");
    activeBranches.insert("recoJetsBtag_0");
    activeBranches.insert("loose_isoTrks_mtw");
    activeBranches.insert("elesMtw");
    activeBranches.insert("loose_isoTrks_iso");
    activeBranches.insert("loose_isoTrksLVec");
    activeBranches.insert("muMatchedJetIdx");
    activeBranches.insert("eleMatchedJetIdx");
    activeBranches.insert("recoJetschargedEmEnergyFraction");
    activeBranches.insert("recoJetsneutralEmEnergyFraction");
    activeBranches.insert("recoJetschargedHadronEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetschargedEmEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetsneutralEmEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetschargedHadronEnergyFraction");
    activeBranches.insert("elesisEB");
    activeBranches.insert("elesMiniIso");
    activeBranches.insert("elesCharge");
}


RegisterFunctionsNTuple::RegisterFunctionsNTuple() : RegisterFunctions()
{            
    AnaFunctions::prepareTopTagger();

    myBLV     = new BaselineVessel("", "");
    blvZinv   = new BaselineVessel("Zinv");
    blvZinv1b = new BaselineVessel("Zinv1b");
    blvZinv2b = new BaselineVessel("Zinv2b");
    blvZinv3b = new BaselineVessel("Zinv3b");
    blvZinvJEUUp = new BaselineVessel("ZinvJEUUp");
    blvZinvJEUDn = new BaselineVessel("ZinvJEUDn");

    weights              = new plotterFunctions::GenerateWeight;
    njWeight             = new plotterFunctions::NJetWeight;
    lepInfo              = new plotterFunctions::LepInfo;
    fakebtagvectors      = new plotterFunctions::Fakebtagvectors;
    getSearchBin         = new plotterFunctions::GetSearchBin;
    triggerInfo          = new plotterFunctions::TriggerInfo;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(true);

    myPDFUnc = new PDFUncertainty();
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(myBLV) delete myBLV;
    if(blvZinv) delete blvZinv;
    if(blvZinv1b) delete blvZinv1b;
    if(blvZinv2b) delete blvZinv2b;
    if(blvZinv3b) delete blvZinv3b;
    if(blvZinvJEUUp) delete blvZinvJEUUp;
    if(blvZinvJEUDn) delete blvZinvJEUDn;
    if(weights) delete weights;
    if(njWeight) delete njWeight;
    if(lepInfo) delete lepInfo;
    if(fakebtagvectors) delete fakebtagvectors;
    if(getSearchBin) delete getSearchBin;
    if(triggerInfo) delete triggerInfo;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
    if(myPDFUnc) delete myPDFUnc;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*myBLV);
    tr.registerFunction(*lepInfo);
    tr.registerFunction(*weights);
    tr.registerFunction(*blvZinv);
    tr.registerFunction(*njWeight);
    tr.registerFunction(*fakebtagvectors);
    tr.registerFunction(*blvZinv1b);
    tr.registerFunction(*blvZinv2b);
    tr.registerFunction(*blvZinv3b);
    tr.registerFunction(*blvZinvJEUUp);
    tr.registerFunction(*blvZinvJEUDn);
    tr.registerFunction(*getSearchBin);
    tr.registerFunction(*triggerInfo);
    tr.registerFunction(*prepareMiniTupleVars);
    //tr.registerFunction(&printInterestingEvents);
    tr.registerFunction(*myPDFUnc);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}


RegisterFunctionsCalcEff::RegisterFunctionsCalcEff() : RegisterFunctions()
{
    AnaFunctions::prepareTopTagger();

    myBLV     = new BaselineVessel;
    lepInfo   = new plotterFunctions::LepInfo;
}

RegisterFunctionsCalcEff::~RegisterFunctionsCalcEff()
{
    if(myBLV) delete myBLV;
    if(lepInfo) delete lepInfo;
}

void RegisterFunctionsCalcEff::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*myBLV);
    tr.registerFunction(*lepInfo);
}

void RegisterFunctionsCalcEff::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}


void RegisterFunctionsSyst::addFunction(std::function<void(NTupleReader&)> func)
{
    funcs_.emplace_back(func);
}

RegisterFunctionsSyst::RegisterFunctionsSyst() : RegisterFunctions()
{
    njWeight             = new plotterFunctions::NJetWeight;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
    systWeights          = new SystWeights;
}

RegisterFunctionsSyst::~RegisterFunctionsSyst()
{
    if(njWeight) delete njWeight;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
    if(systWeights) delete systWeights;
}

void RegisterFunctionsSyst::registerFunctions(NTupleReader& tr)
{
    for(auto& f : funcs_) tr.registerFunction(f);

    tr.registerFunction(*njWeight);
    tr.registerFunction(*prepareMiniTupleVars);
    //tr.registerFunction(*systWeights);
}

RegisterFunctions2Dplot::RegisterFunctions2Dplot()
{
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
}

RegisterFunctions2Dplot::~RegisterFunctions2Dplot()
{
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
}

void RegisterFunctions2Dplot::registerFunctions(NTupleReader& tr)
{
    tr.registerFunction(*prepareMiniTupleVars);
}
