#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"

const std::set<std::string> RegisterFunctions::getMiniTupleSet()
{
    return std::set<std::string>({"HTZinv","cleanMetPt","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv"});
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

    myBLV     = new BaselineVessel("", "SingleMuon_csc2015.txt");
    blvZinv   = new BaselineVessel("Zinv");
    blvZinv1b = new BaselineVessel("Zinv1b");
    blvZinv2b = new BaselineVessel("Zinv2b");
    blvZinv3b = new BaselineVessel("Zinv3b");

    weights              = new plotterFunctions::GenerateWeight;
    lepInfo              = new plotterFunctions::LepInfo;
    fakebtagvectors      = new plotterFunctions::Fakebtagvectors;
    getSearchBin         = new plotterFunctions::GetSearchBin;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars;
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(myBLV) delete myBLV;
    if(blvZinv) delete blvZinv;
    if(blvZinv1b) delete blvZinv1b;
    if(blvZinv2b) delete blvZinv2b;
    if(blvZinv3b) delete blvZinv3b;
    if(weights) delete weights;
    if(lepInfo) delete lepInfo;
    if(fakebtagvectors) delete fakebtagvectors;
    if(getSearchBin) delete getSearchBin;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*myBLV);
    tr.registerFunction(*lepInfo);
    tr.registerFunction(*weights);
    tr.registerFunction(*blvZinv);
    tr.registerFunction(*fakebtagvectors);
    tr.registerFunction(*blvZinv1b);
    tr.registerFunction(*blvZinv2b);
    tr.registerFunction(*blvZinv3b);
    tr.registerFunction(*getSearchBin);
    tr.registerFunction(*prepareMiniTupleVars);
    //tr.registerFunction(&printInterestingEvents);
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
