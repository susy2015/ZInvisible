#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "Systematic.h"
#include "PDFUncertainty.h"
#include "BTagCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"

const std::set<std::string> RegisterFunctions::getMiniTupleSet()
{
    //if you want to not fill the minituple return std::set<std::string>({});
  return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger","genHT","genWeight","bTagSF_EventWeightSimple_Central","isr_Unc_Cent","_PUweightFactor"});
    //return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger","genHT","genWeight"});
}

const std::set<std::string> RegisterFunctions::getMiniTupleSetData()
{
  return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger"});
}

void activateBranches(std::set<std::string>& activeBranches)
{
    activeBranches.insert("run");
}


////////////////////////////////

RegisterFunctionsNTuple::RegisterFunctionsNTuple(bool doSystematics, std::string sbEra, std::string year) : RegisterFunctions()
{            
    // Important: create objects!!
    myBLV                                       = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "");
    myBLV_jetpt30                               = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30");
    myBLV_jetpt30_jesTotalUp                    = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_jesTotalUp");
    myBLV_jetpt30_jesTotalDown                  = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_jesTotalDown");
    blv_drLeptonCleaned                         = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned");
    blv_drLeptonCleaned_jetpt30                 = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30");
    blv_drLeptonCleaned_jetpt30_jesTotalUp      = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30_jesTotalUp");
    blv_drLeptonCleaned_jetpt30_jesTotalDown    = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30_jesTotalDown");
    blv_drPhotonCleaned                         = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned");
    blv_drPhotonCleaned_jetpt30                 = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30");
    blv_drPhotonCleaned_jetpt30_jesTotalUp      = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30_jesTotalUp");
    blv_drPhotonCleaned_jetpt30_jesTotalDown    = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30_jesTotalDown");
    //blvZinv      = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "Zinv");
    //blvZinvJEUUp = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "ZinvJEUUp");
    //blvZinvJEUDn = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "ZinvJEUDn");
    //blvZinvMEUUp = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "ZinvMEUUp");
    //blvZinvMEUDn = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "ZinvMEUDn");

    std::string TopTaggerCfg     = "TopTagger_" + year +".cfg";
    getVectors                   = new GetVectors;
    cleanedJets                  = new CleanedJets;
    // RunTopTagger(std::string taggerCfg = "TopTagger.cfg", std::string suffix = "", bool doLeptonCleaning = false, bool doPhotonCleaning = false)
    runTopTagger                 = new RunTopTagger(TopTaggerCfg);
    runTopTagger_drLeptonCleaned = new RunTopTagger(TopTaggerCfg,"_drLeptonCleaned",true,false);
    runTopTagger_drPhotonCleaned = new RunTopTagger(TopTaggerCfg,"_drPhotonCleaned",false,true);
    gamma                        = new plotterFunctions::Gamma(year);
    weights                      = new plotterFunctions::GenerateWeight;
    generatePhotonEfficiency     = new plotterFunctions::GeneratePhotonEfficiency;
    njWeight                     = new plotterFunctions::NJetWeight;
    lepInfo                      = new plotterFunctions::LepInfo(year);
    basicLepton                  = new plotterFunctions::BasicLepton;
    getSearchBin                 = new plotterFunctions::GetSearchBin;
    prepareMiniTupleVars         = new plotterFunctions::PrepareMiniTupleVars(true);
    systematicPrep               = new plotterFunctions::SystematicPrep;
    systematicCalc               = new plotterFunctions::SystematicCalc(sbEra);
    shapeNJets                   = new plotterFunctions::ShapeNJets;

    doSystematics_ = doSystematics;
    myPDFUnc = new PDFUncertainty();
    bTagCorrector = nullptr;
    ISRcorrector = nullptr;
    pileup = nullptr;
    
}
RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(getVectors)                                  delete getVectors;
    if(cleanedJets)                                 delete cleanedJets;
    if(runTopTagger)                                delete runTopTagger;
    if(runTopTagger_drLeptonCleaned)                delete runTopTagger_drLeptonCleaned;
    if(runTopTagger_drPhotonCleaned)                delete runTopTagger_drPhotonCleaned;
    if(myBLV)                                       delete myBLV;
    if(myBLV_jetpt30)                               delete myBLV_jetpt30;
    if(myBLV_jetpt30_jesTotalUp)                    delete myBLV_jetpt30_jesTotalUp;
    if(myBLV_jetpt30_jesTotalDown)                  delete myBLV_jetpt30_jesTotalDown;
    if(blv_drLeptonCleaned)                         delete blv_drLeptonCleaned;
    if(blv_drLeptonCleaned_jetpt30)                 delete blv_drLeptonCleaned_jetpt30;
    if(blv_drLeptonCleaned_jetpt30_jesTotalUp)      delete blv_drLeptonCleaned_jetpt30_jesTotalUp;
    if(blv_drLeptonCleaned_jetpt30_jesTotalDown)    delete blv_drLeptonCleaned_jetpt30_jesTotalDown;
    if(blv_drPhotonCleaned)                         delete blv_drPhotonCleaned;
    if(blv_drPhotonCleaned_jetpt30)                 delete blv_drPhotonCleaned_jetpt30;
    if(blv_drPhotonCleaned_jetpt30_jesTotalUp)      delete blv_drPhotonCleaned_jetpt30_jesTotalUp;
    if(blv_drPhotonCleaned_jetpt30_jesTotalDown)    delete blv_drPhotonCleaned_jetpt30_jesTotalDown;
    //if(blvZinv)                    delete blvZinv;
    //if(blvZinv1b)                  delete blvZinv1b;
    //if(blvZinv2b)                  delete blvZinv2b;
    //if(blvZinv3b)                  delete blvZinv3b;
    //if(blvZinvJEUUp)               delete blvZinvJEUUp;
    //if(blvZinvJEUDn)               delete blvZinvJEUDn;
    //if(blvZinvMEUUp)               delete blvZinvMEUUp;
    //if(blvZinvMEUDn)               delete blvZinvMEUDn;
    if(weights)                      delete weights;
    if(generatePhotonEfficiency)     delete generatePhotonEfficiency;
    if(njWeight)                     delete njWeight;
    if(lepInfo)                      delete lepInfo;
    if(basicLepton)                  delete basicLepton;
    if(getSearchBin)                 delete getSearchBin;
    if(prepareMiniTupleVars)         delete prepareMiniTupleVars;
    if(myPDFUnc)                     delete myPDFUnc;
    if(systematicPrep)               delete systematicPrep;
    if(systematicCalc)               delete systematicCalc;
    if(bTagCorrector)                delete bTagCorrector;
    if(ISRcorrector)                 delete ISRcorrector;
    if(pileup)                       delete pileup;
    if(gamma)                        delete gamma;
    if(shapeNJets)                   delete shapeNJets;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    
    // order matters
    // get photons and leptons
    // use photons and leptons to clean jets
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    //tr.registerFunction(*generatePhotonEfficiency);
    tr.registerFunction(*cleanedJets);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*runTopTagger_drLeptonCleaned);
    tr.registerFunction(*runTopTagger_drPhotonCleaned);
    tr.registerFunction(*myBLV);
    tr.registerFunction(*myBLV_jetpt30);
    tr.registerFunction(*lepInfo);
    tr.registerFunction(*blv_drLeptonCleaned);
    tr.registerFunction(*blv_drLeptonCleaned_jetpt30);
    tr.registerFunction(*blv_drPhotonCleaned);
    tr.registerFunction(*blv_drPhotonCleaned_jetpt30);
    // apply JEC to MC only
    if (doSystematics_ && tr.checkBranch("GenJet_pt"))
    {
        tr.registerFunction(*myBLV_jetpt30_jesTotalUp);
        tr.registerFunction(*myBLV_jetpt30_jesTotalDown);
        tr.registerFunction(*blv_drLeptonCleaned_jetpt30_jesTotalUp);
        tr.registerFunction(*blv_drLeptonCleaned_jetpt30_jesTotalDown);
        tr.registerFunction(*blv_drPhotonCleaned_jetpt30_jesTotalUp);
        tr.registerFunction(*blv_drPhotonCleaned_jetpt30_jesTotalDown);
    }
    tr.registerFunction(*shapeNJets);
    tr.registerFunction(*getSearchBin);
    
    // old version including JEC and systematics
    // here for reference only
    //tr.registerFunction(*weights);
    //tr.registerFunction(*blvZinv);
    //tr.registerFunction(*blv_drLeptonCleaned);
    //tr.registerFunction(*blv_drPhotonCleaned);
    //tr.registerFunction(*njWeight);
    //tr.registerFunction(*blvZinv1b);
    //tr.registerFunction(*blvZinv2b);
    //tr.registerFunction(*blvZinv3b);
    //tr.registerFunction(*systematicPrep);
    //tr.registerFunction(*blvZinvJEUUp);
    //tr.registerFunction(*blvZinvJEUDn);
    //tr.registerFunction(*blvZinvMEUUp);
    //tr.registerFunction(*blvZinvMEUDn);
    //tr.registerFunction(*systematicCalc);
    //tr.registerFunction(*prepareMiniTupleVars);
    //tr.registerFunction(&printInterestingEvents);
    //tr.registerFunction(*myPDFUnc);
    //tr.registerFunction(*bTagCorrector);
    //tr.registerFunction(*ISRcorrector);
    //tr.registerFunction(*pileup);
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

////////////////////////////////

RegisterFunctionsMiniTuple::RegisterFunctionsMiniTuple() : RegisterFunctions()
{            
    njWeight             = new plotterFunctions::NJetWeight;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
}

RegisterFunctionsMiniTuple::~RegisterFunctionsMiniTuple()
{
    if(njWeight) delete njWeight;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
}
        
void RegisterFunctionsMiniTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    tr.registerFunction(*njWeight);
    tr.registerFunction(*prepareMiniTupleVars);
}

void RegisterFunctionsMiniTuple::activateBranches(std::set<std::string>& activeBranches)
{
}

////////////////////////////////

RegisterFunctionsCalcEff::RegisterFunctionsCalcEff() : RegisterFunctions()
{
    myBLV   = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "2016");
    blvZinv = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "2016", "Zinv");
    getVectors                = new GetVectors;
    cleanedJets               = new CleanedJets;
    runTopTagger              = new RunTopTagger;
    gamma                     = new plotterFunctions::Gamma;
    basicLepton               = new plotterFunctions::BasicLepton;
    generatePhotonEfficiency  = new plotterFunctions::GeneratePhotonEfficiency;
}

RegisterFunctionsCalcEff::~RegisterFunctionsCalcEff()
{
    if(myBLV) delete myBLV;
    if(blvZinv) delete blvZinv;
    if(getVectors) delete getVectors;
    if(cleanedJets) delete cleanedJets;
    if(runTopTagger) delete runTopTagger;
    if(gamma) delete gamma;
    if(basicLepton) delete basicLepton;
    if(generatePhotonEfficiency) delete generatePhotonEfficiency;
}

void RegisterFunctionsCalcEff::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    //tr.registerFunction(*generatePhotonEfficiency);
    tr.registerFunction(*cleanedJets);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*myBLV);
    tr.registerFunction(*blvZinv);

}

void RegisterFunctionsCalcEff::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}

////////////////////////////////

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

////////////////////////////////

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

////////////////////////////////

RegisterFunctionsTopStudy::RegisterFunctionsTopStudy()
{
    myBLV = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), "2016", "TopTag", "");
}

RegisterFunctionsTopStudy::~RegisterFunctionsTopStudy()
{
    if(myBLV) delete myBLV;
}

void RegisterFunctionsTopStudy::registerFunctions(NTupleReader& tr)
{
    tr.registerFunction(*myBLV);
}

/////////////////////////////////

void drawSBregionDefCopy(const double ymin_Yields, const double ymax_Yields, const bool logscale)
{
    SearchBins::drawSBregionDef(ymin_Yields, ymax_Yields, logscale);
}
