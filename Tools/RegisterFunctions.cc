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
    myBLV                                           = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "");
    myBLV_jetpt30                                   = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30");
    myBLV_jetpt30_jesTotalUp                        = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_jesTotalUp");
    myBLV_jetpt30_jesTotalDown                      = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_jesTotalDown");
    myBLV_jetpt30_METUnClustUp                      = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_METUnClustUp");
    myBLV_jetpt30_METUnClustDown                    = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_jetpt30_METUnClustDown");
    blv_drLeptonCleaned_jetpt30                     = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30");
    blv_drLeptonCleaned_jetpt30_jesTotalUp          = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30_jesTotalUp");
    blv_drLeptonCleaned_jetpt30_jesTotalDown        = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drLeptonCleaned_jetpt30_jesTotalDown");
    blv_drPhotonCleaned_jetpt30                     = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30");
    blv_drPhotonCleaned_jetpt30_jesTotalUp          = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30_jesTotalUp");
    blv_drPhotonCleaned_jetpt30_jesTotalDown        = new BaselineVessel(*static_cast<NTupleReader*>(nullptr), year, "_drPhotonCleaned_jetpt30_jesTotalDown");

    //std::string TopTaggerCfg                      = "TopTagger_" + year +".cfg";
    std::string TopTaggerCfg_Tensorflow             = "TopTagger_Tensorflow_"           + year + ".cfg";
    std::string TopTaggerCfg_DiscriminatorFilter    = "TopTagger_DiscriminatorFilter_"  + year + ".cfg";
    getVectors                                      = new GetVectors;
    cleanedJets                                     = new CleanedJets;
    // RunTopTagger(std::string taggerCfg = "TopTagger.cfg", std::string suffix = "", bool doLeptonCleaning = false, bool doPhotonCleaning = false)
    runTopTagger                                    = new RunTopTagger(TopTaggerCfg_Tensorflow,          "_jetpt30");
    runTopTagger_jesTotalUp                         = new RunTopTagger(TopTaggerCfg_Tensorflow,          "_jetpt30_jesTotalUp");
    runTopTagger_jesTotalDown                       = new RunTopTagger(TopTaggerCfg_Tensorflow,          "_jetpt30_jesTotalDown");
    //runTopTagger                                    = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_jetpt30");
    //runTopTagger_jesTotalUp                         = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_jetpt30_jesTotalUp");
    //runTopTagger_jesTotalDown                       = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_jetpt30_jesTotalDown");
    runTopTagger_drLeptonCleaned                    = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drLeptonCleaned_jetpt30",                true,  false);
    runTopTagger_drLeptonCleaned_jesTotalUp         = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drLeptonCleaned_jetpt30_jesTotalUp",     true,  false);
    runTopTagger_drLeptonCleaned_jesTotalDown       = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drLeptonCleaned_jetpt30_jesTotalDown",   true,  false);
    runTopTagger_drPhotonCleaned                    = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drPhotonCleaned_jetpt30",                false, true);
    runTopTagger_drPhotonCleaned_jesTotalUp         = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drPhotonCleaned_jetpt30_jesTotalUp",     false, true);
    runTopTagger_drPhotonCleaned_jesTotalDown       = new RunTopTagger(TopTaggerCfg_DiscriminatorFilter, "_drPhotonCleaned_jetpt30_jesTotalDown",   false, true);
    gamma                                           = new plotterFunctions::Gamma(year);
    weights                                         = new plotterFunctions::GenerateWeight;
    generatePhotonEfficiency                        = new plotterFunctions::GeneratePhotonEfficiency;
    njWeight                                        = new plotterFunctions::NJetWeight;
    lepInfo                                         = new plotterFunctions::LepInfo(year);
    basicLepton                                     = new plotterFunctions::BasicLepton;
    getSearchBin                                    = new plotterFunctions::GetSearchBin("_jetpt30");
    getSearchBin_jesTotalUp                         = new plotterFunctions::GetSearchBin("_jetpt30_jesTotalUp");
    getSearchBin_jesTotalDown                       = new plotterFunctions::GetSearchBin("_jetpt30_jesTotalDown");
    getSearchBin_METUnClustUp                       = new plotterFunctions::GetSearchBin("_jetpt30_METUnClustUp");
    getSearchBin_METUnClustDown                     = new plotterFunctions::GetSearchBin("_jetpt30_METUnClustDown");
    getSearchBin_drPhotonCleaned                    = new plotterFunctions::GetSearchBin("_drPhotonCleaned_jetpt30");
    getSearchBin_drPhotonCleaned_jesTotalUp         = new plotterFunctions::GetSearchBin("_drPhotonCleaned_jetpt30_jesTotalUp");
    getSearchBin_drPhotonCleaned_jesTotalDown       = new plotterFunctions::GetSearchBin("_drPhotonCleaned_jetpt30_jesTotalDown");
    prepareMiniTupleVars                            = new plotterFunctions::PrepareMiniTupleVars(true);
    systematicPrep                                  = new plotterFunctions::SystematicPrep;
    systematicCalc                                  = new plotterFunctions::SystematicCalc(sbEra);
    shapeWeights                                    = new plotterFunctions::ShapeWeights;

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
    if(runTopTagger_jesTotalUp)                     delete runTopTagger_jesTotalUp;
    if(runTopTagger_jesTotalDown)                   delete runTopTagger_jesTotalDown;
    if(runTopTagger_drLeptonCleaned)                delete runTopTagger_drLeptonCleaned;
    if(runTopTagger_drLeptonCleaned_jesTotalUp)     delete runTopTagger_drLeptonCleaned_jesTotalUp;
    if(runTopTagger_drLeptonCleaned_jesTotalDown)   delete runTopTagger_drLeptonCleaned_jesTotalDown;
    if(runTopTagger_drPhotonCleaned)                delete runTopTagger_drPhotonCleaned;
    if(runTopTagger_drPhotonCleaned_jesTotalUp)     delete runTopTagger_drPhotonCleaned_jesTotalUp;
    if(runTopTagger_drPhotonCleaned_jesTotalDown)   delete runTopTagger_drPhotonCleaned_jesTotalDown;
    if(myBLV)                                       delete myBLV;
    if(myBLV_jetpt30)                               delete myBLV_jetpt30;
    if(myBLV_jetpt30_jesTotalUp)                    delete myBLV_jetpt30_jesTotalUp;
    if(myBLV_jetpt30_jesTotalDown)                  delete myBLV_jetpt30_jesTotalDown;
    if(myBLV_jetpt30_METUnClustUp)                  delete myBLV_jetpt30_METUnClustUp;
    if(myBLV_jetpt30_METUnClustDown)                delete myBLV_jetpt30_METUnClustDown;
    if(blv_drLeptonCleaned_jetpt30)                 delete blv_drLeptonCleaned_jetpt30;
    if(blv_drLeptonCleaned_jetpt30_jesTotalUp)      delete blv_drLeptonCleaned_jetpt30_jesTotalUp;
    if(blv_drLeptonCleaned_jetpt30_jesTotalDown)    delete blv_drLeptonCleaned_jetpt30_jesTotalDown;
    if(blv_drPhotonCleaned_jetpt30)                 delete blv_drPhotonCleaned_jetpt30;
    if(blv_drPhotonCleaned_jetpt30_jesTotalUp)      delete blv_drPhotonCleaned_jetpt30_jesTotalUp;
    if(blv_drPhotonCleaned_jetpt30_jesTotalDown)    delete blv_drPhotonCleaned_jetpt30_jesTotalDown;
    if(weights)                                     delete weights;
    if(generatePhotonEfficiency)                    delete generatePhotonEfficiency;
    if(njWeight)                                    delete njWeight;
    if(lepInfo)                                     delete lepInfo;
    if(basicLepton)                                 delete basicLepton;
    if(getSearchBin)                                delete getSearchBin;
    if(getSearchBin_jesTotalUp)                     delete getSearchBin_jesTotalUp;
    if(getSearchBin_jesTotalDown)                   delete getSearchBin_jesTotalDown;
    if(getSearchBin_METUnClustUp)                   delete getSearchBin_METUnClustUp;
    if(getSearchBin_METUnClustDown)                 delete getSearchBin_METUnClustDown;
    if(getSearchBin_drPhotonCleaned)                delete getSearchBin_drPhotonCleaned;
    if(getSearchBin_drPhotonCleaned_jesTotalUp)     delete getSearchBin_drPhotonCleaned_jesTotalUp;
    if(getSearchBin_drPhotonCleaned_jesTotalDown)   delete getSearchBin_drPhotonCleaned_jesTotalDown;
    if(prepareMiniTupleVars)                        delete prepareMiniTupleVars;
    if(myPDFUnc)                                    delete myPDFUnc;
    if(systematicPrep)                              delete systematicPrep;
    if(systematicCalc)                              delete systematicCalc;
    if(bTagCorrector)                               delete bTagCorrector;
    if(ISRcorrector)                                delete ISRcorrector;
    if(pileup)                                      delete pileup;
    if(gamma)                                       delete gamma;
    if(shapeWeights)                                delete shapeWeights;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    // order matters
    // use photons and leptons to clean jets
    tr.registerFunction(*getVectors);
    tr.registerFunction(*gamma);
    tr.registerFunction(*basicLepton);
    tr.registerFunction(*cleanedJets);
    tr.registerFunction(*runTopTagger);
    tr.registerFunction(*runTopTagger_drLeptonCleaned);
    tr.registerFunction(*runTopTagger_drPhotonCleaned);
    tr.registerFunction(*myBLV_jetpt30);
    tr.registerFunction(*lepInfo);
    tr.registerFunction(*blv_drLeptonCleaned_jetpt30);
    tr.registerFunction(*blv_drPhotonCleaned_jetpt30);
    tr.registerFunction(*shapeWeights);
    // run getSearchBin for both Z nu nu MC (search region) and photon Data and MC (control region)
    tr.registerFunction(*getSearchBin);
    tr.registerFunction(*getSearchBin_drPhotonCleaned);
    
    if (doSystematics_ && tr.checkBranch("GenJet_pt"))
    {
        tr.registerFunction(*runTopTagger_jesTotalUp);
        tr.registerFunction(*runTopTagger_jesTotalDown);
        //tr.registerFunction(*runTopTagger_drLeptonCleaned_jesTotalUp);
        //tr.registerFunction(*runTopTagger_drLeptonCleaned_jesTotalDown);
        tr.registerFunction(*runTopTagger_drPhotonCleaned_jesTotalUp);
        tr.registerFunction(*runTopTagger_drPhotonCleaned_jesTotalDown);
        tr.registerFunction(*myBLV_jetpt30_jesTotalUp);
        tr.registerFunction(*myBLV_jetpt30_jesTotalDown);
        tr.registerFunction(*myBLV_jetpt30_METUnClustUp);
        tr.registerFunction(*myBLV_jetpt30_METUnClustDown);
        //tr.registerFunction(*blv_drLeptonCleaned_jetpt30_jesTotalUp);
        //tr.registerFunction(*blv_drLeptonCleaned_jetpt30_jesTotalDown);
        tr.registerFunction(*blv_drPhotonCleaned_jetpt30_jesTotalUp);
        tr.registerFunction(*blv_drPhotonCleaned_jetpt30_jesTotalDown);
        tr.registerFunction(*getSearchBin_jesTotalUp);
        tr.registerFunction(*getSearchBin_jesTotalDown);
        tr.registerFunction(*getSearchBin_METUnClustUp);
        tr.registerFunction(*getSearchBin_METUnClustDown);
        tr.registerFunction(*getSearchBin_drPhotonCleaned_jesTotalUp);
        tr.registerFunction(*getSearchBin_drPhotonCleaned_jesTotalDown);
    }
    
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
