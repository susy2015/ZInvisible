#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "CountWLep.h"
#include "JetSort.h"
#include "RJet.h"
#include "TopPt.h"
#include "SusyAnaTools/Tools/customize.h"

void activateBranches(std::set<std::string>& activeBranches)
{
    for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
    for(auto& bn : AnaConsts::activatedBranchNames_DataOnly) activeBranches.insert(bn);
    activeBranches.insert("run");
}

////////////////////////////////

RegisterFunctionsNTuple::RegisterFunctionsNTuple(bool isCondor, const std::string &year, std::map<std::string, std::string> var_name_map) : RegisterFunctions(), year_(year)
{
    // Important: create objects!!
	countWLep    = new plotterFunctions::CountWLep;
	jetSort      = new plotterFunctions::JetSort;
	rJet         = new plotterFunctions::RJet(var_name_map);
	toppt        = new plotterFunctions::TopPt;

	sample_name_map = {
		{"QCD_Smear_HT_100to200_2016",   "QCD_HT_100to200_2016"},
		{"QCD_Smear_HT_200to300_2016",   "QCD_HT_200to300_2016"},
		{"QCD_Smear_HT_300to500_2016",   "QCD_HT_300to500_2016"},
		{"QCD_Smear_HT_500to700_2016",   "QCD_HT_500to700_2016"},
		{"QCD_Smear_HT_700to1000_2016",  "QCD_HT_700to1000_2016"},
		{"QCD_Smear_HT_1000to1500_2016", "QCD_HT_1000to1500_2016"},
		{"QCD_Smear_HT_1500to2000_2016", "QCD_HT_1500to2000_2016"},
		{"QCD_Smear_HT_2000toInf_2016",  "QCD_HT_2000toInf_2016"},
		{"QCD_Smear_HT_100to200_2017",   "QCD_HT_100to200_2017"},
		{"QCD_Smear_HT_200to300_2017",   "QCD_HT_200to300_2017"},
		{"QCD_Smear_HT_300to500_2017",   "QCD_HT_300to500_2017"},
		{"QCD_Smear_HT_500to700_2017",   "QCD_HT_500to700_2017"},
		{"QCD_Smear_HT_700to1000_2017",  "QCD_HT_700to1000_2017"},
		{"QCD_Smear_HT_1000to1500_2017", "QCD_HT_1000to1500_2017"},
		{"QCD_Smear_HT_1500to2000_2017", "QCD_HT_1500to2000_2017"},
		{"QCD_Smear_HT_2000toInf_2017",  "QCD_HT_2000toInf_2017"},
		{"QCD_Smear_HT_100to200_2018",   "QCD_HT_100to200_2018"},
		{"QCD_Smear_HT_200to300_2018",   "QCD_HT_200to300_2018"},
		{"QCD_Smear_HT_300to500_2018",   "QCD_HT_300to500_2018"},
		{"QCD_Smear_HT_500to700_2018",   "QCD_HT_500to700_2018"},
		{"QCD_Smear_HT_700to1000_2018",  "QCD_HT_700to1000_2018"},
		{"QCD_Smear_HT_1000to1500_2018", "QCD_HT_1000to1500_2018"},
		{"QCD_Smear_HT_1500to2000_2018", "QCD_HT_1500to2000_2018"},
		{"QCD_Smear_HT_2000toInf_2018",  "QCD_HT_2000toInf_2018"}};

}

void RegisterFunctionsNTuple::setSampleName(const std::string &sampleName)
{
	std::string sname;
	auto search = sample_name_map.find(sampleName);
	if(search != sample_name_map.end())
		sname = search->second;
	else
		sname = sampleName;
	topweightcalculator = std::make_unique<TopWeightCalculator>("TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_" + year_ + "_v1.0.6/tTagEff_" + year_ + ".root", sname, year_);
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
	if(countWLep)    delete countWLep;
	if(jetSort)      delete jetSort;
	if(rJet)         delete rJet;
	if(toppt)        delete toppt;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //register functions with NTupleReader
    // order matters
	tr.registerFunction(*countWLep);
	tr.registerFunction(*jetSort);
	tr.registerFunction(*rJet);
	tr.registerFunction(*toppt);
	if(topweightcalculator && tr.checkBranch("GenJet_pt")) tr.registerFunction(*topweightcalculator);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}
