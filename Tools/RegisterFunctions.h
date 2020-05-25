#ifndef REGISTERFUNCTIONS_H
#define REGISTERFUNCTIONS_H

#include <set>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <memory>
#include "SusyAnaTools/Tools/TopWeightCalculator.h"

class NTupleReader;

namespace plotterFunctions
{
	class CountWLep;
	class JetSort;
	class RJet;
	class TopPt;
}

class SystWeights;

class RegisterFunctions
{
private:

public:
    RegisterFunctions() {}
        
    virtual void registerFunctions(NTupleReader& tr) {};
    virtual void activateBranches(std::set<std::string>& activeBranches) {};
	virtual void setSampleName(const std::string &sampleName) {};
    
};

class RegisterFunctionsNTuple : public RegisterFunctions
{
private:
	plotterFunctions::CountWLep    *countWLep;
	plotterFunctions::JetSort      *jetSort;
	plotterFunctions::RJet         *rJet;
	plotterFunctions::TopPt        *toppt;
	std::unique_ptr<TopWeightCalculator> topweightcalculator;
	std::string year_;
	std::map<std::string, std::string> sample_name_map;

public:
    RegisterFunctionsNTuple(bool isCondor, const std::string &year, std::map<std::string, std::string> var_name_map);
    ~RegisterFunctionsNTuple();
    void registerFunctions(NTupleReader& tr);
    void activateBranches(std::set<std::string>& activeBranches);
	void setSampleName(const std::string &sampleName);
};

#endif
