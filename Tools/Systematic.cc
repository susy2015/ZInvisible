#include "Systematic.h"

void Systematic::bookHist(std::vector<Plotter::HistSummary>& vh, std::vector<AnaSamples::FileSummary>& vfs)
{
    Plotter::DatasetSummary dsDY_nunu_syst(  name_+"_"+var_+"_var", vfs, "", name_);
    Plotter::DatasetSummary dsDY_nunu_nosyst(name_+"_"+var_+"_Nom", vfs, "", "");
    Plotter::DataCollection dcDY_nunu( "single", "nSearchBin", {dsDY_nunu_nosyst, dsDY_nunu_syst});
    vh.emplace_back(PHS(name_+"_"+var_, {dcDY_nunu}, {2, 1}, "", 45, 0, 45, true, true, "Search Bin", "Events"));
}

void Systematic::modifyParameters(NTupleReader& tr)
{
    std::string type;
    tr.getType(var_, type);

    double val = -999.9;

    if(type.find("vector") != std::string::npos)
    {
        throw "Systematic::modifyParameters(...): Variables cannot be vectors. Var: " + var_; 
    }
    else
    {
        if     (type.find("double")       != std::string::npos) val = tr.getVar<double>(var_);
        else if(type.find("unsigned int") != std::string::npos) val = static_cast<double>(tr.getVar<unsigned int>(var_));
        else if(type.find("int")          != std::string::npos) val = static_cast<double>(tr.getVar<int>(var_));
        else if(type.find("float")        != std::string::npos) val = static_cast<double>(tr.getVar<float>(var_));
        else if(type.find("char")         != std::string::npos) val = static_cast<double>(tr.getVar<char>(var_));
        else if(type.find("short")        != std::string::npos) val = static_cast<double>(tr.getVar<short>(var_));
        else if(type.find("long")         != std::string::npos) val = static_cast<double>(tr.getVar<long>(var_));
        else
        {
            throw "Systematic::modifyParameters(...): variable not defined. Var: " + var_;
        }
    }

    double weight = func_->Eval(val);
    tr.registerDerivedVar(name_, weight);
}
