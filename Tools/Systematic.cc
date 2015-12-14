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

SystWeights::SystWeights()
{
    TH1::AddDirectory(false);
    
    njWTTbar_0b  = nullptr;
    njWDYZ_0b    = nullptr;
    njWTTbar_g1b = nullptr;
    njWDYZ_g1b   = nullptr;
    
    TFile *f = new TFile("dataMCweights.root");
    if(f)
    {
        njWTTbar_0b  = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_0b_ht200_dphi")->Clone());
        njWDYZ_0b    = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_0b_ht200_dphi")->Clone());
        njWTTbar_g1b = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_g1b_ht200_dphi")->Clone());
        njWDYZ_g1b   = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_g1b_ht200_dphi")->Clone());
        f->Close();
        delete f;
    }
    else
    {
        std::cout << "Failed to open: dataMCweights.root" << std::endl;
    }

    tr3 = new TRandom3(187419);
}

SystWeights::~SystWeights()
{
    if(tr3) delete tr3;
}

void SystWeights::getWeights(NTupleReader& tr)
{
    const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");

    
}
