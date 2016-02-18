#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "Plotter.h"
#include "SusyAnaTools/Tools/samples.h"
#include "NTupleReader.h"

#include "TF1.h"
#include "TH2.h"
#include "TRandom3.h"

class Systematic
{
public:
    Systematic(const std::string name, const std::string var, TF1 *f) : name_(name), var_(var), func_(f), hist_(nullptr), hist2_(nullptr) {};
    Systematic(const std::string name, const std::string var, TH1 *h) : name_(name), var_(var), func_(nullptr), hist_{h}, hist2_(nullptr) {};
    Systematic(const std::string name, const std::string var, const std::string var2, TH2 *h) : name_(name), var_(var), var2_(var2), func_(nullptr), hist_(nullptr), hist2_{h} {};

    void bookHist(std::vector<Plotter::HistSummary>& vh, std::vector<AnaSamples::FileSummary>& vfs);

    void operator()(NTupleReader& tr)
    {
        modifyParameters(tr);
    }

private:
    std::string name_, var_, var2_;
    TF1 * const func_;
    TH1 * const hist_;
    TH2 * const hist2_;

    void modifyParameters(NTupleReader& tr);
};

class SystWeights
{
private:
    TH1* njWTTbar_0b;
    TH1* njWDYZ_0b;
    TH1* njWTTbar_g1b;
    TH1* njWDYZ_g1b;

    TRandom3 tr3;

    void getWeights(NTupleReader& tr);

public:
    SystWeights();
    ~SystWeights();

    void operator()(NTupleReader& tr)
    {
        getWeights(tr);
    }        
};

#endif
