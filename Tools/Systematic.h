#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "Plotter.h"
#include "SusyAnaTools/Tools/samples.h"
#include "NTupleReader.h"

#include "TF1.h"
#include "TRandom3.h"

class Systematic
{
public:
    Systematic(const std::string name, const std::string var, TF1 *f) : name_(name), var_(var), func_(f) {};

    void bookHist(std::vector<Plotter::HistSummary>& vh, std::vector<AnaSamples::FileSummary>& vfs);

    void operator()(NTupleReader& tr)
    {
        modifyParameters(tr);
    }

private:
    std::string name_, var_;
    TF1 * const func_;

    void modifyParameters(NTupleReader& tr);
};

class SystWeights
{
private:
    TH1* njWTTbar_0b;
    TH1* njWDYZ_0b;
    TH1* njWTTbar_g1b;
    TH1* njWDYZ_g1b;

    TRandom3 *tr3;

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
