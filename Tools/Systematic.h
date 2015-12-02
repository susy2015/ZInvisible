#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "Plotter.h"
#include "SusyAnaTools/Tools/samples.h"
#include "NTupleReader.h"

#include "TF1.h"

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

#endif
