#ifndef NJETAK8_H
#define NJETAK8_H

#include "TypeDefinitions.h"
#include "PhotonTools.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"
#include <vector>
#include <iostream>
#include <string>
#include <set>

namespace plotterFunctions
{
    class NJetAk8
    {
    private:
        void generateNJetAk8(NTupleReader& tr)
        {
            const auto& puppiJetsLVec  = tr.getVec<TLorentzVector>("puppiJetsLVec"); 
            // const int& nJetsAk8 = puppiJetsLVec.size();
            //tr.registerDerivedVar("nJetsAk8", nJetsAk8);
            //tr.registerDerivedVar("nJetsPuppi", nJetsPuppi);
        }
    public:
        NJetAk8(){}
        ~NJetAk8(){}
        void operator()(NTupleReader& tr)
        {
            generateNJetAk8(tr);
        }
    };   
}

#endif
