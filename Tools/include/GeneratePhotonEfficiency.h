#ifndef GENERATEPHOTONEFFICIENCY_H
#define GENERATEPHOTONEFFICIENCY_H

#include "TypeDefinitions.h"
#include "PhotonTools.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"

#include "TopTagger.h"
#include "TTModule.h"
#include "TopTaggerUtilities.h"
#include "TopTaggerResults.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"

#include "TopTagger/TopTagger/include/TopObject.h"

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
    class GeneratePhotonEfficiency
    {
    private:
        TH1* hPhotonAccPt_num;
        TH1* hPhotonAccPt_den;
        TH1* hPhotonEffPt_num;
        TH1* hPhotonEffPt_den;

        void generatePhotonEfficiency(NTupleReader& tr)
        {
            float photonEfficienyPt;
            tr.registerDerivedVar("photonEfficienyPt", photonEfficienyPt);
        }
    
    public:
        GeneratePhotonEfficiency()
        {
            hPhotonAccPt_num = nullptr;
            hPhotonAccPt_den = nullptr;
            hPhotonEffPt_num = nullptr;
            hPhotonEffPt_den = nullptr;
            TH1::AddDirectory(false);
            std::string histFile = "effhists_GJets.root";
            TFile *f = new TFile(histFile);
            if(f)
            {
                hPhotonAccPt_num = static_cast<TH1*>(f->Get("hPhotonAccPt_num"));
                hPhotonAccPt_den = static_cast<TH1*>(f->Get("hPhotonAccPt_den"));
                hPhotonEffPt_num = static_cast<TH1*>(f->Get("hPhotonEffPt_num"));
                hPhotonEffPt_den = static_cast<TH1*>(f->Get("hPhotonEffPt_den"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open the file " << histFile << std::endl;
            }
        }
        ~GeneratePhotonEfficiency() {}

        void operator()(NTupleReader& tr)
        {
            generatePhotonEfficiency(tr);
        }
    };
}

#endif
