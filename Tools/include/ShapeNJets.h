#ifndef SHAPENJETS_H
#define SHAPENJETS_H

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
#include <sys/stat.h>

namespace plotterFunctions
{
    class ShapeNJets
    {
    private:
        std::map <std::string, TH1D*> hmap;
        bool file_exists;

        void shapeNJets(NTupleReader& tr)
        {
            // don't run module if file does not exist
            if (! file_exists) return;
            const auto& njets  = tr.getVar<int>("nJets_drLeptonCleaned");
            // check that histogram is exists
            
            for (auto const& h : hmap){
                float njetWeight = 1.0;
                if (h.second)
                {
                    int bin_n = h.second->GetXaxis()->FindBin(njets);
                    njetWeight = h.second->GetBinContent(bin_n); 
                    // std::cout<<"The key: "<<h.first<<std::endl;
                    // std::cout<<"bin number: " <<bin_n<<" ,njets: "<<njets <<" ,bin value: "<<njetWeight<<std::endl;
                } 
                if (njetWeight==0.0){
                    njetWeight = 1.0;
                    // std::cout<<"values changed from 0.0 to "<<njetWeight<<std::endl;
                }

                tr.registerDerivedVar("njetWeight_" + h.first,  njetWeight);
            }
        }
    
    public:
        ShapeNJets()
        {
            TH1::AddDirectory(false);
            std::string histFile = "shapes_njets.root";
            // check if file exists
            struct stat buffer;  
            file_exists = bool(stat(histFile.c_str(), &buffer) == 0);
            TFile *f = new TFile(histFile.c_str());
            if(file_exists && f)
            {
               hmap["Electron_HighDM"]  = static_cast<TH1D*>(f->Get("njets_shape_2016_Electron_HighDM_Loose"));
               hmap["Electron_LowDM"]   = static_cast<TH1D*>(f->Get("njets_shape_2016_Electron_LowDM_Loose"));
               hmap["Muon_HighDM"]      = static_cast<TH1D*>(f->Get("njets_shape_2016_Muon_HighDM_Loose"));
               hmap["Muon_LowDM"]       = static_cast<TH1D*>(f->Get("njets_shape_2016_Muon_LowDM_Loose"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open the file " << histFile << ". The ShapeNJets.h module will not be run."<< std::endl;
            }
        }
        ~ShapeNJets() {}

        void operator()(NTupleReader& tr)
        {
            shapeNJets(tr);
        }
    };
}

#endif
