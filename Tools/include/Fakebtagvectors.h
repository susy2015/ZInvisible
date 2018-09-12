#ifndef FAKEBTAGVECTORS_H 
#define FAKEBTAGVECTORS_H 

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
    class Fakebtagvectors
    {
    private:
        void fakebtagvectors(NTupleReader& tr)
        {
            const auto& prodJetsNoLep_jetsLVec = tr.getVec<TLorentzVector>("prodJetsNoLep_jetsLVec");
            const auto& cleanJetpt30ArrBTag = tr.getVec<data_t>("prodJetsNoLep_recoJetsCSVv2");

            double maxCSV = 0.0;
            double secCSV = 0.0;
            double tenCSV = 0.0;
            int iMaxCSV = -1;
            int iSecCSV = -1;
            int iTenCSV = -1;

            if(prodJetsNoLep_jetsLVec.size() != cleanJetpt30ArrBTag.size()) std::cout << "fakebtagvectors(...): Vector size missmatch!!!!" << std::endl;

            int njet = 0;

            //find index of 3 highest CSV values
            for(int i = 0; i < cleanJetpt30ArrBTag.size(); ++i)
            {
                //Skip jets which cannot pass bTag Acceptance requirements
                if(!AnaFunctions::jetPassCuts(prodJetsNoLep_jetsLVec[i], AnaConsts::bTagArr)) continue;

                //count possible fake b-jets
                njet++;

                if(cleanJetpt30ArrBTag[i] > maxCSV)
                {
                    tenCSV = secCSV;
                    secCSV = maxCSV;
                    maxCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = iMaxCSV;
                    iMaxCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > secCSV)
                {
                    tenCSV = secCSV;
                    secCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > tenCSV)
                {
                    tenCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = i;
                }
            }

            auto* cleanJetpt30ArrBTag1fake = new std::vector<data_t>(cleanJetpt30ArrBTag);
            auto* cleanJetpt30ArrBTag2fake = new std::vector<data_t>(cleanJetpt30ArrBTag);
            auto* cleanJetpt30ArrBTag3fake = new std::vector<data_t>(cleanJetpt30ArrBTag);
            auto* fakedCSVValues = new std::vector<data_t>();

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag1fake)[iMaxCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iSecCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iSecCSV] = 0.99;
            if(iTenCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iTenCSV] = 0.99;

            if(iMaxCSV >= 0) fakedCSVValues->push_back(maxCSV);
            if(iSecCSV >= 0) fakedCSVValues->push_back(secCSV);
            if(iTenCSV >= 0) fakedCSVValues->push_back(tenCSV);

            //Calculate the combinatoric weights for b-jet faking
            data_t weight1fakeb = TMath::Binomial(njet, 1);
            data_t weight2fakeb = TMath::Binomial(njet, 2);
            data_t weight3fakeb = TMath::Binomial(njet, 3);
            //check for nans
            if(weight1fakeb != weight1fakeb) weight1fakeb = 0.0;
            if(weight2fakeb != weight2fakeb) weight2fakeb = 0.0;
            if(weight3fakeb != weight3fakeb) weight3fakeb = 0.0;

            tr.registerDerivedVar("weight1fakeb", weight1fakeb);
            tr.registerDerivedVar("weight2fakeb", weight2fakeb);
            tr.registerDerivedVar("weight3fakeb", weight3fakeb);

            tr.registerDerivedVec("cleanJetpt30ArrBTag1fake", cleanJetpt30ArrBTag1fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag2fake", cleanJetpt30ArrBTag2fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag3fake", cleanJetpt30ArrBTag3fake);
            tr.registerDerivedVec("fakedCSVValues", fakedCSVValues);
            tr.registerDerivedVar("maxCSV", maxCSV);
        }
    public:

        void operator()(NTupleReader& tr)
        {
            fakebtagvectors(tr);
        }
    };

}

#endif
