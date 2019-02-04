#ifndef METSMEAR 
#define METSMEAR 

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

#include "TopTagger/TopTagger/interface/TopObject.h"

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
    class MetSmear
    {
    private:

        TRandom3 *trand;

        data_t calcB(data_t x1, data_t x2, data_t e1, data_t e2)
        {
            return (x2*x2*(e1-1) - x1*x1*(e2-1)) / (x1*x2*(x2-x1));
        }

        data_t calcC(data_t x1, data_t x2, data_t e1, data_t e2)
        {
            return (x2*(e1-1) - x1*(e2-1)) / (x1*x2*(x1-x2));
        }

        data_t calcQuad(data_t x, data_t x1, data_t x2, data_t e1, data_t e2)
        {
            return 1 + x*calcB(x1, x2, e1, e2) + x*x*calcC(x1, x2, e1, e2);
        }

        data_t logistical(data_t met, data_t A, data_t B, data_t C)
        {
            return 1 - A/(1 + exp(-(B*(met - C))));
        }
        
        void metSmear(NTupleReader& tr)
        {
            const auto& met   = tr.getVar<data_t>("cleanMetPt");
            
            // Logistical smearing 
            data_t met_logi_1 = met * logistical(met, 0.15, 0.01, 300);
            data_t met_logi_2 = met * logistical(met, 0.20, 0.01, 400);
            data_t met_logi_3 = met * logistical(met, 0.25, 0.01, 500);
            data_t met_logi_4 = met * logistical(met, 0.20, 0.01, 400);
            data_t met_logi_5 = met * logistical(met, 0.20, 0.02, 400);
            data_t met_logi_6 = met * logistical(met, 0.20, 0.03, 400);
            data_t met_logi_7 = met * logistical(met, 0.20, 0.02, 300);
            data_t met_logi_8 = met * logistical(met, 0.20, 0.02, 400);
            data_t met_logi_9 = met * logistical(met, 0.20, 0.02, 500);

            // gaussian smearing 
            //data_t met_gaus_5  = trand->Gaus(met, 5);
            //data_t met_gaus_10 = trand->Gaus(met, 10);
            //data_t met_gaus_15 = trand->Gaus(met, 15);
            data_t met_gaus_20 = trand->Gaus(met, 20);
            //data_t met_gaus_25 = trand->Gaus(met, 25);
            data_t met_gaus_30 = trand->Gaus(met, 30);
            data_t met_gaus_40 = trand->Gaus(met, 40);
            data_t met_gaus_50 = trand->Gaus(met, 50);

            tr.registerDerivedVar("met_gaus_20", met_gaus_20);
            //tr.registerDerivedVar("met_gaus_25", met_gaus_25);
            tr.registerDerivedVar("met_gaus_30", met_gaus_30);
            tr.registerDerivedVar("met_gaus_40", met_gaus_40);
            tr.registerDerivedVar("met_gaus_50", met_gaus_50);

            tr.registerDerivedVar("met_logi_1", met_logi_1);
            tr.registerDerivedVar("met_logi_2", met_logi_2);
            tr.registerDerivedVar("met_logi_3", met_logi_3);
            tr.registerDerivedVar("met_logi_4", met_logi_4);
            tr.registerDerivedVar("met_logi_5", met_logi_5);
            tr.registerDerivedVar("met_logi_6", met_logi_6);
            tr.registerDerivedVar("met_logi_7", met_logi_7);
            tr.registerDerivedVar("met_logi_8", met_logi_8);
            tr.registerDerivedVar("met_logi_9", met_logi_9);
        }

        void mt2Smear(NTupleReader& tr)
        {
            const data_t& metphi       = tr.getVar<data_t>("cleanMetPhi");
            const data_t& met_logi_1   = tr.getVar<data_t>("met_logi_1");
            const data_t& met_gaus_30  = tr.getVar<data_t>("met_gaus_30");
            
            const std::vector<TLorentzVector>& jetsLVec_forTagger  = tr.getVec<TLorentzVector>("jetsLVec_forTaggerZinv");
            const std::vector<data_t>&     recoJetsBtag_forTagger  = tr.getVec<data_t>("recoJetsBtag_forTaggerZinv");

            //We choose 30 GeV gaussian smearing and logi_1 for the study

            // Form TLorentzVector of MET
            TLorentzVector metLVec_Logi;
            metLVec_Logi.SetPtEtaPhiM(met_logi_1, 0, metphi, 0);
            
            //type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Logi);
            data_t MT2_Logi = 0.0;//type3Ptr->best_had_brJet_MT2;

            TLorentzVector metLVec_Gaus;
            metLVec_Gaus.SetPtEtaPhiM(met_gaus_30, 0, metphi, 0);
            
            //type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Gaus); 
            data_t MT2_Gaus = 0.0;//type3Ptr->best_had_brJet_MT2;

            tr.registerDerivedVar("mt2_logi_1",  MT2_Logi);
            tr.registerDerivedVar("mt2_gaus_30", MT2_Gaus);
        }

    public:
        MetSmear()
        {
            trand = new TRandom3(452147);
        }

        //~MetSmear()
        // {
        //    if(trand) delete trand;
        // }

        void operator()(NTupleReader& tr)
        {
            metSmear(tr);
            mt2Smear(tr);
        }
    };
}

#endif
