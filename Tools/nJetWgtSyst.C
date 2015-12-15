#include "../../SusyAnaTools/Tools/NTupleReader.h"

#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"

int main()
{
    TH1::AddDirectory(false);

    TH1* njWTTbar_0b;
    TH1* njWDYZ_0b;
    TH1* njWTTbar_g1b;
    TH1* njWDYZ_g1b;
    
    TRandom3 tr3(153474);

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

    std::cout << njWTTbar_0b << "\t" << njWDYZ_0b << "\t" << njWTTbar_g1b << "\t" << njWDYZ_g1b << std::endl;


    f = TFile::Open("condor/minituple.root");

    std::vector<std::string> treeNames_DY = {"DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
                                             "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
                                             "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
                                             "DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    };

    std::vector<std::string> treeNames_tt = {"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"};

    std::vector<std::string> treeNames_other = {};

    TCanvas c("c","c",800, 800);
    
    TLegend *leg = new TLegend(0.60, 0.60, 0.90, 0.90);

    TH1 *h[45];

    for(int i = 0; i < 45; ++i)
    {
        char name[128];
        sprintf(name, "hSB_%d", i);
        h[i] = new TH1D(name, name, 45, 0, 45);
    }

    //double wgt_0b_DY  = tr3.Gaus(1.0, rms_0b_DY/mean_0b_DY);
    //double wgt_g1b_DY = tr3.Gaus(1.0, rms_g1b_DY/mean_g1b_DY);
    //double wgt_0b_tt  = tr3.Gaus(1.0, rms_0b_tt/mean_0b_tt);
    //double wgt_g1b_tt = tr3.Gaus(1.0, rms_g1b_tt/mean_g1b_tt);

    for(auto& treeName : treeNames_DY)
    {
        TTree * t = (TTree*)f->Get(treeName.c_str());

        NTupleReader tr(t);

        while(tr.getNextEvent())
        {
            if(tr.getEvtNum() % 10000 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

            const int& nSearchBin = tr.getVar<int>("nSearchBin");
            const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");

            double mean_0b_DY  = njWDYZ_0b   ->GetBinContent(njWDYZ_0b->FindBin(cntNJetsPt30Eta24Zinv));
            double mean_g1b_DY = njWDYZ_g1b  ->GetBinContent(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
            double mean_0b_tt  = njWTTbar_0b ->GetBinContent(njWTTbar_0b->FindBin(cntNJetsPt30Eta24Zinv));
            double mean_g1b_tt = njWTTbar_g1b->GetBinContent(njWTTbar_g1b->FindBin(cntNJetsPt30Eta24Zinv));
            
            double rms_0b_DY  = njWDYZ_0b   ->GetBinError(njWDYZ_0b->FindBin(cntNJetsPt30Eta24Zinv));
            double rms_g1b_DY = njWDYZ_g1b  ->GetBinError(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
            double rms_0b_tt  = njWTTbar_0b ->GetBinError(njWTTbar_0b->FindBin(cntNJetsPt30Eta24Zinv));
            double rms_g1b_tt = njWTTbar_g1b->GetBinError(njWTTbar_g1b->FindBin(cntNJetsPt30Eta24Zinv));

            if(nSearchBin >= 0 && nSearchBin < 45)
            {
                for(int iTrial = 0; iTrial < 500; ++iTrial)
                {
                    h[nSearchBin]->Fill(tr3.Gaus(1.0, rms_g1b_DY/mean_g1b_DY));
                }
            }
        }
    }
}
