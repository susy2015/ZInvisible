#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"

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

    AnaSamples::SampleSet        ss("", 3.0);
    AnaSamples::SampleCollection sc(ss);

    f = TFile::Open("condor/minituple.root");

    TH1 *h[45];

    for(int i = 0; i < 45; ++i)
    {
        char name[128];
        sprintf(name, "hSB_%d", i);
        h[i] = new TH1D(name, name, 1000, 0, 2);
    }

    //for(auto& fs : sc["DYJetsToLL"])
    for(auto& fs : sc["ZJetsToNuNu"])
    {
        size_t start = fs.filePath.rfind('/');
        size_t stop  = fs.filePath.rfind('.');
        std::string treeName = fs.filePath.substr(start + 1, stop - start - 1);

        TTree * t = (TTree*)f->Get(treeName.c_str());

        plotterFunctions::PrepareMiniTupleVars pmt(false);

        NTupleReader tr(t);
        tr.registerFunction(pmt);

        while(tr.getNextEvent())
        {
            if(tr.getEvtNum() % 10000 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

            const int& nSearchBin = tr.getVar<int>("nSearchBin");
            const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
            const bool& passBaseline = tr.getVar<bool>("passBaseline");

            if(passBaseline)
            {
                double mean_g1b_DY = njWDYZ_g1b  ->GetBinContent(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
            
                double rms_g1b_DY = njWDYZ_g1b  ->GetBinError(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));

                if(nSearchBin >= 0 && nSearchBin < 45)
                {
                    for(int iTrial = 0; iTrial < 15000; ++iTrial)
                    {
                        h[nSearchBin]->Fill(tr3.Gaus(1.0, rms_g1b_DY/mean_g1b_DY), fs.getWeight());
                    }
                }
            }
        }
    }

    TH1 *syst68 = new TH1D("syst68", "syst68", 45, 0, 45);
    TH1 *systRMS = new TH1D("systRMS", "systRMS", 45, 0, 45);

    TFile fout("syst_nJetWgt.root", "RECREATE");
    for(int i = 0; i < 45; ++i) 
    {
        h[i]->Write();
        h[i]->Scale(1/h[i]->Integral(0, h[i]->GetNbinsX() + 1));
        double ll = -999.9, ul = -999.9;
        TH1* hint = (TH1*)h[i]->Clone((std::string(h[i]->GetName())+"_int").c_str());
        for(int iBin = 1; iBin <= h[i]->GetNbinsX(); ++iBin)
        {
            hint->SetBinContent(iBin, h[i]->Integral(0, iBin));
            if     (ll < 0 && h[i]->Integral(0, iBin) > 0.16) ll = h[i]->GetBinCenter(iBin);
            else if(ul < 0 && h[i]->Integral(0, iBin) > 0.84) ul = h[i]->GetBinCenter(iBin);
        }
        syst68->SetBinContent(i + 1, (ul - ll) / 2.0);
        systRMS->SetBinContent(i + 1, h[i]->GetRMS());
        //std::cout << "bin: " << i << "\tll: " << ll << "\tul: " << ul << "\tsym err: " << (ul - ll) / 2.0 << "\t:RMS: " << h[i]->GetRMS() << std::endl;
        hint->Write();
    }
    syst68->Write();
    systRMS->Write();
}
