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

    TH1* shapeMET;
    TH1* shapeMT2;
    TH1* shapeNT;
    TH1* shapeNB;
    
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

    f = new TFile("syst_shape.root");
    if(f)
    {
        shapeMET = static_cast<TH1*>(f->Get("ShapeRatio_met")->Clone());
        shapeMT2 = static_cast<TH1*>(f->Get("ShapeRatio_mt2")->Clone());
        shapeNT  = static_cast<TH1*>(f->Get("ShapeRatio_nt")->Clone());
        shapeNB  = static_cast<TH1*>(f->Get("ShapeRatio_nb")->Clone());
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

    TH1 *h[5][45];
    std::vector<std::string> hnames = {"njet", "met", "mt2", "nt", "nb"};

    for(int ih = 0; ih < 5; ++ih)
    {
        for(int i = 0; i < 45; ++i)
        {
            char name[128];
            sprintf(name, "hSB_%s_%d", hnames[ih].c_str(), i);
            h[ih][i] = new TH1D(name, name, 1000, 0, 2);
        }
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
            const double& cleanMetPt = tr.getVar<double>("cleanMetPt");
            const double& best_had_brJet_MT2Zinv = tr.getVar<double>("best_had_brJet_MT2Zinv");
            const int& cntCSVSZinv            = tr.getVar<int>("cntCSVSZinv");
            const int& nTopCandSortedCntZinv  = tr.getVar<int>("nTopCandSortedCntZinv");
            const bool& passBaseline = tr.getVar<bool>("passBaseline");

            if(passBaseline)
            {
                double mean_g1b_DY = njWDYZ_g1b  ->GetBinContent(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
                double rms_g1b_DY = njWDYZ_g1b  ->GetBinError(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));

                double mean_met = shapeMET->GetBinContent(shapeMET->FindBin(cleanMetPt));
                double rms_met =  shapeMET->GetBinError(  shapeMET->FindBin(cleanMetPt));
                
                double mean_mt2 = shapeMT2->GetBinContent(shapeMT2->FindBin(best_had_brJet_MT2Zinv));
                double rms_mt2 =  shapeMT2->GetBinError(  shapeMT2->FindBin(best_had_brJet_MT2Zinv));
                
                double mean_nt = shapeNT->GetBinContent(shapeNT->FindBin(nTopCandSortedCntZinv));
                double rms_nt =  shapeNT->GetBinError(  shapeNT->FindBin(nTopCandSortedCntZinv));
                
                double mean_nb = shapeNB->GetBinContent(shapeNB->FindBin(cntCSVSZinv));
                double rms_nb =  shapeNB->GetBinError(  shapeNB->FindBin(cntCSVSZinv));

                if(nSearchBin >= 0 && nSearchBin < 45)
                {
                    for(int iTrial = 0; iTrial < 15000; ++iTrial)
                    {
                        h[0][nSearchBin]->Fill(tr3.Gaus(1.0, rms_g1b_DY/mean_g1b_DY), fs.getWeight());
                        h[1][nSearchBin]->Fill(tr3.Gaus(1.0, rms_met/mean_met),       fs.getWeight());
                        h[2][nSearchBin]->Fill(tr3.Gaus(1.0, rms_mt2/mean_mt2),       fs.getWeight());
                        h[3][nSearchBin]->Fill(tr3.Gaus(1.0, rms_nt/mean_nt),         fs.getWeight());
                        h[4][nSearchBin]->Fill(tr3.Gaus(1.0, rms_nb/mean_nb),         fs.getWeight());
                    }
                }
            }
        }
    }

    TFile fout("syst_nJetWgt.root", "RECREATE");
    TH1 *syst68Max = new TH1D("syst68Max", "syst68Max", 45, 0, 45);
    for(int ih = 0; ih < 5; ++ih)
    {
        char name[128];
        sprintf(name, "syst68_%s", hnames[ih].c_str());
        TH1 *syst68 = new TH1D(name, name, 45, 0, 45);
        sprintf(name, "systRMS_%s", hnames[ih].c_str());
        TH1 *systRMS = new TH1D(name, name, 45, 0, 45);

        for(int i = 0; i < 45; ++i) 
        {
            h[ih][i]->Write();
            h[ih][i]->Scale(1/h[ih][i]->Integral(0, h[ih][i]->GetNbinsX() + 1));
            double ll = -999.9, ul = -999.9;
            TH1* hint = (TH1*)h[ih][i]->Clone((std::string(h[ih][i]->GetName())+"_int").c_str());
            for(int iBin = 1; iBin <= h[ih][i]->GetNbinsX(); ++iBin)
            {
                hint->SetBinContent(iBin, h[ih][i]->Integral(0, iBin));
                if     (ll < 0 && h[ih][i]->Integral(0, iBin) > 0.16) ll = h[ih][i]->GetBinCenter(iBin);
                else if(ul < 0 && h[ih][i]->Integral(0, iBin) > 0.84) ul = h[ih][i]->GetBinCenter(iBin);
            }
            syst68->SetBinContent(i + 1, (ul - ll) / 2.0);
            syst68Max->SetBinContent(i + 1, std::max(syst68Max->GetBinContent(i + 1), (ul - ll) / 2.0));
            systRMS->SetBinContent(i + 1, h[ih][i]->GetRMS());
            //std::cout << "bin: " << i << "\tll: " << ll << "\tul: " << ul << "\tsym err: " << (ul - ll) / 2.0 << "\t:RMS: " << h[ih][i]->GetRMS() << std::endl;
            hint->Write();
        }
        syst68->Write();
        systRMS->Write();
    }
    syst68Max->Write();
}
