#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <cstdio>
#include <string>

TH2* rebin2d(TH2* old, double xbins[], int nx, double ybins[], int ny)
{
    //create a new TH2 with your bin arrays spec
    //int nx = sizeof(xbins)/sizeof(double) - 1;

    //int ny = sizeof(ybins)/sizeof(double) - 1;

    char hname[128];
    sprintf(hname, "%s_rebin", old->GetName());
    TH2 *h = new TH2F(hname, old->GetTitle(), nx, xbins, ny, ybins);
    TAxis *xaxis = old->GetXaxis();
    TAxis *yaxis = old->GetYaxis();
    for (int i=1; i<=xaxis->GetNbins(); i++) 
    {
        for(int j=1; j<=yaxis->GetNbins(); j++) 
        {
            h->Fill(xaxis->GetBinCenter(i), yaxis->GetBinCenter(j), old->GetBinContent(i, j));
        }
    }

    return h;
}

static double dzbt[] = {0.0};

void makeRatio1D(std::string label, TFile *fin, TFile *fout, double bins[] = dzbt, int nbins = 0)
{
    fin->cd();
    TH1 *num = (TH1*)fin->Get((label+"_num").c_str());
    TH1 *den = (TH1*)fin->Get((label+"_den").c_str());

    if(nbins > 1)
    {
        num = num->Rebin(nbins - 1, "", bins);
        den = den->Rebin(nbins - 1, "", bins);
    }

    TH1 *ratio = (TH1*)num->Clone((label + "_ratio").c_str());
    ratio->Divide(den);

    fout->cd();
    ratio->Write();

    std::cout << "HELLO!" << std::endl;
}

void makeRatio2D(std::string label, TFile *fin, TFile *fout, double binsX[] = dzbt, int nx = 0, double binsY[] = dzbt, int ny = 0)
{
    fin->cd();
    TH2 *num = (TH2*)fin->Get((label+"_num").c_str());
    TH2 *den = (TH2*)fin->Get((label+"_den").c_str());

    if(nx > 1 && ny > 1)
    {
        num = rebin2d(num, binsX, nx - 1, binsY, ny - 1);
        den = rebin2d(den, binsX, nx - 1, binsY, ny - 1);
    }

    TH2 *ratio = (TH2*)num->Clone((label + "_ratio").c_str());
    ratio->Divide(den);

    fout->cd();
    ratio->Write();
    std::cout << label << "\tTHERE!" << std::endl;
}

int main ()
{
    TH1::AddDirectory(false);

    double xbins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 340.0, 380.0, 420.0, 460.0, 500.0, 560.0, 2000.0};
    int nxbins = sizeof(xbins)/sizeof(double);
    double ybins[] = {0.0, 10.0, 40.0, 120.0, 200.0, 300.0, 400.0, 650.0, 800.0, 1000.0, 5000.0};
    int nybins = sizeof(ybins)/sizeof(double);

    double zptbins[] = {0.0, 100.0, 140.0, 180.0, 220.0, 300.0, 400.0, 500.0, 600.0, 800.0, 2000.0};
    int nzptbins = sizeof(zptbins)/sizeof(double);

    double muptbins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 180.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 2000.0};
    int nmuptbins = sizeof(muptbins)/sizeof(double);

    double muptbins2[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 2000.0};
    int nmuptbins2 = sizeof(muptbins2)/sizeof(double);

    double muptbins3[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 180.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1200.0, 2000.0};
    int nmuptbins3 = sizeof(muptbins3)/sizeof(double);

    double atbins[] = {0.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 500.0, 3000.0};
    int nactbins = sizeof(atbins)/sizeof(double);

    TFile  *fin = new TFile("effhists.root");
    TFile *fout = new TFile("muEffHists2.root", "RECREATE");

    makeRatio1D("hMuEffPt", fin, fout, muptbins, nmuptbins);
    makeRatio1D("hMuAccPt", fin, fout);
    
    makeRatio1D("hMuEffHt", fin, fout);
    makeRatio1D("hMuAccHt", fin, fout);
    
    makeRatio2D("hMuEff", fin, fout);
    makeRatio2D("hMuAcc", fin, fout);
    
    makeRatio1D("hZEffPt", fin, fout, zptbins, nzptbins);
    makeRatio1D("hZAccPt", fin, fout, zptbins, nzptbins);
    
    makeRatio1D("hZAccPtSmear", fin, fout);//, zptbins, nzptbins);
    
    makeRatio2D("hZEff", fin, fout);
    makeRatio2D("hZAcc", fin, fout);
    
    makeRatio2D("hZEff_jActR1", fin, fout, zptbins, nzptbins, ybins, nybins);
    makeRatio2D("hZEff_jActR2", fin, fout);
    
    makeRatio2D("hMuEffPtActReco", fin, fout, muptbins2, nmuptbins2, atbins, nactbins);
    makeRatio2D("hMuEffPtActIso", fin, fout, muptbins2, nmuptbins2, atbins, nactbins);
    
    makeRatio1D("hMuEffPtReco", fin, fout, muptbins3, nmuptbins3);
    makeRatio1D("hMuEffPtIso", fin, fout, muptbins3, nmuptbins3);

    fout->Close();

    //Derive N(b) scale factors 
    TFile *fin2 = new TFile("histoutput.root");
    TH1 *h_1b =      (TH1*)fin2->Get("fake1b_baselineNoTag_nTopnTopCandSortedCntZinvnTopCandSortedCntZinvZ#rightarrow#nu#nu N(b) = 1single");
    TH1 *h_1b_fake = (TH1*)fin2->Get("fake1b_baselineNoTag_nTopnTopCandSortedCntZinv1bnTopCandSortedCntZinv1bZ#rightarrow#nu#nu N(b) = 0, 1 fake bsingle");
    TH1 *h_2b =      (TH1*)fin2->Get("fake2b_baselineNoTag_nTopnTopCandSortedCntZinvnTopCandSortedCntZinvZ#rightarrow#nu#nu N(b) = 2single");
    TH1 *h_2b_fake = (TH1*)fin2->Get("fake2b_baselineNoTag_nTopnTopCandSortedCntZinv2bnTopCandSortedCntZinv2bZ#rightarrow#nu#nu N(b) = 0, 2 fake bsingle");
    TH1 *h_3b =      (TH1*)fin2->Get("fake3b_baselineNoTag_nTopnTopCandSortedCntZinvnTopCandSortedCntZinvZ#rightarrow#nu#nu N(b) > 2single");
    TH1 *h_3b_fake = (TH1*)fin2->Get("fake3b_baselineNoTag_nTopnTopCandSortedCntZinv3bnTopCandSortedCntZinv3bZ#rightarrow#nu#nu N(b) = 0, 3 fake bsingle");
    
    printf("N(b) extrapolation scale factors\n");
    printf("N(b) = 0 -> 3: %e\n", h_1b->Integral()/h_1b_fake->Integral());
    printf("N(b) = 0 -> 3: %e\n", h_2b->Integral()/h_2b_fake->Integral());
    printf("N(b) = 0 -> 3: %e\n", h_3b->Integral()/h_3b_fake->Integral());
}
