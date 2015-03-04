#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <cstdio>

TH2* rebin2d(TH2* old)
{
    //create a new TH2 with your bin arrays spec
    double xbins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 340.0, 380.0, 420.0, 460.0, 500.0, 560.0, 2000.0};
    int nx = sizeof(xbins)/sizeof(double) - 1;

    double ybins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 340.0, 380.0, 420.0, 460.0, 500.0, 580.0, 660.0, 800.0, 940.0, 3000.0};
    int ny = sizeof(ybins)/sizeof(double) - 1;

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

int main ()
{
    TFile  *fin = new TFile("effhists.root");

    TH1 *hMuEffPt_num = (TH1*)fin->Get("hMuEffPt_num");
    TH1 *hMuEffPt_den = (TH1*)fin->Get("hMuEffPt_den");
    TH1 *hMuAccPt_num = (TH1*)fin->Get("hMuAccPt_num");
    TH1 *hMuAccPt_den = (TH1*)fin->Get("hMuAccPt_den");

    TH1 *hMuEffHt_num = (TH1*)fin->Get("hMuEffHt_num");
    TH1 *hMuEffHt_den = (TH1*)fin->Get("hMuEffHt_den");
    TH1 *hMuAccHt_num = (TH1*)fin->Get("hMuAccHt_num");
    TH1 *hMuAccHt_den = (TH1*)fin->Get("hMuAccHt_den");

    TH2 *hMuEff_num = (TH2*)fin->Get("hMuEff_num");
    TH2 *hMuEff_den = (TH2*)fin->Get("hMuEff_den");
    TH2 *hMuAcc_num = (TH2*)fin->Get("hMuAcc_num");
    TH2 *hMuAcc_den = (TH2*)fin->Get("hMuAcc_den");

    TH2 *hMuEff_num_pp = (TH2D*)fin->Get("hMuEff_num_pp");
    TH2 *hMuEff_num_fp = (TH2D*)fin->Get("hMuEff_num_fp");
    TH2 *hMuEff_num_pf = (TH2D*)fin->Get("hMuEff_num_pf");

    TH2 *hMuEff_num_rand = (TH2D*)fin->Get("hMuEff_num_rand");

    //Rebin
    hMuEff_num = rebin2d(hMuEff_num);
    hMuEff_den = rebin2d(hMuEff_den);
    hMuAcc_num = rebin2d(hMuAcc_num);
    hMuAcc_den = rebin2d(hMuAcc_den);

    hMuEff_num_pp = rebin2d(hMuEff_num_pp);
    hMuEff_num_fp = rebin2d(hMuEff_num_fp);
    hMuEff_num_pf = rebin2d(hMuEff_num_pf);

    hMuEff_num_rand = rebin2d(hMuEff_num_rand);


    //Clone numerators
    TH1 *hMuEffPt = (TH1*)hMuEffPt_num->Clone("hMuEffPt");
    TH1 *hMuAccPt = (TH1*)hMuAccPt_num->Clone("hMuAccPt");

    TH1 *hMuEffHt = (TH1*)hMuEffHt_num->Clone("hMuEffHt");
    TH1 *hMuAccHt = (TH1*)hMuAccHt_num->Clone("hMuAccHt");

    TH2 *hMuEff = (TH2*)hMuEff_num->Clone("hMuEff");
    TH2 *hMuAcc = (TH2*)hMuAcc_num->Clone("hMuAcc");

    TH2 *hMuEff_pp = (TH2D*)hMuEff_num_pp->Clone("hMuEff_pp");
    TH2 *hMuEff_pf = (TH2D*)hMuEff_num_pf->Clone("hMuEff_pf");
    TH2 *hMuEff_fp = (TH2D*)hMuEff_num_fp->Clone("hMuEff_fp");

    TH2 *hMuEff_rand = (TH2D*)hMuEff_num_rand->Clone("hMuEff_rand");    

    hMuEffPt->Divide(hMuEffPt_den);
    hMuAccPt->Divide(hMuAccPt_den);

    hMuEffHt->Divide(hMuEffHt_den);
    hMuAccHt->Divide(hMuAccHt_den);

    hMuEff->Divide(hMuEff_den);
    hMuAcc->Divide(hMuAcc_den);

    hMuEff_pp->Divide(hMuEff_den);
    hMuEff_pf->Divide(hMuEff_den);
    hMuEff_fp->Divide(hMuEff_den);

    hMuEff_rand->Divide(hMuEff_den);

    TFile *fout = new TFile("muEffHists2.root", "RECREATE");
    fout->cd();

    hMuEffPt->Write();
    hMuAccPt->Write();

    hMuEffHt->Write();
    hMuAccHt->Write();

    hMuEff->Write();
    hMuAcc->Write();

    hMuEff_pp->Write();
    hMuEff_pf->Write();
    hMuEff_fp->Write();

    hMuEff_rand->Write();


    fin->Close();
}
