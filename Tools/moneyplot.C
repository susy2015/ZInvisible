#include "tdrstyle.h"
#include "ScaleFactors.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

void smartMax(const TH1* const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}


int main(int argc, char* argv[])
{

    // Get the relevant information
    TFile* f1 = TFile::Open("/uscms_data/d3/nstrobbe/HadronicStop/DataTest/CMSSW_7_4_8/src/ZInvisible/Tools/condor/dataplots_muon_Dec15.root");
    TH1D* h1 = (TH1D*)f1->Get("nSearchBin/NJetWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nusingle");
    // Scale the prediction by the normalization factor by hand for now
    h1->Scale(ScaleFactors::sf_norm0b());

    //TFile* f2 = TFile::Open("/uscms/home/pastika/nobackup/zinv/dev/CMSSW_7_4_8/src/ZInvisible/Tools/syst_nJetWgt.root");
    TFile* f2 = TFile::Open("syst_all.root");
    TH1D* h2 = (TH1D*)f2->Get("shape_central");
    TH1D* h3 = (TH1D*)f2->Get("shape_stat");
    TH1D* h4 = (TH1D*)f2->Get("MC_stats");

    //TGraphAsymmErrors* g1 = (TGraphAsymmErrors*)f2->Get("");
    TGraphAsymmErrors* g3 = new TGraphAsymmErrors();
    TGraphAsymmErrors* g4 = new TGraphAsymmErrors();
    const int n = h2->GetNbinsX();
    double x[n];
    double y[n];
    double exl[n];
    double exh[n];
    double eyl_1[n];
    double eyh_1[n];
    double eyl_2[n];
    double eyh_2[n];
    double rel_unc_2 = ScaleFactors::sfunc_norm0b()/ScaleFactors::sf_norm0b();
    double e2_temp;
    for(int i = 1; i < n+1; ++i)
    {
	x[i-1] = h1->GetBinCenter(i);
	y[i-1] = h1->GetBinContent(i);
	exl[i-1] = 0.5;
	exh[i-1] = 0.5;
	eyl_1[i-1] = y[i-1]*rel_unc_2;
	eyh_1[i-1] = y[i-1]*rel_unc_2;
	e2_temp = y[i-1]*h2->GetBinContent(i);
	eyl_2[i-1] = sqrt(eyl_1[i-1]*eyl_1[i-1] + e2_temp*e2_temp);
	eyh_2[i-1] = sqrt(eyh_1[i-1]*eyh_1[i-1] + e2_temp*e2_temp);
      
        double err = h3->GetBinContent(i) * h1->GetBinContent(i);
        err = sqrt(err*err + eyl_2[i-1]*eyl_2[i-1]);
        g3->SetPoint(i - 1, h1->GetBinCenter(i), h1->GetBinContent(i));
        g3->SetPointError(i - 1, 0.5, 0.5, err, err);

        err = sqrt(err*err + pow(h4->GetBinContent(i) * h1->GetBinContent(i), 2));
        g4->SetPoint(i - 1, h1->GetBinCenter(i), h1->GetBinContent(i));
        g4->SetPointError(i - 1, 0.5, 0.5, err, err);
	std::cout << "bin " << i << ", rel unc (njet): " << h2->GetBinContent(i) << ", rel unc (norm): " << rel_unc_2  << "" << std::endl;
    }

    TGraphAsymmErrors* g1 = new TGraphAsymmErrors(n,x,y,exl,exh,eyl_1,eyh_1);
    TGraphAsymmErrors* g2 = new TGraphAsymmErrors(n,x,y,exl,exh,eyl_2,eyh_2);


    // Set the style
    setTDRStyle();

    // Prepare canvas
    TCanvas *c;
    double fontScale;
    bool showRatio = false;
    if(showRatio)
    {
	c = new TCanvas("c1", "c1", 800, 900);
	c->Divide(1, 2);
	c->cd(1);
	gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
	gPad->SetBottomMargin(0.01);
	fontScale = 1.0;
    }
    else
    {
	c = new TCanvas("c1", "c1", 1200, 800);
	c->Divide(1, 1);
	c->cd(1);
	gPad->SetPad("p1", "p1", 0, 0, 1, 1, kWhite, 0, 0);
	gPad->SetBottomMargin(0.12);
	fontScale = 6.5 / 8;
    }
    
    
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.06 * (8.0 / 6.5) * fontScale);
    
    
    // Use dummy histogram to set the axes and stuff 
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()) + h1->GetBinWidth(h1->GetNbinsX()));
    dummy->SetStats(0);
    dummy->SetTitle(0);
    dummy->GetYaxis()->SetTitle("Events");
    dummy->GetYaxis()->SetTitleOffset(.9*1.05 / (fontScale));
    dummy->GetXaxis()->SetTitleOffset(1.05);
    if(showRatio) dummy->GetXaxis()->SetTitle("");
    else          dummy->GetXaxis()->SetTitle("Search Bin");
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    TLegend *leg = new TLegend(0.55, 0.73, 0.89, 0.89);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    double max = 0.0, lmax = 0.0, min = 1.0e300, minAvgWgt = 1.0e300;
    int iSingle = 0, iRatio = 0;
    char legEntry[128];

    // Format the histogram
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(3);
    iSingle++;
    double integral = h1->Integral(0, h1->GetNbinsX() + 1);
    sprintf(legEntry, "%s", "prediction");
    leg->AddEntry(h1, legEntry);
    smartMax(h1, leg, static_cast<TPad*>(gPad), min, max, lmax);
    minAvgWgt = std::min(minAvgWgt, h1->GetSumOfWeights()/h1->GetEntries());

    // Format errors
    g1->SetFillColor(kRed-7);
    sprintf(legEntry, "%s", "Normalization unc.");
    leg->AddEntry(g1, legEntry);
    g2->SetFillColor(kCyan-6);
    //sprintf(legEntry, "%s", "Njet reweighting unc.");
    sprintf(legEntry, "%s", "Data/MC shape uncertinty");
    leg->AddEntry(g2, legEntry);
    g3->SetFillColor(kMagenta);
    sprintf(legEntry, "%s", "Njet/shape stat. unc.");
    leg->AddEntry(g3, legEntry);
    g4->SetFillColor(kGreen+2);
    sprintf(legEntry, "%s", "MC Stats");
    leg->AddEntry(g4, legEntry);

    bool isLog = false;
    gPad->SetLogy(isLog);
    if(isLog)
    {
	double locMin = std::min(0.2*minAvgWgt, std::max(0.0001, 0.05 * min));
	double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
	double legMin = legSpan + log10(locMin);
	if(log10(lmax) > legMin)
	{
	    double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
	    max = pow(max/locMin, scale)*locMin;
	}
	//dummy->GetYaxis()->SetRangeUser(locMin, 50*max);
        dummy->GetYaxis()->SetRangeUser(0.001, 15*max);
    }
    else
    {
	double locMin = 0.0;
	double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
	if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
	dummy->GetYaxis()->SetRangeUser(0.0, max*1.4);
    }

    dummy->Draw();
    
    g4->Draw("2 same");
    g3->Draw("2 same");
    g2->Draw("2 same");
    g1->Draw("2 same");
    h1->Draw("histsame");
    leg->Draw();

    fixOverlay();

    c->cd(1);
    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 13 TeV", 2156. / 1000.0);
    TLatex mark;
    mark.SetTextSize(0.042 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    mark.DrawLatex(gPad->GetLeftMargin(), 0.95, "CMS Preliminary");
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 0.95, lumistamp);
    
    fixOverlay();
    c->Print("moneyplot.png");
    c->Print("moneyplot.pdf");

    f1->Close();
    
    return 0;
}
