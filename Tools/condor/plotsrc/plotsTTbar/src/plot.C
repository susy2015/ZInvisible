#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

#include <sstream>
#include <iomanip>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

//This is a helper function which will keep the plot from overlapping with the legend
void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = 0.0;
        if(error) bin = h->GetBinContent(i) + h->GetBinError(i);
        else      bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}

//Class to hold TH1* with various helper functions 
class histInfo
{
public:
    std::string legEntry, histFile, histName, drawOptions;
    int color, rebin;
    std::shared_ptr<TH1> h;
    std::shared_ptr<TH1> b; //We will save background histogram here, and automatically subtract it.
    double yield = 0;

    //Subtract a histogram (This is useful if we are using background subtraction)
    void subtract(TH1* bkg){

        b.reset(static_cast<TH1*>(bkg->Clone()));
        if(b){ std::cout << "Subtracting histogram has yield " << b->Integral() << std::endl; }
        //If h has already been retrieved, let's go ahead and do the subtraction
        if(b&&h){
            std::cout << "Performing background subtraction." << std::endl;
            TH1D *btemp = (TH1D*)b->Clone();
            h->Add(btemp,-1);
            std::cout << "New Yield: " << h->Integral() << std::endl;
        }

    }

    //Formay legEntry w/ yield
    std::string getlegEntry(){

    return legEntry+" ("+to_string_with_precision<double>(yield,3)+")";

    }

    //helper function to get histogram from file and configure its optional settings
    void retrieveHistogram()
    {
        //Open the file for this histogram
        TFile *f = TFile::Open(histFile.c_str());

        //check that the file was opened successfully
        if(!f)
        {
            printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
            h = nullptr;
            return;
        }

        //get the histogram from the file
        h.reset(static_cast<TH1*>(f->Get(histName.c_str())));

        //with the histogram retrieved, close the file
        f->Close();
        delete f;

        //check that the histogram was retireved from the file successfully
        if(!h)
        {
            printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
            return;
        }



        //Let's make the name nicer.
        std::string hName = (histName.find("/") != std::string::npos ? histName.substr(histName.find("/")+1) : histName);
        h->SetNameTitle(hName.c_str(),hName.c_str());


        //set the histogram color
        h->SetLineColor(color);
        h->SetLineWidth(3);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);

        //Get yield for legend
        yield = h->Integral();

        std::cout << yield << std::endl;

        // rebin the histogram if desired
        if(rebin >0) h->Rebin(rebin);

        //Now let's do a background subtraction if it has been specified
        if(h&&b){
            std::cout << "Performing background subtraction" << std::endl;
            TH1D *btemp = (TH1D*)b->Clone();
            h->Add(btemp,-1);
            std::cout << "New Yield is " << h->Integral() << std:: endl;
        }

    }

    //helper function for axes
    void setupAxes(float percent)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->SetTitleSize(0.045/percent);
        h->GetXaxis()->SetLabelSize(0.045/percent);
        h->GetYaxis()->SetTitleSize(0.045/percent);
        h->GetYaxis()->SetLabelSize(0.045/percent);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    void setupAxesR(float percent)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->SetTitleSize(0.045/percent);
        h->GetXaxis()->SetLabelSize(0.045/percent);
        h->GetYaxis()->SetTitleSize(0.045/percent);
        h->GetYaxis()->SetLabelSize(0.045/percent);
        h->GetYaxis()->SetRangeUser(0.,2.0);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
        h->GetYaxis()->SetNdivisions(5,4,0);
    }


    //wrapper to draw histogram
    void draw(const std::string& additionalOptions = "", bool noSame = false) const
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    void setFillColor(int newColor = -1)
    {
        if(newColor >= 0) h->SetFillColor(newColor);
        else              h->SetFillColor(color);
    }

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr)
    {
    }

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class Plotter
{
private:
    //entry for data
    histInfo data_;
    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries_;
    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries_;
    
public:
    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}

    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries) : data_(data), bgEntries_(bgEntries) {}

    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {

        bool bkgSub = sigEntries_.size() > 0;

        if(bkgSub){ 
            std::cout << "Background Subtraction is enabled." << std::endl;  
        }

        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);

        float percent = 0.2;

        //Divide the canvas in 2 so that we have room for a ratio plot as well (data/bkg)
        TPad *pad1 = new TPad("pad1","pad1", 0.0, percent, 1.0, 1.0, 0);
        TPad *pad2 = new TPad("pad2","pad2", 0.0, 0.0, 1.0, percent, 0);

        pad1->Draw();
        pad2->Draw();

        //switch to the canvas to ensure it is the active object
        pad1->cd();

        //Create TLegend
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        //Create the THStack for the background ... warning, THStacks are terrible and must be filled "backwards"
        THStack *bgStack = new THStack();
        //Make seperate histogram from sum of BG histograms because I don't know how to make a THStack give me this 
        TH1* hbgSum = nullptr;
        TH1* hsgSum = nullptr;

        for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
        {
            //Get new histogram
            bgEntries_[iBG].histName = histName;
            bgEntries_[iBG].rebin = rebin;
            bgEntries_[iBG].retrieveHistogram();

            //We are really making this code a frakenstein's monster, only fill this here if we aren't doing background subtraction
            if(!bkgSub) bgStack->Add(bgEntries_[iBG].h.get(), bgEntries_[iBG].drawOptions.c_str());

            if(!hbgSum) hbgSum = static_cast<TH1*>(bgEntries_[iBG].h->Clone());
            else        hbgSum->Add(bgEntries_[iBG].h.get());
        }

        //If we are doing background substration, setup the bgStack here:
        if(bkgSub){
            for(int iSG = 0; iSG < sigEntries_.size(); iSG++){
                sigEntries_[iSG].histName = histName;
                sigEntries_[iSG].rebin = rebin;
                sigEntries_[iSG].retrieveHistogram();
                bgStack->Add(sigEntries_[iSG].h.get(), sigEntries_[iSG].drawOptions.c_str());
                if(!hsgSum) hsgSum = static_cast<TH1*>(sigEntries_[iSG].h->Clone());
                else        hsgSum->Add(sigEntries_[iSG].h.get());
            } 
        }

        //data
        //get new histogram from file
        data_.histName = histName;
        data_.rebin = rebin;
        data_.retrieveHistogram();
        if(bkgSub) data_.subtract(hbgSum);
        leg->AddEntry(data_.h.get(), data_.getlegEntry().c_str(), data_.drawOptions.c_str());

        if(!bkgSub) {smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);}
        else {smartMax(hsgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);}

        //background
        if(!bkgSub){
            for(auto& entry : bgEntries_)
            {
                //set fill color so BG will have solid fill
                entry.setFillColor();

                //add histograms to TLegend
                leg->AddEntry(entry.h.get(), entry.getlegEntry().c_str(), "F");
            }
            smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }else{
            for(auto& entry : sigEntries_)
            {
                //get new histogram
                entry.setFillColor();;

                //add histograms to TLegend
                leg->AddEntry(entry.h.get(), entry.getlegEntry().c_str(), "F");
            }
            smartMax(hsgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }
        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        //gPad->SetBottomMargin(0.12);
        gPad->SetBottomMargin(0);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, data_.h->GetBinLowEdge(1), data_.h->GetBinLowEdge(data_.h->GetNbinsX()) + data_.h->GetBinWidth(data_.h->GetNbinsX())));
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //Here is the frakenstein monster, this is either all the MC, or it is just the wanted
        bgStack->Draw("same");

        //plot data histogram
        data_.draw();

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        //Now we work on the ratio plot
        TH1* ratio_hist;
        TH1* ratio_denom;        

        ratio_hist = static_cast<TH1*>(data_.h->Clone());
        if(!bkgSub){ratio_denom = (TH1D*)hbgSum->Clone();}
        else{ratio_denom = (TH1D*)hsgSum->Clone();}
        
        ratio_hist->Sumw2();
        ratio_denom->Sumw2();

        std::cout << "ratio_hist ";
        ratio_hist->Print();

        std::cout << "ratio_denom ";
        ratio_denom->Print();

        ratio_hist->Divide(ratio_denom); // data/bkg

        std::cout << "Ratio ";
        ratio_hist->Print();

        c->cd(2); // Let's move to the bottom pad
        pad2->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.25);

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //Draw a dummy
        //create a dummy histogram to act as the axes
        std::cout << "Setting up dummy ratio histogram" << std::endl;
        histInfo dummyR(new TH1D("dummyR", "dummyR", 1000, data_.h->GetBinLowEdge(1), data_.h->GetBinLowEdge(data_.h->GetNbinsX()) + data_.h->GetBinWidth(data_.h->GetNbinsX())));
        dummyR.setupAxesR(.4);
        if(bkgSub){dummyR.h->GetYaxis()->SetTitle("Ratio (data/bkg)");}
        else{dummyR.h->GetYaxis()->SetTitle("Ratio (data/MC)");}
        dummyR.h->GetXaxis()->SetTitle(xAxisLabel.c_str());

        //draw dummy axes
        std::cout << "Drawing ratio axes" << std::endl;
        dummyR.draw();

        //draw the ratio plot
        std::cout << "Now drawing ratio histogram" << std::endl;
        ratio_hist->Draw("same");
        std::cout << "Ratio histogram was drawn" << std::endl;

        std::string hName = (histName.find("/") != std::string::npos ? histName.substr(histName.find("/")+1) : histName);

        //save new plot to file
        std::cout << "Saving plot." << std::endl;
        if(!bkgSub){c->Print((hName + ".png").c_str());}
        else{c->Print(("bkgSub_" + hName + ".png").c_str());}

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;
        delete hsgSum;
        delete ratio_denom;
        delete ratio_hist;
    }


    void plotEff(bool data, const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {

        std::cout << "Called plotEff " << histNameNum << " " << histNameDen << std::endl;

        bool bkgSub = sigEntries_.size() > 0;

        if(bkgSub){ 
            std::cout << "Background Subtraction is enabled." << std::endl;  
        }

        if(data){
            std::cout << "Will plot from data histogram" << std::endl;
        } else {
            std::cout << "Will plot from background histogram" << std::endl;
        }

        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);

        float percent = 0.2;

        //Divide the canvas in 2 so that we have room for a ratio plot as well (data/bkg)
        TPad *pad1 = new TPad("pad1","pad1", 0.0, percent, 1.0, 1.0, 0);
        TPad *pad2 = new TPad("pad2","pad2", 0.0, 0.0, 1.0, percent, 0);

        pad1->Draw();
        pad2->Draw();

        //switch to the canvas to ensure it is the active object
        pad1->cd();

        //Create TLegend
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        TH1D* hNum;
        TH1D* hDen;

        //Get histograms from background if not plotting data
        if(!data){
            std::cout << "Loading background histograms" << std::endl;
            TH1D* hbgSumNum = nullptr;
            TH1D* hbgSumDen = nullptr;

            std::vector<histInfo> MCEntries_;

            if(bkgSub){ MCEntries_ = sigEntries_; }
            else{ MCEntries_ = bgEntries_; }

            for(int iBG = MCEntries_.size() - 1; iBG >= 0; --iBG)
            {
                //Get numerator histogram
                MCEntries_[iBG].histName = histNameNum;
                MCEntries_[iBG].rebin = rebin;
                MCEntries_[iBG].retrieveHistogram();

                if(!hbgSumNum){
                    hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                    hbgSumNum->SetDirectory(0);
                }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

                //Get denominator histogram
                MCEntries_[iBG].histName = histNameDen;
                MCEntries_[iBG].rebin = rebin;
                MCEntries_[iBG].retrieveHistogram();

                if(!hbgSumDen){
                    hbgSumDen = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                    hbgSumDen->SetDirectory(0);
                } else        hbgSumDen->Add(MCEntries_[iBG].h.get());
            }

            hNum = hbgSumNum;
            hDen = hbgSumDen;

            std::cout << "Finished background sum" << std::endl;

        } else { //We want to draw from the data histograms.

            std::cout << "loading data histograms" << std::endl;

            //We implemented subtraction in the HistInfo class, so we need to
            //Calculate the background to be subtracted before extracting
            //the histogram from histInfo.
            TH1D* hbgSumNum = nullptr; //Putting these here for scope reasons
            TH1D* hbgSumDen = nullptr;

            if(bkgSub){ //We need to subtract out the MC backgrounds
                std::cout << "Loading background histograms" << std::endl;
                for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
                {
                    //Get numerator histogram
                    bgEntries_[iBG].histName = histNameNum;
                    bgEntries_[iBG].rebin = rebin;
                    bgEntries_[iBG].retrieveHistogram();

                    if(!hbgSumNum){
                        hbgSumNum = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                        hbgSumNum->SetDirectory(0);
                    }else        hbgSumNum->Add(bgEntries_[iBG].h.get());

                    //Get denominator histogram
                    bgEntries_[iBG].histName = histNameDen;
                    bgEntries_[iBG].rebin = rebin;
                    bgEntries_[iBG].retrieveHistogram();

                    if(!hbgSumDen){
                        hbgSumDen = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                        hbgSumDen->SetDirectory(0);
                    } else        hbgSumDen->Add(bgEntries_[iBG].h.get());
                }
            }

            //data
            //get numerator histogram from file
            data_.histName = histNameNum;
            data_.rebin = rebin;
            if(bkgSub) data_.subtract(hbgSumNum);
            data_.retrieveHistogram();
            hNum = (static_cast<TH1D*>(data_.h.get()->Clone()));
            hNum->SetDirectory(0);
            std::cout << "hNum info" << std::endl;
            hNum->Print();

            //data
            //get denominator histogram from file
            data_.histName = histNameDen;
            data_.rebin = rebin;
            if(bkgSub) data_.subtract(hbgSumDen);
            data_.retrieveHistogram();
            hDen = static_cast<TH1D*>(data_.h.get()->Clone());
            hDen->SetDirectory(0);
            std::cout << "hDen info" << std::endl;
            hDen->Print();

            std::cout << "Loaded data histograms" << std::endl;

        }

        hNum->SetLineColor(kBlack);
        hNum->SetLineWidth(3);
        hNum->SetMarkerColor(kBlack);
        hNum->SetMarkerStyle(20);

        hDen->SetLineColor(kGray);
        hDen->SetLineWidth(3);
        hDen->SetMarkerColor(kGray);
        hDen->SetMarkerStyle(20);
        hDen->SetFillColor(kGray);


        TH1D* ratio_hist;
        
        std::cout << "Preparing ratio histogram" << std::endl;
        //hNum->ResetBit(kCanDelete);
        hNum->SetDirectory(0);
        hNum->Print();
        ratio_hist = (TH1D*)(hNum->Clone("ratio_hist"));

        std::cout << "Cloning of hNum complete" << std::endl;


        leg->AddEntry(hNum, ("Passes Top Tagger ("+to_string_with_precision<double>(hNum->Integral(),3)+")").c_str(), "PEX0");
        smartMax(hNum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        leg->AddEntry(hDen, ("All (" + to_string_with_precision<double>(hDen->Integral(),3)+")").c_str(), "F");
        smartMax(hDen, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        //background
//        for(auto& entry : bgEntries_)
//        {
            //set fill color so BG will have solid fill
//            entry.setFillColor();

            //add histograms to TLegend
            //leg->AddEntry(entry.h.get(), entry.getlegEntry().c_str(), "F");
//        }
        //smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        //gPad->SetBottomMargin(0.12);
        gPad->SetBottomMargin(0);

        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hDen->GetBinLowEdge(1), hDen->GetBinLowEdge(hDen->GetNbinsX()) + hDen->GetBinWidth(hDen->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //plot background stack
        hDen->Draw("same hist");

        //plot data histogram
        hNum->Draw("same PEX0");

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        //Now we work on the ratio plot

        ratio_hist->Sumw2();
        hDen->Sumw2();
        hNum->Sumw2();

        std::cout << "Printing numerator histogram" << std::endl;
        ratio_hist->Print();
        std::cout << "Printing denominator histogram" << std::endl;
        hDen->Print();

        ratio_hist->Divide(hDen); // data/bkg

        c->cd(2); // Let's move to the bottom pad
        pad2->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.25);

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //Draw a dummy
        //create a dummy histogram to act as the axes
        histInfo dummyR(new TH1D("dummyR", "dummyR", 1000, hDen->GetBinLowEdge(1), hDen->GetBinLowEdge(hDen->GetNbinsX()) + hDen->GetBinWidth(data_.h->GetNbinsX())));
        dummyR.setupAxesR(.4);
        dummyR.h->GetYaxis()->SetTitle("Tagging Rate");
        dummyR.h->GetXaxis()->SetTitle(xAxisLabel.c_str());

        //draw dummy axes
        dummyR.draw();

        //draw the ratio plot
        ratio_hist->Draw("same");

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        if(bkgSub){c->Print(("bkgSub_"+hNameDen + "TaggingRate"+(data ? "Data" : "MC")+".png").c_str());}
        else{c->Print((hNameDen + "TaggingRate"+(data ? "Data" : "MC")+".png").c_str());}

        //clean up dynamic memory
        delete c;
        delete leg;
    }

    void plotSF(const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {

        std::cout << "Called plotSF " << histNameNum << " " << histNameDen << std::endl;

        bool bkgSub = sigEntries_.size() > 0;

        if(bkgSub){ 
            std::cout << "Background Subtraction is enabled." << std::endl;  
        }

        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);

        float percent = 0.2;

        //Divide the canvas in 2 so that we have room for a ratio plot as well (data/bkg)
        TPad *pad1 = new TPad("pad1","pad1", 0.0, percent, 1.0, 1.0, 0);
        TPad *pad2 = new TPad("pad2","pad2", 0.0, 0.0, 1.0, percent, 0);

        pad1->Draw();
        pad2->Draw();

        //switch to the canvas to ensure it is the active object
        pad1->cd();

        //Create TLegend
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        TH1D* hNumData;
        TH1D* hDenData;

        TH1D* hNumMC;
        TH1D* hDenMC;

        TH1D* hNumSub = nullptr;
        TH1D* hDenSub = nullptr;


        //Get histograms from background if not plotting data

        std::cout << "Loading background histograms" << std::endl;
        TH1D* hbgSumNum = nullptr;
        TH1D* hbgSumDen = nullptr;
        std::vector<histInfo> MCEntries_;

        if(bkgSub){ MCEntries_ = sigEntries_; }
        else { MCEntries_ = bgEntries_; }

        for(int iBG = MCEntries_.size() - 1; iBG >= 0; --iBG)
        {
            //Get numerator histogram
            MCEntries_[iBG].histName = histNameNum;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumNum){
                hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumNum->SetDirectory(0);
            }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

            //Get denominator histogram
            MCEntries_[iBG].histName = histNameDen;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumDen){
                hbgSumDen = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumDen->SetDirectory(0);
            } else        hbgSumDen->Add(MCEntries_[iBG].h.get());
        }

        if(bkgSub){
            for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
            {
                //Get numerator histogram
                bgEntries_[iBG].histName = histNameNum;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].retrieveHistogram();

                if(!hNumSub){
                    hNumSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hNumSub->SetDirectory(0);
                }else        hNumSub->Add(bgEntries_[iBG].h.get());
 
                //Get denominator histogram
                bgEntries_[iBG].histName = histNameDen;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].retrieveHistogram();

                if(!hDenSub){
                    hDenSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hDenSub->SetDirectory(0);
                } else        hDenSub->Add(bgEntries_[iBG].h.get());
            }
        }

        hNumMC = hbgSumNum;
        hDenMC = hbgSumDen;

        hNumMC->Sumw2();
        hDenMC->Sumw2();

        std::cout << "Finished background sum" << std::endl;


        std::cout << "loading data histograms" << std::endl;

        //data
        //get numerator histogram from file
        data_.histName = histNameNum;
        data_.rebin = rebin;
        if(bkgSub)data_.subtract(hNumSub);
        data_.retrieveHistogram();
        hNumData = (static_cast<TH1D*>(data_.h.get()->Clone()));
        hNumData->Sumw2();
        hNumData->SetDirectory(0);
        std::cout << "hNum info" << std::endl;
        hNumData->Print();

        //data
        //get denominator histogram from file
        data_.histName = histNameDen;
        data_.rebin = rebin;
        if(bkgSub) data_.subtract(hDenSub);
        data_.retrieveHistogram();
        hDenData = static_cast<TH1D*>(data_.h.get()->Clone());
        hDenData->Sumw2();
        hDenData->SetDirectory(0);
        std::cout << "hDen info" << std::endl;
        hDenData->Print();

        std::cout << "Loaded data histograms" << std::endl;


        TH1D *hMC = (TH1D*)hNumMC->Clone();
        hMC->Sumw2();
        hMC->Divide(hDenMC);
        hMC->Sumw2();

        TH1D *hData = (TH1D*)hNumData->Clone();
        hData->Sumw2();
        hData->Divide(hDenData);
        hData->Sumw2();

        hMC->SetLineColor(kRed);
        hMC->SetLineWidth(3);
        hMC->SetMarkerColor(kRed);
        hMC->SetMarkerStyle(20);

        hData->SetLineColor(kBlue);
        hData->SetLineWidth(3);
        hData->SetMarkerColor(kBlue);
        hData->SetMarkerStyle(20);
        hData->SetFillColor(kGray);


        TH1D* ratio_hist;
        
        std::cout << "Preparing ratio histogram" << std::endl;
        //hMC->ResetBit(kCanDelete);
        hMC->SetDirectory(0);
        hMC->Print();
        ratio_hist = (TH1D*)(hData->Clone("ratio_hist"));
        ratio_hist->Sumw2();

        std::cout << "Cloning of hMC complete" << std::endl;


        leg->AddEntry(hMC, "MC Tagging Rate", "PEX0");
        smartMax(hMC, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        leg->AddEntry(hData, "Data Tagging Rate", "PEX0");
        smartMax(hData, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        //background
//        for(auto& entry : bgEntries_)
//        {
            //set fill color so BG will have solid fill
//            entry.setFillColor();

            //add histograms to TLegend
            //leg->AddEntry(entry.h.get(), entry.getlegEntry().c_str(), "F");
//        }
        //smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        //gPad->SetBottomMargin(0.12);
        gPad->SetBottomMargin(0);

        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hData->GetBinLowEdge(1), hData->GetBinLowEdge(hData->GetNbinsX()) + hData->GetBinWidth(hData->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //plot background stack
        hData->Draw("same PEX0");

        //plot data histogram
        hMC->Draw("same PEX0");

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        //Now we work on the ratio plot

        ratio_hist->Sumw2();
        hData->Sumw2();

        std::cout << "Printing numerator histogram" << std::endl;
        ratio_hist->Print();
        std::cout << "Printing denominator histogram" << std::endl;
        hData->Print();

        ratio_hist->Divide(hMC); // data/MC
        ratio_hist->Sumw2();

        ratio_hist->SetLineColor(kBlack);
        ratio_hist->SetLineWidth(3);
        ratio_hist->SetMarkerColor(kBlack);
        ratio_hist->SetMarkerStyle(20);

        c->cd(2); // Let's move to the bottom pad
        pad2->cd();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.25);

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //Draw a dummy
        //create a dummy histogram to act as the axes
        histInfo dummyR(new TH1D("dummyR", "dummyR", 1000, hData->GetBinLowEdge(1), hData->GetBinLowEdge(hData->GetNbinsX()) + hData->GetBinWidth(data_.h->GetNbinsX())));
        dummyR.setupAxesR(.4);
        dummyR.h->GetYaxis()->SetTitle("Tagging Rate");
        dummyR.h->GetXaxis()->SetTitle(xAxisLabel.c_str());

        //draw dummy axes
        dummyR.draw();

        //draw the ratio plot
        ratio_hist->Draw("same");

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        if(bkgSub){c->Print(("bkgSub_" + hNameDen + "SF.png").c_str());}
        else{c->Print((hNameDen + "SF.png").c_str());}

        //clean up dynamic memory
        delete c;
        delete leg;
    }

    void plotOnlySF(const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {

        std::cout << "Called plotOnlySF " << histNameNum << " " << histNameDen << std::endl;

        bool bkgSub = sigEntries_.size() > 0;

        if(bkgSub){ 
            std::cout << "Background Subtraction is enabled." << std::endl;  
        }


        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);

        float percent = 0.2;

        //Create TLegend
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        TH1D* hNumData;
        TH1D* hDenData;

        TH1D* hNumMC;
        TH1D* hDenMC;

        TH1D* hNumSub = nullptr;
        TH1D* hDenSub = nullptr;


        std::cout << "Loading background histograms" << std::endl;
        TH1D* hbgSumNum = nullptr;
        TH1D* hbgSumDen = nullptr;

        std::vector<histInfo> MCEntries_;

        if(bkgSub){ MCEntries_ = sigEntries_; }
        else{ MCEntries_ = bgEntries_; }


        for(int iBG = MCEntries_.size() - 1; iBG >= 0; --iBG)
        {
            //Get numerator histogram
            MCEntries_[iBG].histName = histNameNum;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumNum){
                hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumNum->SetDirectory(0);
            }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

            //Get denominator histogram
            MCEntries_[iBG].histName = histNameDen;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumDen){
                hbgSumDen = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumDen->SetDirectory(0);
            } else        hbgSumDen->Add(MCEntries_[iBG].h.get());
        }

        if(bkgSub){
            for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
            {
                //Get numerator histogram
                bgEntries_[iBG].histName = histNameNum;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].retrieveHistogram();

                if(!hNumSub){
                    hNumSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hNumSub->SetDirectory(0);
                }else        hNumSub->Add(bgEntries_[iBG].h.get());

                //Get denominator histogram
                bgEntries_[iBG].histName = histNameDen;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].retrieveHistogram();

                if(!hDenSub){
                    hDenSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hDenSub->SetDirectory(0);
                } else        hDenSub->Add(bgEntries_[iBG].h.get());
            }
        }


        hNumMC = hbgSumNum;
        hDenMC = hbgSumDen;

        hNumMC->Sumw2();
        hDenMC->Sumw2();

        std::cout << "Finished background sum" << std::endl;


        std::cout << "loading data histograms" << std::endl;

        //data
        //get numerator histogram from file
        data_.histName = histNameNum;
        data_.rebin = rebin;
        if(bkgSub) data_.subtract(hNumSub);
        data_.retrieveHistogram();
        hNumData = (static_cast<TH1D*>(data_.h.get()->Clone()));
        hNumData->Sumw2();
        hNumData->SetDirectory(0);
        std::cout << "hNum info" << std::endl;
        hNumData->Print();

        //data
        //get denominator histogram from file
        data_.histName = histNameDen;
        data_.rebin = rebin;
        if(bkgSub) data_.subtract(hDenSub);
        data_.retrieveHistogram();
        hDenData = static_cast<TH1D*>(data_.h.get()->Clone());
        hDenData->Sumw2();
        hDenData->SetDirectory(0);
        std::cout << "hDen info" << std::endl;
        hDenData->Print();

        std::cout << "Loaded data histograms" << std::endl;


        TH1D *hMC = (TH1D*)hNumMC->Clone();
        hMC->Sumw2();
        hMC->Divide(hDenMC);
        hMC->Sumw2();

        TH1D *hData = (TH1D*)hNumData->Clone();
        hData->Sumw2();
        hData->Divide(hDenData);
        hData->Sumw2();

        hMC->SetLineColor(kRed);
        hMC->SetLineWidth(3);
        hMC->SetMarkerColor(kRed);
        hMC->SetMarkerStyle(20);

        hData->SetLineColor(kBlue);
        hData->SetLineWidth(3);
        hData->SetMarkerColor(kBlue);
        hData->SetMarkerStyle(20);
        hData->SetFillColor(kGray);


        TH1D* ratio_hist;
        
        std::cout << "Preparing ratio histogram" << std::endl;
        //hMC->ResetBit(kCanDelete);
        hMC->SetDirectory(0);
        hMC->Print();
        ratio_hist = (TH1D*)(hData->Clone("ratio_hist"));
        ratio_hist->Sumw2();
        hData->Sumw2();

        std::cout << "Printing numerator histogram" << std::endl;
        ratio_hist->Print();
        std::cout << "Printing denominator histogram" << std::endl;
        hData->Print();

        ratio_hist->Divide(hMC); // data/MC
        ratio_hist->Sumw2();

        ratio_hist->SetLineColor(kBlack);
        ratio_hist->SetLineWidth(3);
        ratio_hist->SetMarkerColor(kBlack);
        ratio_hist->SetMarkerStyle(20);



        leg->AddEntry(hMC, "MC Tagging Rate", "PEX0");
        smartMax(hMC, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        leg->AddEntry(hData, "Data Tagging Rate", "PEX0");
        smartMax(hData, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        //background
//        for(auto& entry : bgEntries_)
//        {
            //set fill color so BG will have solid fill
//            entry.setFillColor();

            //add histograms to TLegend
            //leg->AddEntry(entry.h.get(), entry.getlegEntry().c_str(), "F");
//        }
        //smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);
        //gPad->SetBottomMargin(0);



        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hData->GetBinLowEdge(1), hData->GetBinLowEdge(hData->GetNbinsX()) + hData->GetBinWidth(hData->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        //Set the y-range of the histogram
/*        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        } */
        double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 2;
        ymax = (ymax < 1.5 ? 1.5 : ymax);
        dummy.h->GetYaxis()->SetRangeUser(0.0, ymax);


        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        //gPad->SetLeftMargin(0.12);
        //gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        //gPad->SetTopMargin(0);
        //gPad->SetBottomMargin(0.25);

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //draw the ratio plot
        ratio_hist->Draw("same");

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        if(bkgSub){c->Print(("bkgSub_" + hNameDen + "OnlySF.png").c_str());}
        else{c->Print((hNameDen + "OnlySF.png").c_str());}

        //plot legend
        //leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);


        //clean up dynamic memory
        delete c;
        delete leg;
    }
};

int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_MET-combined.root", "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
//        {"GJets",         "GJets-combined.root",      "hist", kViolet},   
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
//        {"DYJets",        "DYJets-combined.root",     "hist", kCyan},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring}, 
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //These are the background that we want to subtract from the data
    std::vector<histInfo> unwantedEntries = {
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring}, 
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3}
    };

    //These are the entires that have the objects we want (hadronically decaying tops)
    std::vector<histInfo> wantedEntries = {
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //vector summarizing signal histograms to include in the plot
//    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
//    };

    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));
    Plotter pltSub(std::move(data), std::move(unwantedEntries), std::move(wantedEntries)); //Plotter to make background subtracted plots

    plt.plot("ttbar/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbar/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/nJets", "Jets/Event", "Events", true);
    plt.plot("ttbar/nVertices", "primary vertices/event", "Events", true);
    plt.plot("ttbar/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("ttbar/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);

    pltSub.plot("ttbar/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/nJets", "Jets/Event", "Events", true);
    pltSub.plot("ttbar/nVertices", "primary vertices/event", "Events", true);
    pltSub.plot("ttbar/photon", "Photon p_{T} [GeV]", "Events", true);
    pltSub.plot("ttbar/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("ttbar/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbar/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    pltSub.plot("ttbar/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    plt.plotEff(true, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    pltSub.plotEff(true, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotSF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    pltSub.plotSF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    pltSub.plotOnlySF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(true, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(false, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(false, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    pltSub.plotSF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    pltSub.plotSF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

    pltSub.plotOnlySF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    pltSub.plotOnlySF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
