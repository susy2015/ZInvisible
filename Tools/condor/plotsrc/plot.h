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
    double merge; //This provides a value, any bin that includes this value or is greater, will be merged into a single bin
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
        std::cout << "Attempting to retrieve " << histName << std::endl;

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

        std::cout << hName << " has a yield of " << yield << "." << std::endl;
        std::cout << hName << " has " << h->GetNbinsX() << " bins." << std::endl;

        // rebin the histogram if desired
        if(rebin >0){
            std::cout << "Attempting to rebin " << hName << "." << std::endl;
            h->Rebin(rebin);
            std::cout << "After rebinning, " << hName << " has " << h->GetNbinsX() << "bins." << std::endl;
        }

        //merge last bin if necessary
        if(merge > -999){ //i'm not using -1 for the null case :)
            std::cout << "Merging bins for values greater than " << merge << std::endl;
            int nbins = h->GetNbinsX();
            int firstbin = h->GetXaxis()->FindBin(merge);
            int newbins = firstbin; //the bin number of this last bin is the number of bins the reshaped histogram will have
            double xbins[newbins+1]; //With 'newbins' number of bins, we need to specify 'newbins+1' number of edges.
            std::cout << h->GetName() << " has " << nbins << " bins. After merging it will have " << newbins << " bins." << std::endl;

            if(firstbin > -1){
                for(int i = 1; i <= newbins; i++){ //the first bin is '1'
                    xbins[i-1] = h->GetXaxis()->GetBinLowEdge(i);
                }
                xbins[newbins] = h->GetXaxis()->GetBinUpEdge(nbins); // Get the upper edge of the last bin in the original
                std::string hName = std::string("rebinned_")+h->GetName();
                h = (std::shared_ptr<TH1>)h->Rebin(newbins, hName.c_str(), xbins);
            }else{
               std::cout << "Merging bin was requested, but no action taken since there are no bins above " << merge << std::endl;
            }

            std::cout << "After merging " << h->GetName() << " has " << h->GetNbinsX() << " bins." << std::endl;
        }

        //Now let's do a background subtraction if it has been specified
        if(h&&b){
            std::cout << "Performing background subtraction on " << hName << std::endl;
            TH1D *btemp = (TH1D*)b->Clone();
            h->Add(btemp,-1);
            std::cout << "New Yield for " << hName << " is " << h->Integral() << std:: endl;
        }

        std::cout << "Attempted Retrieval of " << histName << " complete." << std::endl << std::endl;

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

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr), merge(-1000)
    {
    }

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), h(h), merge(-1000)
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

    //Identifying information, especially useful when we do multiplots
    std::string prefix = "";
    std::string label  = "";
    
public:
    void editPrefix(std::string pre){prefix = pre;}
    std::string getPrefix(){return prefix;}

    void editLabel(std::string labl){label = labl;}
    std::string getLabel(){return label;}

    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}

    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries) : data_(data), bgEntries_(bgEntries) {}

    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
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
            bgEntries_[iBG].merge = merge;
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
                sigEntries_[iSG].merge = merge;
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
        data_.merge = merge;
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

        //std::cout << "ratio_hist ";
        //ratio_hist->Print();

        //std::cout << "ratio_denom ";
        //ratio_denom->Print();

        ratio_hist->Divide(ratio_denom); // data/bkg

        //std::cout << "Ratio ";
        //ratio_hist->Print();

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
        std::string fName = hName + ".png";
        if(bkgSub) fName = "bkgSub_"+fName;
        
        c->Print(fName.c_str());

        std::string fNamePDF = hName + ".pdf";
        if(bkgSub) fNamePDF = "bkgSub_"+fNamePDF;
        
        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;
        delete hsgSum;
        delete ratio_denom;
        delete ratio_hist;
    }


    void plotEff(bool data, const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
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
                MCEntries_[iBG].merge = merge;
                MCEntries_[iBG].retrieveHistogram();

                if(!hbgSumNum){
                    hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                    hbgSumNum->SetDirectory(0);
                }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

                //Get denominator histogram
                MCEntries_[iBG].histName = histNameDen;
                MCEntries_[iBG].rebin = rebin;
                MCEntries_[iBG].merge = merge;
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
                    bgEntries_[iBG].merge = merge;
                    bgEntries_[iBG].retrieveHistogram();

                    if(!hbgSumNum){
                        hbgSumNum = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                        hbgSumNum->SetDirectory(0);
                    }else        hbgSumNum->Add(bgEntries_[iBG].h.get());

                    //Get denominator histogram
                    bgEntries_[iBG].histName = histNameDen;
                    bgEntries_[iBG].rebin = rebin;
                    bgEntries_[iBG].merge = merge;
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
            data_.merge = merge;
            if(bkgSub) data_.subtract(hbgSumNum);
            data_.retrieveHistogram();
            hNum = (static_cast<TH1D*>(data_.h.get()->Clone()));
            hNum->SetDirectory(0);
            //std::cout << "hNum info" << std::endl;
            //hNum->Print();

            //data
            //get denominator histogram from file
            data_.histName = histNameDen;
            data_.rebin = rebin;
            data_.merge = merge;
            if(bkgSub) data_.subtract(hbgSumDen);
            data_.retrieveHistogram();
            hDen = static_cast<TH1D*>(data_.h.get()->Clone());
            hDen->SetDirectory(0);
            //std::cout << "hDen info" << std::endl;
            //hDen->Print();

            //std::cout << "Loaded data histograms" << std::endl;

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
        
        //std::cout << "Preparing ratio histogram" << std::endl;
        //hNum->ResetBit(kCanDelete);
        hNum->SetDirectory(0);
        //hNum->Print();
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
        //std::cout << "ratio_hist has " << ratio_hist->GetNbinsX() << " bins." << std::endl;
        hDen->Sumw2();
        //std::cout << "hDen has " << hDen->GetNbinsX() << " bins." << std::endl;
        hNum->Sumw2();
        //std::cout << "hNum has " << hNum->GetNbinsX() << " bins." << std::endl;

        //std::cout << "Printing numerator histogram" << std::endl;
        ratio_hist->Print();
        //std::cout << "Printing denominator histogram" << std::endl;
        //hDen->Print();

        ratio_hist->Divide(hDen); // data/bkg

        //std::cout << "Successfully completed ratio division." << std::endl;

        c->cd(2); // Let's move to the bottom pad
        pad2->cd();

        //std::cout << "Working with bottom pad." << std::endl;

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.25);

        std::cout << "gPad set margins." << std::endl;

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        std::cout << "gPad set tick and grid" << std::endl;

        //Draw a dummy
        //create a dummy histogram to act as the axes
        histInfo dummyR(new TH1D("dummyR", "dummyR", 1000, hDen->GetBinLowEdge(1), hDen->GetBinLowEdge(hDen->GetNbinsX()) + hDen->GetBinWidth(hDen->GetNbinsX())));
        std::cout << "Created histInfo object for dummyR" << std::endl;
        dummyR.setupAxesR(.4);
        dummyR.h->GetYaxis()->SetTitle("Tagging Rate");
        dummyR.h->GetXaxis()->SetTitle(xAxisLabel.c_str());

        double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 1.5;
        ymax = (ymax < 1.5 ? ymax : 2);
        dummyR.h->GetYaxis()->SetRangeUser(0.0, ymax);

        std::cout << "Created dummy axes histogram for ratio plot." << std::endl;

        //draw dummy axes
        dummyR.draw();

        std::cout << "Drew dummyR histogram." << std::endl;

        //draw the ratio plot
        ratio_hist->Draw("same");

        std::cout << "Drew ratio histogram." << std::endl;

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        //std::cout << "Trying to save tagging rate plot." << std::endl;
        std::string fName = "TaggingRate"+hNameDen;
        if(data){fName = fName+"Data.png";}else{fName = fName+"MC.png";}
        if(bkgSub) fName = "bkgSub_" + fName;
        c->Print(fName.c_str());

        std::string fNamePDF = "TaggingRate+hNameDen";
        if(data){fNamePDF = fNamePDF+"Data.pdf";}else{fNamePDF = fNamePDF+"MC.pdf";}
        if(bkgSub) fNamePDF = "bkgSub_" + fNamePDF;
        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
    }

    void plotSF(const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double combineLast = -1, double lumi = 36100)
    {

        std::cout << "Called plotSF " << histNameNum << " " << histNameDen << std::endl;

        bool bkgSub = sigEntries_.size() > 0;

        bool combine = combineLast > 0;

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
            MCEntries_[iBG].merge = merge;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumNum){
                hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumNum->SetDirectory(0);
            }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

            //Get denominator histogram
            MCEntries_[iBG].histName = histNameDen;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].merge = merge;
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
                bgEntries_[iBG].merge = merge;
                bgEntries_[iBG].retrieveHistogram();

                if(!hNumSub){
                    hNumSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hNumSub->SetDirectory(0);
                }else        hNumSub->Add(bgEntries_[iBG].h.get());
 
                //Get denominator histogram
                bgEntries_[iBG].histName = histNameDen;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].merge = merge;
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
        data_.merge = merge;
        if(bkgSub)data_.subtract(hNumSub);
        data_.retrieveHistogram();
        hNumData = (static_cast<TH1D*>(data_.h.get()->Clone()));
        hNumData->Sumw2();
        hNumData->SetDirectory(0);
        //std::cout << "hNum info" << std::endl;
        //hNumData->Print();

        //data
        //get denominator histogram from file
        data_.histName = histNameDen;
        data_.rebin = rebin;
        data_.merge = merge;
        if(bkgSub) data_.subtract(hDenSub);
        data_.retrieveHistogram();
        hDenData = static_cast<TH1D*>(data_.h.get()->Clone());
        hDenData->Sumw2();
        hDenData->SetDirectory(0);
        //std::cout << "hDen info" << std::endl;
        //hDenData->Print();

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
        
        //std::cout << "Preparing ratio histogram" << std::endl;
        //hMC->ResetBit(kCanDelete);
        hMC->SetDirectory(0);
        //hMC->Print();
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

        //std::cout << "Printing numerator histogram" << std::endl;
        //ratio_hist->Print();
        //std::cout << "Printing denominator histogram" << std::endl;
        //hData->Print();

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
        histInfo dummyR(new TH1D("dummyR", "dummyR", 1000, hData->GetBinLowEdge(1), hData->GetBinLowEdge(hData->GetNbinsX()) + hData->GetBinWidth(hData->GetNbinsX())));
        dummyR.setupAxesR(.4);
        dummyR.h->GetYaxis()->SetTitle("Tagging Rate");
        dummyR.h->GetXaxis()->SetTitle(xAxisLabel.c_str());

        //double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 1.5;
        //ymax = (ymax < 1.5 ? ymax : 2);
        dummyR.h->GetYaxis()->SetRangeUser(0.0, 2.0);

        //draw dummy axes
        dummyR.draw();

        //draw the ratio plot
        ratio_hist->Draw("same");

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        std::string fName = hNameDen + "SF.png";
        if(bkgSub) fName = "bkgSub_" + fName;
        c->Print(fName.c_str());

        std::string fNamePDF = hNameDen + "SF.pdf";
        if(bkgSub) fNamePDF = "bkgSub_" + fNamePDF;
        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
    }

    TH1D* plotOnlySF(const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
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
            MCEntries_[iBG].merge = merge;
            MCEntries_[iBG].retrieveHistogram();

            if(!hbgSumNum){
                hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                hbgSumNum->SetDirectory(0);
            }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

            //Get denominator histogram
            MCEntries_[iBG].histName = histNameDen;
            MCEntries_[iBG].rebin = rebin;
            MCEntries_[iBG].merge = merge;
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
                bgEntries_[iBG].merge = merge;
                bgEntries_[iBG].retrieveHistogram();

                if(!hNumSub){
                    hNumSub = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                    hNumSub->SetDirectory(0);
                }else        hNumSub->Add(bgEntries_[iBG].h.get());

                //Get denominator histogram
                bgEntries_[iBG].histName = histNameDen;
                bgEntries_[iBG].rebin = rebin;
                bgEntries_[iBG].merge = merge;
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
        data_.merge = merge;
        if(bkgSub) data_.subtract(hNumSub);
        data_.retrieveHistogram();
        hNumData = (static_cast<TH1D*>(data_.h.get()->Clone()));
        hNumData->Sumw2();
        hNumData->SetDirectory(0);
        //std::cout << "hNum info" << std::endl;
        //hNumData->Print();

        //data
        //get denominator histogram from file
        data_.histName = histNameDen;
        data_.rebin = rebin;
        data_.merge = merge;
        if(bkgSub) data_.subtract(hDenSub);
        data_.retrieveHistogram();
        hDenData = static_cast<TH1D*>(data_.h.get()->Clone());
        hDenData->Sumw2();
        hDenData->SetDirectory(0);
        //std::cout << "hDen info" << std::endl;
        //hDenData->Print();

        //std::cout << "Loaded data histograms" << std::endl;


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
        
        //std::cout << "Preparing ratio histogram" << std::endl;
        //hMC->ResetBit(kCanDelete);
        hMC->SetDirectory(0);
        //hMC->Print();
        ratio_hist = (TH1D*)(hData->Clone("ratio_hist"));
        ratio_hist->Sumw2();
        hData->Sumw2();

        //std::cout << "Printing numerator histogram" << std::endl;
        //ratio_hist->Print();
        //std::cout << "Printing denominator histogram" << std::endl;
        //hData->Print();

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
        dummy.h->GetYaxis()->SetRangeUser(0.,2.0);
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
        //double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 1.5;
        //ymax = (ymax < 1.5 ? ymax : 2);
        dummy.h->GetYaxis()->SetRangeUser(0.0, 2.0);


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

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        std::string fName = hNameDen + "OnlySF.png";
        if(bkgSub) fName = "bkgSub_"+fName;
        if(prefix.length() > 0) fName = prefix + "_" + fName;
        c->Print(fName.c_str());

        std::string fNamePDF = hNameDen + "OnlySF.pdf";
        if(bkgSub) fNamePDF = "bkgSub_"+fNamePDF;
        if(prefix.length() > 0) fNamePDF = prefix + "_" + fNamePDF;
        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;

        return ratio_hist;
    }

    TH1D* plotOnlyEff(bool data, const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
    {

        std::cout << "Called plotOnlyEff " << histNameNum << " " << histNameDen << std::endl;

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
                MCEntries_[iBG].merge = merge;
                MCEntries_[iBG].retrieveHistogram();

                if(!hbgSumNum){
                    hbgSumNum = static_cast<TH1D*>(MCEntries_[iBG].h->Clone());
                    hbgSumNum->SetDirectory(0);
                }else        hbgSumNum->Add(MCEntries_[iBG].h.get());

                //Get denominator histogram
                MCEntries_[iBG].histName = histNameDen;
                MCEntries_[iBG].rebin = rebin;
                MCEntries_[iBG].merge = merge;
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
                    bgEntries_[iBG].merge = merge;
                    bgEntries_[iBG].retrieveHistogram();

                    if(!hbgSumNum){
                        hbgSumNum = static_cast<TH1D*>(bgEntries_[iBG].h->Clone());
                        hbgSumNum->SetDirectory(0);
                    }else        hbgSumNum->Add(bgEntries_[iBG].h.get());

                    //Get denominator histogram
                    bgEntries_[iBG].histName = histNameDen;
                    bgEntries_[iBG].rebin = rebin;
                    bgEntries_[iBG].merge = merge;
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
            data_.merge = merge;
            if(bkgSub) data_.subtract(hbgSumNum);
            data_.retrieveHistogram();
            hNum = (static_cast<TH1D*>(data_.h.get()->Clone()));
            hNum->SetDirectory(0);
            //std::cout << "hNum info" << std::endl;
            //hNum->Print();

            //data
            //get denominator histogram from file
            data_.histName = histNameDen;
            data_.rebin = rebin;
            data_.merge = merge;
            if(bkgSub) data_.subtract(hbgSumDen);
            data_.retrieveHistogram();
            hDen = static_cast<TH1D*>(data_.h.get()->Clone());
            hDen->SetDirectory(0);
            //std::cout << "hDen info" << std::endl;
            //hDen->Print();

            //std::cout << "Loaded data histograms" << std::endl;

        }

        TH1D* ratio_hist;
        
        //std::cout << "Preparing ratio histogram" << std::endl;
        //hMC->ResetBit(kCanDelete);
        hDen->SetDirectory(0);
        //hMC->Print();
        ratio_hist = (TH1D*)(hNum->Clone("ratio_hist"));
        ratio_hist->Sumw2();
        hNum->Sumw2();

        //std::cout << "Printing numerator histogram" << std::endl;
        //ratio_hist->Print();
        //std::cout << "Printing denominator histogram" << std::endl;
        //hData->Print();

        ratio_hist->Divide(hDen); // data/MC
        ratio_hist->Sumw2();

        ratio_hist->SetLineColor(kBlack);
        ratio_hist->SetLineWidth(3);
        ratio_hist->SetMarkerColor(kBlack);
        ratio_hist->SetMarkerStyle(20);


        leg->AddEntry(ratio_hist, (data ? "Data Tagging Rate" : "MC Tagging Rate"), "PEX0");
        smartMax(ratio_hist, leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);
        //gPad->SetBottomMargin(0);



        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hNum->GetBinLowEdge(1), hNum->GetBinLowEdge(hNum->GetNbinsX()) + hNum->GetBinWidth(hNum->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetYaxis()->SetRangeUser(0.,2.0);
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        
        double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 1.5;
        ymax = (ymax < 1.5 ? ymax : 2);
        if(ymax < 1.9) std::cout << "The yaxis is allowed to be less than 2 in an efficiency plot." << std::endl;
        dummy.h->GetYaxis()->SetRangeUser(0.0, ymax);


        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        //gPad->SetLogy(isLogY);

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
        std::cout << "Drew Efficiency Only." << std::endl;
        ratio_hist->Print();
        std::cout << "Did the ratio_hist info print?" << std::endl;

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

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        std::string fName = hNameDen + "OnlyEff";
        if(data){ fName = fName + "Data.png"; }else{ fName = fName + "MC.png"; }
        if(bkgSub) fName = "bkgSub_"+fName;
        if(prefix.length() > 0) fName = prefix + "_" + fName;
        c->Print(fName.c_str());

        std::string fNamePDF = hNameDen + "OnlyEff";
        if(data){ fNamePDF = fNamePDF + "Data.pdf"; }else{ fNamePDF = fNamePDF + "MC.pdf"; }
        if(bkgSub) fNamePDF = "bkgSub_"+fNamePDF;
        if(prefix.length() > 0) fNamePDF = prefix + "_" + fNamePDF;
        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;

        return ratio_hist;
    }
};

int colorList[13] = { kRed+1, kBlack, kBlue+1, kSpring, kBlue+2, kOrange, kTeal, kCyan, kAzure, kPink, kViolet, kMagenta+3, kMagenta };

class MultiPlotter{
private:
    std::vector<Plotter*> plots;

public:
    void addPlot(Plotter* plot){
        plots.push_back(plot);
    }

    MultiPlotter(std::vector<Plotter*> tmpPlots){
        plots = tmpPlots;
    }

    //Loop over the Plotter objects and plot them in a single plot
    void plotOnlySF(const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
    {

        std::cout << "Called plotOnlySF (MultiPlotter)" << histNameNum << " " << histNameDen << std::endl;
        std::vector<TH1D*> sfHistos;
        TH1D *sfHisto;

        for(int i = 0; i < plots.size(); i++){
            sfHisto = plots[i]->plotOnlySF(histNameNum,histNameDen,xAxisLabel,yAxisLabel,isLogY,xmin,xmax,rebin,merge,lumi);

            sfHisto->SetLineColor(colorList[i]);
            sfHisto->SetLineWidth(3);
            sfHisto->SetMarkerColor(colorList[i]);
            sfHisto->SetMarkerStyle(20);

            sfHistos.push_back(sfHisto);

            std::cout << "Loaded histogram for " << plots[i]->getLabel() << std::endl;
        }


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

        for(int i = 0; i < sfHistos.size(); i++){
            leg->AddEntry(sfHistos[i], (plots[i]->getLabel()).c_str(), "PEX0");
            smartMax(sfHistos[i], leg, static_cast<TPad*>(gPad), min, max, lmax, true);
        }

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);


        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, sfHistos[0]->GetBinLowEdge(1), sfHistos[0]->GetBinLowEdge(sfHistos[0]->GetNbinsX()) + sfHistos[0]->GetBinWidth(sfHistos[0]->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetYaxis()->SetRangeUser(0.,2.0);
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        //Set the y-range of the histogram
        //double ymax = ratio_hist->GetBinContent(ratio_hist->GetMaximumBin()) * 1.5;
        //ymax = (ymax < 1.5 ? ymax : 2);
        dummy.h->GetYaxis()->SetRangeUser(0.0, 2.0);

        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //draw the histograms
        for(int i = 0; i < sfHistos.size(); i++){sfHistos[i]->Draw("same");}

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

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        std::string fName = hNameDen + "OnlySF.png";
        for(int i = 0; i < sfHistos.size(); i++){ fName = plots[i]->getPrefix()+"_"+fName; }
        fName = "Multi_" + fName;

        c->Print(fName.c_str());

        std::string fNamePDF = hNameDen + "OnlySF.pdf";
        for(int i = 0; i < sfHistos.size(); i++){ fNamePDF = plots[i]->getPrefix()+"_"+fNamePDF; }
        fNamePDF = "Multi_" + fNamePDF;

        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
    }

    //Loop over the Plotter objects and plot them in a single plot
    void plotOnlyEff(bool data, const std::string& histNameNum, const std::string& histNameDen, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double merge = -1000, double lumi = 36100)
    {

        std::cout << "Called plotOnlyEff (MultiPlotter)" << histNameNum << " " << histNameDen << std::endl;
        std::vector<TH1D*> sfHistos;
        TH1D *sfHisto;

        for(int i = 0; i < plots.size(); i++){
            sfHisto = plots[i]->plotOnlyEff(data,histNameNum,histNameDen,xAxisLabel,yAxisLabel,isLogY,xmin,xmax,rebin,merge,lumi);

            sfHisto->SetLineColor(colorList[i]);
            sfHisto->SetLineWidth(3);
            sfHisto->SetMarkerColor(colorList[i]);
            sfHisto->SetMarkerStyle(20);

            sfHistos.push_back(sfHisto);

            std::cout << "Loaded histogram for " << plots[i]->getLabel() << std::endl;
        }


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

        for(int i = 0; i < sfHistos.size(); i++){
            leg->AddEntry(sfHistos[i], (plots[i]->getLabel()).c_str(), "PEX0");
            smartMax(sfHistos[i], leg, static_cast<TPad*>(gPad), min, max, lmax, true);
        }

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas "c")
        std::cout << "Setting up pad" << std::endl;

        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);


        //create a dummy histogram to act as the axes
        std::cout << "Creating Dummy histogram" << std::endl;
        histInfo dummy(new TH1D("dummy", "dummy", 1000, sfHistos[0]->GetBinLowEdge(1), sfHistos[0]->GetBinLowEdge(sfHistos[0]->GetNbinsX()) + sfHistos[0]->GetBinWidth(sfHistos[0]->GetNbinsX())));
        std::cout << "Configuring dummy histogram" << std::endl;
        dummy.setupAxes(1);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetYaxis()->SetRangeUser(0.,2.0);
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        std::cout << "Configured dummy histogram" << std::endl;
        //Set the y-range of the histogram
        double ymax = 0;
        for(int i = 0; i < sfHistos.size(); i++){
            double ytest = sfHistos[i]->GetBinContent(sfHistos[i]->GetMaximumBin()) * 1.5;
            std::cout << plots[i]->getPrefix() << "suggests ymax of " << ytest << "." << std::endl;
            ymax = (ymax < ytest ? ytest : ymax);
        }
        ymax = (ymax < 1.5 ? ymax : 2);
        if(ymax < 1.9) std::cout << "For efficiency plots, we allow the y-axis to float" << std::endl;
        dummy.h->GetYaxis()->SetRangeUser(0.0, ymax);

        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        //gPad->SetLogy(isLogY);

        //Tick marks on both sides
        gPad->SetTicky();

        //Tick marks on both sides
        gPad->SetTicky();
        gPad->SetGridy();

        //draw the histograms
        std::cout << "Drawing the histograms" << std::endl;
        for(int i = 0; i < sfHistos.size(); i++){sfHistos[i]->Draw("same"); sfHistos[i]->Print(); std::cout << "Do we have to do this?" << std::endl;}

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

        //save new plot to file
        std::string hNameDen = (histNameDen.find("/") != std::string::npos ? histNameDen.substr(histNameDen.find("/")+1) : histNameDen);
        std::string fName = hNameDen + "OnlyEff";
        if(data){ fName = fName + "Data.png";}else{ fName = fName + "MC.png"; }
        for(int i = 0; i < sfHistos.size(); i++){ fName = plots[i]->getPrefix()+"_"+fName; }
        fName = "Multi_" + fName;

        c->Print(fName.c_str());

        std::string fNamePDF = hNameDen + "OnlyEff";
        if(data){ fNamePDF = fNamePDF + "Data.pdf";}else{ fNamePDF = fNamePDF + "MC.pdf"; }
        for(int i = 0; i < sfHistos.size(); i++){ fNamePDF = plots[i]->getPrefix()+"_"+fNamePDF; }
        fNamePDF = "Multi_" + fNamePDF;

        c->Print(fNamePDF.c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
    }
};

void makeplots(std::string CR, histInfo& data, std::vector<histInfo>& bgEntries, std::vector<histInfo>& wantedEntries, std::vector<histInfo>& unwantedEntries){
    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));
    Plotter pltSub(std::move(data), std::move(unwantedEntries), std::move(wantedEntries));

    plt.plot(CR+"/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/nJets", "Jets/Event", "Events", true);
    plt.plot(CR+"/nVertices", "primary vertices/event", "Events", true);
    plt.plot(CR+"/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot(CR+"/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot(CR+"/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot(CR+"/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotSF(CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF(CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF(CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF(CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF(CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF(CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, CR+"/METTagged",        CR+"/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/HTTagged",         CR+"/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/nJetsTagged",      CR+"/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/nVerticesTagged",  CR+"/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, CR+"/METTagged",        CR+"/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/HTTagged",         CR+"/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/nJetsTagged",      CR+"/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, CR+"/nVerticesTagged",  CR+"/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF(CR+"/METTagged",       CR+"/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF(CR+"/HTTagged",        CR+"/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF(CR+"/nJetsTagged",     CR+"/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF(CR+"/nVerticesTagged", CR+"/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF(CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF(CR+"/METTagged",       CR+"/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF(CR+"/HTTagged",        CR+"/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF(CR+"/nJetsTagged",     CR+"/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF(CR+"/nVerticesTagged", CR+"/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF(CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);
	
    pltSub.plot(CR+"/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/nJets", "Jets/Event", "Events", true);
    pltSub.plot(CR+"/nVertices", "primary vertices/event", "Events", true);
    pltSub.plot(CR+"/photon", "Photon p_{T} [GeV]", "Events", true);
    pltSub.plot(CR+"/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    pltSub.plot(CR+"/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot(CR+"/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    pltSub.plotSF(CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF(CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF(CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    pltSub.plotOnlySF(CR+"/bestTopPt", CR+"/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF(CR+"/bestTopMass", CR+"/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF(CR+"/bestTopEta", CR+"/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    pltSub.plotEff(true, CR+"/METTagged",        CR+"/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/HTTagged",         CR+"/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/nJetsTagged",      CR+"/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/nVerticesTagged",  CR+"/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(false, CR+"/METTagged",        CR+"/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/HTTagged",         CR+"/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/nJetsTagged",      CR+"/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(false, CR+"/nVerticesTagged",  CR+"/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotSF(CR+"/METTagged",       CR+"/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF(CR+"/HTTagged",        CR+"/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF(CR+"/nJetsTagged",     CR+"/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    pltSub.plotSF(CR+"/nVerticesTagged", CR+"/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF(CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    pltSub.plotOnlySF(CR+"/METTagged",       CR+"/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF(CR+"/HTTagged",        CR+"/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF(CR+"/nJetsTagged",     CR+"/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    pltSub.plotOnlySF(CR+"/nVerticesTagged", CR+"/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF(CR+"/photonTagged",     CR+"/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
