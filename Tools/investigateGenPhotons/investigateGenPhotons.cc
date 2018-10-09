// investigateGenPhotons.cc
// investigate gen photons
// Caleb J. Smith
// October 1, 2018

#include <iostream>
#include <string>
#include <vector>
#include "stdio.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TLegend.h"

// use this command to copy the required root file before running this script
// xrdcp root://xrootd.unl.edu//store/mc/RunIISummer16MiniAODv2/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_qcut19_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/FC894077-2DCA-E611-8008-002590DE6E3C.root .

TChain* combineTrees(std::vector<std::string> inputFiles, int n_lines)
{ 
    // create new TChain
    TChain *chain = new TChain("Events");

    for (const auto & infile : inputFiles)
    {
        std::ifstream in;
        in.open(infile.c_str());
        if(in.is_open()){
            std::cout << "- Successfully opened file: " << infile << std::endl;
        }
        else
        {
            std::cout << "- ERROR: Cannot open file: " << infile << std::endl;
            return nullptr;  
        }
          
        std::string line;
        int i = 0;
        while(in.good()){
            if (i >= n_lines) break;                // only read n_lines
            if( ! std::getline(in,line) ) break;    // read a line from the file
            if( chain->Add( line.c_str() ) )        // add to the chain 
            {
                std::cout << "- Loaded tree " << i <<  ": " << line << std::endl;
            }
            else
            {
                std::cout << "- Problem loading tree " << i <<  ": " << line << std::endl;
            }
            // increase line index
            i++; 
        }
          
        in.close();
    }

    // return TChain
    return chain;
}

// simple class
class DrawOptions
{
    public:
        std::string plotName;
        std::string varexp;
        std::string selection;
        int color;
};

void plot(TChain* chain, DrawOptions p)
{
    printf("- Creating plot %s\n", p.plotName.c_str());
    std::string plotDir = "plots/";
    TCanvas* c1 = new TCanvas("c","c", 600.0, 600.0);
    TH1F* htemp = nullptr;
    
    // draw
    chain->SetLineColor(p.color);
    chain->Draw(p.varexp.c_str(), p.selection.c_str(), "hist E");
    // use htemp to set x-axis range
    // chain->Draw() creates htemp
    // this must be done after chain->Draw()
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    if (htemp)
    {
        htemp->GetXaxis()->SetRangeUser(-5.0, 5.0);
        htemp->SetTitle(p.plotName.c_str());
    }
    else
    {
        printf("ERROR: htemp is nullptr\n");
        exit(1);
    }
    
    c1->Update();
    c1->SaveAs((plotDir + p.plotName + ".png").c_str());
    c1->SaveAs((plotDir + p.plotName + ".pdf").c_str());
    delete c1;
}

void multiplot(TChain* chain, std::vector<DrawOptions> vp, std::string plotName, bool normalize=false)
{
    printf("- Creating multiplot %s\n", plotName.c_str());
    std::string plotDir = "plots/";
    TCanvas* c1 = new TCanvas("c","c", 600.0, 600.0);
    TLegend* legend = new TLegend(0.68, 0.55, 0.98, 0.75); // x1, y1, x2, y2
    legend->SetHeader("prunedGenParticles","C"); // option "C" allows to center the header
    int i = 0;
    double h_max = -999.0;
    double h_integral = -999.0;
    TH1F* hfirst = nullptr;
    for (const auto & p : vp)
    {
        // draw
        chain->SetLineColor(p.color);
        if (i == 0)
        {
            chain->Draw(p.varexp.c_str(), p.selection.c_str(), "hist E");
        }
        else
        {
            chain->Draw(p.varexp.c_str(), p.selection.c_str(), "same hist E");
        }
        // this must be done after chain->Draw()
        // use htemp to set x-axis range and to put in legend with correct color
        // chain->Draw() creates a new htemp every time
        // every time we do chain->Draw(), we get another htemp
        // use GetListOfPrimitives and index to get the correct htemp
        
        // print what's available
        printf("plotName: %s, i: %d, color: %d\n", p.plotName.c_str(), i, p.color);
        gPad->GetListOfPrimitives()->Print();
        
        TH1F* htemp = nullptr;
        htemp = (TH1F*)gPad->GetListOfPrimitives()->At(i);
        if (htemp)
        {
            // save first histogram to set ranges
            if (i == 0) hfirst = htemp;
            legend->AddEntry(htemp, p.plotName.c_str());
            htemp->SetTitle(plotName.c_str());
            // get max integral
            // do this before applying normalization
            h_integral = std::max(h_integral, htemp->Integral());
            if (normalize)
            {
                htemp->Scale(1.0/htemp->Integral());
            }
            // do this after applying normalization
            // note that htemp->GetMaximum() only give unscale max
            // this method gives the max after scaling
            h_max = std::max(h_max, htemp->GetBinContent(htemp->GetMaximumBin()));
        }
        else
        {
            printf("ERROR: htemp is nullptr\n");
            exit(1);
        }
        i++;
    }

    printf("h_max = %f\n", h_max);
    printf("h_integral = %f\n", h_integral);
    hfirst->GetXaxis()->SetRangeUser(-5.0, 5.0);
    hfirst->GetYaxis()->SetRangeUser(0.0, 1.25*h_max);

    legend->Draw("hist E");
    c1->Update();
    c1->SaveAs((plotDir + plotName + ".png").c_str());
    c1->SaveAs((plotDir + plotName + ".pdf").c_str());
    delete c1;
}

void investigate(std::vector<std::string> inputFiles)
{
    TChain* chain = combineTrees(inputFiles, 2); 
    if (!chain)
    {
        printf("ERROR: Did not create chain from input files\n");
        exit(1);
    }
    
    // make plots
    
    std::vector<DrawOptions> plotOptions;
    std::string plotName = "";
    std::string varexp = "";
    std::string selection = "";
    
    plotName = "genEta_eta<5";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0";
    DrawOptions p1 = {plotName, varexp, selection, kRed};
    plotOptions.push_back(p1);
    
    plotName = "genEta_eta<5_pt>10";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0";
    DrawOptions p2 = {plotName, varexp, selection, kBlue};
    plotOptions.push_back(p2);
    
    plotName = "genEta_eta<5_pt>10_pdgId=22";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0 && \
                 recoGenParticles_prunedGenParticles__PAT.obj.pdgId() == 22";
    DrawOptions p3 = {plotName, varexp, selection, kGreen-2};
    plotOptions.push_back(p3);
    
    plotName = "genEta_eta<5_pt>10_pdgId=22_status=(1||2||20--29)";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0 && \
                     recoGenParticles_prunedGenParticles__PAT.obj.pdgId() == 22 && \
                    (recoGenParticles_prunedGenParticles__PAT.obj.status() == 1 || \
                     recoGenParticles_prunedGenParticles__PAT.obj.status() == 2 || \
                    (recoGenParticles_prunedGenParticles__PAT.obj.status())/10 == 2)";
    DrawOptions p4 = {plotName, varexp, selection, kMagenta+1};
    plotOptions.push_back(p4);
    
    // loop over plots
    for (const auto& p : plotOptions)
    {
        plot(chain, p);
    }
    
    // multiplot (original)
    multiplot(chain, plotOptions, "gen_with_cuts_original");
    // multiplot (normalized)
    multiplot(chain, plotOptions, "gen_with_cuts_normalized", true);

    delete chain;
}

int main()
{
    std::vector<std::string> inputFiles;
    inputFiles.push_back("RunIISummer16MiniAODv2/GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt");
    inputFiles.push_back("RunIISummer16MiniAODv2/GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt");
    inputFiles.push_back("RunIISummer16MiniAODv2/GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt");
    inputFiles.push_back("RunIISummer16MiniAODv2/GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt");
    inputFiles.push_back("RunIISummer16MiniAODv2/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt");
    investigate(inputFiles);
}




