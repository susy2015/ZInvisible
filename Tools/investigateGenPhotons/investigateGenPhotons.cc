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

TChain* chainFiles(const char *infile, int n_lines)
{ 
    std::ifstream in;
    in.open(infile);
    if(!in.is_open()){
        std::cout << "Cannot open list file: " << infile << std::endl;
        return nullptr;  
    }

    TChain *chain = new TChain("Events");
      
    std::string line;
    int i = 0;
    while(in.good()){
        if (i >= n_lines) break;                // only read n_lines
        if( ! std::getline(in,line) ) break;    // read a line from the file
        if( chain->Add( line.c_str() ) )        // add to the chain 
        {
            std::cout << "Added file " << i <<  ": " << line << std::endl;
        }
        else
        {
            std::cout << "Problem loading tree from " << line << std::endl;
        }
        // increase line index
        i++; 
    }
      
    in.close();

    //chain->Merge(outfile);
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

void plot(TTree* tree, DrawOptions p)
{
    printf("- Creating plot %s\n", p.plotName.c_str());
    std::string plotDir = "plots/";
    TCanvas* c1 = new TCanvas("c","c", 600.0, 600.0);
    TH1F* htemp = nullptr;
    
    // draw
    tree->SetLineColor(p.color);
    tree->Draw(p.varexp.c_str(), p.selection.c_str(), "hist E");
    // use htemp to set x-axis range
    // tree->Draw() creates htemp
    // this must be done after tree->Draw()
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

void multiplot(TTree* tree, std::vector<DrawOptions> vp, std::string plotName, bool normalize=false)
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
        tree->SetLineColor(p.color);
        if (i == 0)
        {
            tree->Draw(p.varexp.c_str(), p.selection.c_str(), "hist E");
        }
        else
        {
            tree->Draw(p.varexp.c_str(), p.selection.c_str(), "same hist E");
        }
        // this must be done after tree->Draw()
        // use htemp to set x-axis range and to put in legend with correct color
        // tree->Draw() creates a new htemp every time
        // every time we do tree->Draw(), we get another htemp
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

void investigate(const char* inputFileName)
{
    TChain* chain = chainFiles(inputFileName, 10); 
    if (!chain)
    {
        printf("ERROR: Did not create chain from %s\n", inputFileName);
        exit(1);
    }
    // // open that file please 
    // TFile* file = nullptr;
    // file = TFile::Open(inputFileName);
    // if (!file)
    // {
    //     printf("ERROR: Did not open file %s\n", inputFileName);
    //     exit(1);
    // }
    // get tree
    std::string treeName = "Events";
    TTree* tree = nullptr;
    std::cout << "Number of trees: " << chain->GetNtrees() << std::endl;
    std::cout << "Find branch: " << chain->FindBranch(treeName.c_str()) << std::endl;
    tree = (TTree*)chain->GetTree();
    //tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
    {
        printf("ERROR: Did not open tree %s\n", treeName.c_str());
        exit(1);
    }
    // make plots
    std::vector<DrawOptions> plotOptions;
    std::string plotName = "";
    std::string varexp = "";
    std::string selection = "";
    
    // from Dr. Ken
    // plotName = "genStatus_pdgId=22_pt>10";
    // varexp = "recoGenParticles_prunedGenParticles__PAT.obj.status()";
    // selection = "recoGenParticles_prunedGenParticles__PAT.obj.pdgId() == 22 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0";
    // DrawOptions p1 = {plotName, varexp, selection};
    // plotOptions.push_back(p1);
    
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
        plot(tree, p);
    }
    
    // multiplot (original)
    multiplot(tree, plotOptions, "gen_with_cuts_original");
    // multiplot (normalized)
    multiplot(tree, plotOptions, "gen_with_cuts_normalized", true);

    // close that file please
    //file->Close();
    // delete file
    //delete file;
    delete chain;
}

int main()
{
    //const char* inputFileName = "FC894077-2DCA-E611-8008-002590DE6E3C.root"; 
    const char* inputFileName = "RunIISummer16MiniAODv2/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"; 
    investigate(inputFileName);
}




