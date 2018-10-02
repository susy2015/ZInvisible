// investigate.cc
// investigate gen photons
// Caleb J. Smith
// October 1, 2018

#include <string>
#include <vector>
#include "stdio.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"

void investigate(const char* inputFileName);
void plot(TTree* tree, std::string plotName, const char * varexp, const char * selection);

void plot(TTree* tree, std::string plotName, const char * varexp, const char * selection)
{
    std::string plotDir = "plots/";
    TCanvas* c1 = new TCanvas("c","c", 600.0, 600.0);
    tree->Draw(varexp, selection);
    c1->SaveAs((plotDir + plotName + ".png").c_str());
    c1->SaveAs((plotDir + plotName + ".pdf").c_str());
    delete c1;
}

void investigate(const char* inputFileName)
{
    // open that file please 
    TFile* file = NULL;
    file = TFile::Open(inputFileName);
    if (!file)
    {
        printf("ERROR: Did not open file %s\n", inputFileName);
        exit(1);
    }
    // get tree
    const char* treeName = "Events";
    TTree* tree = NULL;
    tree = (TTree*)file->Get(treeName);
    if (!tree)
    {
        printf("ERROR: Did not open tree %s\n", treeName);
        exit(1);
    }
    // make plots
    std::string plotName = "";
    const char * varexp = "";
    const char * selection = "";
    
    plotName = "genStatus_pdgId=22_pt>10";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.status()";
    selection = "recoGenParticles_prunedGenParticles__PAT.obj.pdgId() == 22 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0";
    plot(tree, plotName, varexp, selection);
    
    plotName = "genEta_eta<5";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0";
    plot(tree, plotName, varexp, selection);
    
    plotName = "genEta_eta<5_pt>10";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0";
    plot(tree, plotName, varexp, selection);
    
    plotName = "genEta_eta<5_pt>10_pdgId=22";
    varexp = "recoGenParticles_prunedGenParticles__PAT.obj.eta()";
    selection = "abs(recoGenParticles_prunedGenParticles__PAT.obj.eta()) < 5.0 && recoGenParticles_prunedGenParticles__PAT.obj.pt() > 10.0 && recoGenParticles_prunedGenParticles__PAT.obj.pdgId() == 22";
    plot(tree, plotName, varexp, selection);
    
    // close that file please
    file->Close();
}

int main()
{
    const char* inputFileName = "FC894077-2DCA-E611-8008-002590DE6E3C.root"; 
    investigate(inputFileName);
}




