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

void investigate(const char* fileName);

void investigate(const char* fileName)
{
    // open file 
    TFile* file = NULL;
    file = TFile::Open(fileName);
    if (!file)
    {
        printf("ERROR: Did not open file %s\n", fileName);
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
    TCanvas* c1 = new TCanvas("c","c", 600.0, 600.0);

    delete c1;
    // close file 
    file->Close();
}

int main()
{
    const char* fileName = "FC894077-2DCA-E611-8008-002590DE6E3C.root"; 
    investigate(fileName);
}


