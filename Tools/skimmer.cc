// script to copy a subset of a Tree to a new Tree

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>


int main(int argc, char *argv[])
{
    //gROOT->Reset();
    
    // TTGJets_2018 files:
    // root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_0.root
    // root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_1.root
    // root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_2.root
    // root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_3.root
    // root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_4.root
    
    // Tree: Events
    
    //Get old file, old tree and set top branch address
    TFile *oldfile = TFile::Open("root://cmseos.fnal.gov//store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PostProcessed_22March2019_v6/TTGJets_2018/TTGJets_2018_0.root");
    TTree *oldtree = (TTree*)oldfile->Get("Events");
    Int_t nentries = (Int_t)oldtree->GetEntries();
    
    std::cout << "nEvents before skimming: " << nentries << std::endl;
    
    // branches
    ULong64_t *event = 0;
    oldtree->SetBranchAddress("event", &event);

    // //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile("TTGJets_2018_0_skimmed.root","recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    //for (Int_t i=0;i<nentries; i++) {
    for (Int_t i=0;i<1000; i++) {
       oldtree->GetEntry(i);
       printf("i = %d, event = %d\n", i, event);
    }
    
    //newtree->Print();
    newtree->AutoSave();
    delete oldfile;
    delete newfile;
}
