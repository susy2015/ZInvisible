// script to copy a subset of a Tree to a new Tree

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>


int main(int argc, char *argv[])
{
    int verbose = 1;
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
    Int_t oldnentries = (Int_t)oldtree->GetEntries();
    
    // branches
    ULong64_t   event              = 0;
    int         Stop0l_nResolved   = 0;
    oldtree->SetBranchAddress("event", &event);
    oldtree->SetBranchAddress("Stop0l_nResolved", &Stop0l_nResolved);

    // //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile("TTGJets_2018_0_skimmed.root","recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    //for (Int_t i=0;i<1000; i++) {
    for (Int_t i=0;i<oldnentries; i++) {
        oldtree->GetEntry(i);
        if (Stop0l_nResolved >= 2)
        if (true)
        {
            newtree->Fill();
            if (verbose > 1) printf("i = %d, event = %d, Stop0l_nResolved = %d\n", i, event, Stop0l_nResolved);
        }
    }
    
    Int_t newnentries = (Int_t)newtree->GetEntries();
    
    if (verbose > 0)
    {
        std::cout << "nEvents before skimming: " << oldnentries << std::endl;
        std::cout << "nEvents after skimming: "  << newnentries << std::endl;
    }
    
    //newtree->Print();
    newtree->AutoSave();
    delete oldfile;
    delete newfile;
}
