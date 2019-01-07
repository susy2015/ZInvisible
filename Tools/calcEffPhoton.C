#include "SusyAnaTools/Tools/samples.h"
#include "RegisterFunctions.h"
#include "NTupleReader.h"

#include <iostream>
#include <getopt.h>

#include "TVectorD.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TChain.h"
#include "Math/VectorUtil.h"

int main(int argc, char* argv[])
{
    using namespace std;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"condor",           no_argument, 0, 'c'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
    };

    bool runOnCondor = false;
    int nFiles = -1, startFile = 0, nEvts = -1;
    string dataSets = "";

    while((opt = getopt_long(argc, argv, "cD:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'c':
            runOnCondor = true;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            nEvts = int(atoi(optarg));
            break;
        }
    }

    string filename = "effhists.root";
    string sampleloc = AnaSamples::fileDir;

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "effhists_%s_%d.root", dataSets.c_str(), startFile);
        filename = thistFile;
        sampleloc = "condor";
    }

  
    // follow this syntax; order matters for your arguments
    
    //SampleSet::SampleSet(std::string file, bool isCondor, double lumi)
    AnaSamples::SampleSet        ss("sampleSets.cfg", runOnCondor, AnaSamples::luminosity);
    
    //SampleCollection::SampleCollection(const std::string& file, SampleSet& samples) : ss_(samples)
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);
    double xsec_ratio = -1.0;
    
    if (ss[dataSets] != ss.null())
    {
        std::vector<std::string> sampleTags1 = {"GJets_HT-200To400", "GJets_HT-400To600", "GJets_HT-600ToInf"};
        std::vector<std::string> sampleTags2 = {"ZJetsToNuNu_HT_200to400", "ZJetsToNuNu_HT_400to600", "ZJetsToNuNu_HT_600to800", "ZJetsToNuNu_HT_800to1200", "ZJetsToNuNu_HT_1200to2500", "ZJetsToNuNu_HT_2500toInf"};
        xsec_ratio = ss.getCrossSectionRatio(sampleTags1, sampleTags2);
    }
    else if (sc[dataSets] != sc.null())
    {
        std::string sampleTag1 = "GJets";
        std::string sampleTag2 = "ZJetsToNuNu";
        xsec_ratio = sc.getCrossSectionRatio(sampleTag1, sampleTag2);
    }


    TFile *f = new TFile(filename.c_str(),"RECREATE");
    f->cd();

    // cross section ratio
    TH1 *hCrossSectionRatio = new TH1D("hCrossSectionRatio", "hCrossSectionRatio", 1, 0, 1);
    hCrossSectionRatio->Fill(0.5, xsec_ratio);
    // photons
    TH1 *hPhotonAccPt_num = new TH1D("hPhotonAccPt_num", "hPhotonAccPt_num", 200, 0, 2000);
    TH1 *hPhotonAccPt_den = new TH1D("hPhotonAccPt_den", "hPhotonAccPt_den", 200, 0, 2000);
    TH1 *hPhotonEffPt_num = new TH1D("hPhotonEffPt_num", "hPhotonEffPt_num", 200, 0, 2000);
    TH1 *hPhotonEffPt_den = new TH1D("hPhotonEffPt_den", "hPhotonEffPt_den", 200, 0, 2000);
    
    RegisterFunctionsCalcEff rt;

    std::set<std::string> activeBranches;
    rt.activateBranches(activeBranches);

    // --- new method: testing ---
    map<string, vector<AnaSamples::FileSummary>> fileMap;
    
    std::cout << "dataSets = " << dataSets << std::endl;
    
    if(ss[dataSets] != ss.null())
    {
        fileMap[dataSets] = {ss[dataSets]};
        for(const auto& colls : ss[dataSets].getCollections())
        {
            fileMap[colls] = {ss[dataSets]};
        }
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = {sc[dataSets]};
        int i = 0;
        for(const auto& fs : sc[dataSets])
        {
            fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }
    
    set<AnaSamples::FileSummary> setFS;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) setFS.insert(fs);
    
    std::cout << "Running over files..." << std::endl;
    int printInterval = 1000;
    for(const AnaSamples::FileSummary& file : setFS)
    {
        int fileCount = 0, startCount = 0;
        int NEvtsTotal = 0;
        file.readFileList(); //get file list
        for(const std::string& fname : file.filelist_)
        {
            if(startCount++ < startFile) continue;
            if(nFiles > 0 && fileCount++ >= nFiles) break;
            if(nFiles > 0) NEvtsTotal = 0;
            else if(nEvts >= 0 && NEvtsTotal > nEvts) break;
            TFile *currentFile = TFile::Open(fname.c_str());
            
            TTree *t = (TTree*)currentFile->Get(file.treePath.c_str());
            if (!t) 
            {
                std::cout << "ERROR: TTree " << file.treePath << " not found in file " << filename << std::endl;
                currentFile->Close();
                delete currentFile;
                f->Close();
                delete f;
                return 0;
            }
            try
            {
                NTupleReader tr(t, activeBranches);
                tr.setReThrow(false);
                rt.registerFunctions(tr);
                while(tr.getNextEvent())
                {
                    if(nEvts > 0 && NEvtsTotal > nEvts) break;
                    if(tr.getEvtNum() % printInterval == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;
                    const auto& gammaLVecReco        = tr.getVec<TLorentzVector>("gammaLVecReco");
                    const auto& gammaLVecRecoEtaPt   = tr.getVec<TLorentzVector>("gammaLVecRecoEtaPt");
                    const auto& gammaLVecPassLooseID = tr.getVec<TLorentzVector>("gammaLVecPassLooseID");
                    for(auto& tlv : gammaLVecReco)
                    {
                        hPhotonAccPt_den->Fill(tlv.Pt(), file.getWeight());
                    }
                    for(auto& tlv : gammaLVecRecoEtaPt)
                    {
                        hPhotonAccPt_num->Fill(tlv.Pt(), file.getWeight());
                        hPhotonEffPt_den->Fill(tlv.Pt(), file.getWeight());
                    }
                    for(auto& tlv : gammaLVecPassLooseID)
                    {
                        hPhotonEffPt_num->Fill(tlv.Pt(), file.getWeight());
                    }
                    ++NEvtsTotal;
                }
            }
            catch(const std::string e)
            {
                std::cout << "ERROR: Exception when creating NTupleReader: " << e << std::endl;
            }
    
            t->Reset();
            delete t;
            currentFile->Close();
            delete currentFile;
        }
    }

    
    f->cd();
    
    hCrossSectionRatio->Write();
    hPhotonAccPt_num->Write();
    hPhotonAccPt_den->Write();
    hPhotonEffPt_num->Write();
    hPhotonEffPt_den->Write();
    
    std::cout << "Quitting..." << std::endl;
    f->Close();
    return 0;

}
