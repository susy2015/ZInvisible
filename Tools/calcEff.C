#include "SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

int main()
{
    AnaSamples::SampleSet        ss("/eos/uscms/store/user/lpcsusyhad/PHYS14_720_Dec23_2014/");
    AnaSamples::SampleCollection sc(ss);

    TFile *f = new TFile("effhists.root","RECREATE");
    f->cd();

    TH1 *hMuEffPt_num = new TH1D("hMuEffPt_num", "hMuEffPt_num", 200, 0, 2000);
    TH1 *hMuEffPt_den = new TH1D("hMuEffPt_den", "hMuEffPt_den", 200, 0, 2000);
    TH1 *hMuAccPt_num = new TH1D("hMuAccPt_num", "hMuAccPt_num", 200, 0, 2000);
    TH1 *hMuAccPt_den = new TH1D("hMuAccPt_den", "hMuAccPt_den", 200, 0, 2000);

    TH1 *hMuEffHt_num = new TH1D("hMuEffHt_num", "hMuEffHt_num", 300, 0, 3000);
    TH1 *hMuEffHt_den = new TH1D("hMuEffHt_den", "hMuEffHt_den", 300, 0, 3000);
    TH1 *hMuAccHt_num = new TH1D("hMuAccHt_num", "hMuAccHt_num", 300, 0, 3000);
    TH1 *hMuAccHt_den = new TH1D("hMuAccHt_den", "hMuAccHt_den", 300, 0, 3000);

    TH2 *hMuEff_num = new TH2D("hMuEff_num", "hMuEff_num", 200, 0, 2000, 300, 0, 3000);
    TH2 *hMuEff_num_pp = new TH2D("hMuEff_num_pp", "hMuEff_num_pp", 200, 0, 2000, 300, 0, 3000);
    TH2 *hMuEff_num_pf = new TH2D("hMuEff_num_pf", "hMuEff_num_pf", 200, 0, 2000, 300, 0, 3000);
    TH2 *hMuEff_num_fp = new TH2D("hMuEff_num_fp", "hMuEff_num_fp", 200, 0, 2000, 300, 0, 3000);

    TH2 *hMuEff_num_rand = new TH2D("hMuEff_num_rand", "hMuEff_num_rand", 200, 0, 2000, 300, 0, 3000);

    TH2 *hMuEff_den = new TH2D("hMuEff_den", "hMuEff_den", 200, 0, 2000, 300, 0, 3000);
    TH2 *hMuAcc_num = new TH2D("hMuAcc_num", "hMuAcc_num", 200, 0, 2000, 300, 0, 3000);
    TH2 *hMuAcc_den = new TH2D("hMuAcc_den", "hMuAcc_den", 200, 0, 2000, 300, 0, 3000);

    std::set<std::string> activeBranches;

    activeBranches.insert("ht");
    activeBranches.insert("genDecayPdgIdVec");
    activeBranches.insert("genDecayLVecOB");
    activeBranches.insert("genDecayLVec");
    activeBranches.insert("muonsLVec");
    activeBranches.insert("muonsRelIso");
    activeBranches.insert("W_emuVec");
    activeBranches.insert("muonsCharge");
    activeBranches.insert("muonsMtw");
    activeBranches.insert("met");
    activeBranches.insert("metphi");
    activeBranches.insert("jetsLVec");
    activeBranches.insert("elesLVec");
    activeBranches.insert("elesRelIso");
    activeBranches.insert("recoJetsBtag_0");

    TRandom3 *trg = new TRandom3(12321);

    for(auto& file : sc["DYJetsToLL"]) 
    {
        TChain *t = new TChain(sc["DYJetsToLL"].front().treePath.c_str());

        file.addFilesToChain(t);

        std::cout << "Processing file(s): " << file.filePath << std::endl;

        NTupleReader tr(t, activeBranches);
        tr.registerFunction(&plotterFunctions::cleanJets);
        tr.registerFunction(&plotterFunctions::muInfo);

        while(tr.getNextEvent())
        {
            const std::vector<const TLorentzVector*>& genMatchMuInAcc = tr.getVec<const TLorentzVector*>("genMatchMuInAcc");
            const std::vector<const TLorentzVector*>& genMuInAcc      = tr.getVec<const TLorentzVector*>("genMuInAcc");
            const std::vector<const TLorentzVector*>& genMu           = tr.getVec<const TLorentzVector*>("genMu");

            const double& cleanHt = tr.getVar<double>("ht");
        
            for(auto& tlv : genMu)
            {
                hMuAccPt_den->Fill(tlv->Pt(), file.getWeight());
                hMuAccHt_den->Fill(cleanHt,   file.getWeight());

                hMuAcc_den->Fill(tlv->Pt(), cleanHt, file.getWeight());
            }


            int count = 0, random = trg->Integer(400000000) & 1;
            bool oneMatch = false, twoMatch = false;

            for(auto& tlv : genMuInAcc)
            {
                hMuAccPt_num->Fill(tlv->Pt(), file.getWeight());
                hMuAccHt_num->Fill(cleanHt,   file.getWeight());

                hMuAcc_num->Fill(tlv->Pt(), cleanHt, file.getWeight());

                hMuEffPt_den->Fill(tlv->Pt(), file.getWeight());
                hMuEffHt_den->Fill(cleanHt,   file.getWeight());

                hMuEff_den->Fill(tlv->Pt(), cleanHt, file.getWeight());

                for(auto& tlv2 : genMatchMuInAcc)
                {
                    if(     count == 0 && tlv == tlv2) oneMatch = true;
                    else if(count == 1 && tlv == tlv2) twoMatch = true;

                    if(count == random && tlv == tlv2) hMuEff_num_rand->Fill(tlv->Pt(), cleanHt, file.getWeight());
                }
                count++;
            }

            for(auto& tlv : genMatchMuInAcc)
            {
                hMuEffPt_num->Fill(tlv->Pt(), file.getWeight());
                hMuEffHt_num->Fill(cleanHt,   file.getWeight());

                hMuEff_num->Fill(tlv->Pt(), cleanHt, file.getWeight());

                if(      oneMatch && twoMatch) hMuEff_num_pp->Fill(tlv->Pt(), cleanHt, file.getWeight());
                else if(!oneMatch && twoMatch) hMuEff_num_fp->Fill(tlv->Pt(), cleanHt, file.getWeight());
                else if(oneMatch && !twoMatch) hMuEff_num_pf->Fill(tlv->Pt(), cleanHt, file.getWeight());
            }
        }
    }

    hMuEffPt_num->Write();
    hMuEffPt_den->Write();
    hMuAccPt_num->Write();
    hMuAccPt_den->Write();
    hMuEffHt_num->Write();
    hMuEffHt_den->Write();
    hMuAccHt_num->Write();
    hMuAccHt_den->Write();
    hMuEff_num->Write();
    hMuEff_den->Write();
    hMuAcc_num->Write();
    hMuAcc_den->Write();

    hMuEff_num_pp->Write();
    hMuEff_num_pf->Write();
    hMuEff_num_fp->Write();

    hMuEff_num_rand->Write();

    f->Close();
}
