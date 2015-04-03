#include "SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

int main()
{
    AnaSamples::SampleSet        ss;
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

    TH1 *hZEffPt_num = new TH1D("hZEffPt_num", "hZEffPt_num", 200, 0, 2000);
    TH1 *hZEffPt_den = new TH1D("hZEffPt_den", "hZEffPt_den", 200, 0, 2000);
    TH1 *hZAccPt_num = new TH1D("hZAccPt_num", "hZAccPt_num", 200, 0, 2000);
    TH1 *hZAccPt_den = new TH1D("hZAccPt_den", "hZAccPt_den", 200, 0, 2000);

    TH2 *hZEff_num = new TH2D("hZEff_num", "hZEff_num", 200, 0, 2000, 300, 0, 3000);
    TH2 *hZEff_den = new TH2D("hZEff_den", "hZEff_den", 200, 0, 2000, 300, 0, 3000);
    TH2 *hZAcc_num = new TH2D("hZAcc_num", "hZAcc_num", 200, 0, 2000, 300, 0, 3000);
    TH2 *hZAcc_den = new TH2D("hZAcc_den", "hZAcc_den", 200, 0, 2000, 300, 0, 3000);

    TH2 *hZEff_jActR1_num = new TH2D("hZEff_jActR1_num", "hZEff_jActR1_num", 200, 0, 2000, 500, 0, 5000);
    TH2 *hZEff_jActR1_den = new TH2D("hZEff_jActR1_den", "hZEff_jActR1_den", 200, 0, 2000, 500, 0, 5000);
    TH2 *hZEff_jActR2_num = new TH2D("hZEff_jActR2_num", "hZEff_jActR2_num", 200, 0, 2000, 500, 0, 5000);
    TH2 *hZEff_jActR2_den = new TH2D("hZEff_jActR2_den", "hZEff_jActR2_den", 200, 0, 2000, 500, 0, 5000);

    std::set<std::string> activeBranches;

    activeBranches.insert("ht");
    activeBranches.insert("genDecayPdgIdVec");
    activeBranches.insert("genDecayLVecOB");
    activeBranches.insert("genDecayLVec");
    activeBranches.insert("muonsLVec");
    activeBranches.insert("muonsRelIso");
    activeBranches.insert("muonsMiniIso");
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
        tr.registerFunction(&plotterFunctions::generateWeight);

        while(tr.getNextEvent())
        {
            const std::vector<const TLorentzVector*>& genMatchMuInAcc = tr.getVec<const TLorentzVector*>("genMatchMuInAcc");
            const std::vector<const TLorentzVector*>& genMuInAcc      = tr.getVec<const TLorentzVector*>("genMuInAcc");
            const std::vector<const TLorentzVector*>& genMu           = tr.getVec<const TLorentzVector*>("genMu");

            const double& recoZPt = tr.getVar<double>("bestRecoZPt");
            const double& genZPt  = tr.getVar<double>("genZPt");
            const double& cleanHt = tr.getVar<double>("ht");
            const int&    pdgIdZDec = tr.getVar<int>("pdgIdZDec");

            if(pdgIdZDec != 13) continue;

            for(auto& tlv : genMu)
            {
                hMuAccPt_den->Fill(tlv->Pt(), file.getWeight());
                hMuAccHt_den->Fill(cleanHt,   file.getWeight());

                hMuAcc_den->Fill(tlv->Pt(), cleanHt, file.getWeight());
            }

            int count = 0, random = 0;//trg->Integer(400000000) & 1;
            bool oneMatch = false, twoMatch = false;

            double modHt = cleanHt;

            for(auto& tlv : genMuInAcc)
            {
                if(tlv->Pt() > 50) modHt -= tlv->Pt();
            }

            const std::vector<double>& jActR1 = tr.getVec<double>("jActR1");
            const std::vector<double>& jActR2 = tr.getVec<double>("jActR2");

            for(auto& tlv : genMuInAcc)
            {
                hMuAccPt_num->Fill(tlv->Pt(), file.getWeight());
                hMuAccHt_num->Fill(cleanHt,   file.getWeight());

                hMuAcc_num->Fill(tlv->Pt(), cleanHt, file.getWeight());

                hMuEffPt_den->Fill(tlv->Pt(), file.getWeight());
                hMuEffHt_den->Fill(cleanHt,   file.getWeight());

                hMuEff_den->Fill(tlv->Pt(), cleanHt, file.getWeight());

                hZEff_jActR1_den->Fill(tlv->Pt(), jActR1[count], file.getWeight());
                hZEff_jActR2_den->Fill(tlv->Pt(), jActR2[count], file.getWeight());

                for(auto& tlv2 : genMatchMuInAcc)
                {
                    if(     count == 0 && tlv == tlv2) oneMatch = true;
                    else if(count == 1 && tlv == tlv2) twoMatch = true;

                    if(count == random && tlv == tlv2) hMuEff_num_rand->Fill(tlv->Pt(), cleanHt, file.getWeight());

                    if(tlv == tlv2)
                    {
                        hZEff_jActR1_num->Fill(tlv->Pt(), jActR1[count], file.getWeight());
                        hZEff_jActR2_num->Fill(tlv->Pt(), jActR2[count], file.getWeight());
                    }
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


            hZAccPt_den->Fill(genZPt, file.getWeight());
            hZAcc_den->Fill(genZPt, modHt, file.getWeight());
            if(genMuInAcc.size() >= 2)
            {

                hZAccPt_num->Fill(genZPt, file.getWeight());
                hZAcc_num->Fill(genZPt, modHt, file.getWeight());
                
                hZEffPt_den->Fill(genZPt, file.getWeight());
                hZEff_den->Fill(genZPt, modHt, file.getWeight());
                const bool& passMuZinvSel = tr.getVar<bool>("passMuZinvSel");
                if(passMuZinvSel)//genMatchMuInAcc.size() >= 2)
                {
                    hZEffPt_num->Fill(genZPt, file.getWeight());
                    hZEff_num->Fill(genZPt, modHt, file.getWeight());
                }
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

    hZEffPt_num->Write();
    hZEffPt_den->Write();
    hZAccPt_num->Write();
    hZAccPt_den->Write();
    hZEff_num->Write();
    hZEff_den->Write();
    hZAcc_num->Write();
    hZAcc_den->Write();

    hZEff_jActR1_num->Write();
    hZEff_jActR1_den->Write();
    hZEff_jActR2_num->Write();
    hZEff_jActR2_den->Write();

    f->Close();
}
