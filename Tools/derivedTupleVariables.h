#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/searchBins.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>

namespace plotterFunctions
{
    //Ugly global variables here
    static TH1* muEff;
    static TH2* muEff_jActR1;
    static TH1* muEffReco;
    static TH2* muEffIso;
    static TH2* muAcc;
    static TH1* hZEff;
    static TH1* hZAcc;
    static TRandom3 *tr3;
    static BaselineVessel *blvZinv;

//    topTagger::type3TopTagger * type3Ptr2;

    void generateWeight(NTupleReader& tr)
    {
        const std::vector<TLorentzVector>& jetsLVec         = tr.getVec<TLorentzVector>("jetsLVec");
        const std::vector<TLorentzVector>& cleanJetVec      = tr.getVec<TLorentzVector>("cleanJetVec");        
        const std::vector<TLorentzVector>& cutMuVec         = tr.getVec<TLorentzVector>("cutMuVec");
        const std::vector<double>& cutMuActivity            = tr.getVec<double>("cutMuActivity");
        const std::vector<TLorentzVector*>& genMu           = tr.getVec<TLorentzVector*>("genMu");
        const std::vector<TLorentzVector*>& genMuInAcc      = tr.getVec<TLorentzVector*>("genMuInAcc");
        const std::vector<TLorentzVector*>& genMatchMuInAcc = tr.getVec<TLorentzVector*>("genMatchMuInAcc");

        const int& pdgIdZDec      = tr.getVar<int>("pdgIdZDec");
        const bool& passMuZinvSel = tr.getVar<bool>("passMuZinvSel");
        const double& ht          = tr.getVar<double>("ht");
        const double& bestRecoZPt = tr.getVar<double>("bestRecoZPt");
        const double& genZPt      = tr.getVar<double>("genZPt");

        // Calculate PU weight

        // Calculate Z-eff weight

        const double zMassMin = 71.0;
        const double zMass    = 91.0;
        const double zMassMax = 111.0;

        double zMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
        for(int i = 0; i < cutMuVec.size(); ++i)
        {
            if(cutMuVec[i].Pt() < 10) continue;
            for(int j = 0; j < i && j < cutMuVec.size(); ++j)
            {
                if(cutMuVec[j].Pt() < 10) continue;
                double zm = (cutMuVec[i] + cutMuVec[j]).M();
                //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                {
                    zMassCurrent = zm;
                    if(muEff && muAcc)
                    {
                        double mu1pt = cutMuVec[i].Pt();
                        double mu2pt = cutMuVec[j].Pt();
                        double Ht = ht;

                        //set to not overflow histograms
                        if(mu1pt >= 2000.0) mu1pt = 1999.9;
                        if(mu2pt >= 2000.0) mu2pt = 1999.9;
                        if(Ht >= 3000.0) Ht = 2999.9;
                        
                        //Get mu efficiencies
                        double muEff1 = 0.0, muEff2 = 0.0;

                        if(muEff && muEffReco && muEffIso)
                        {
                            //Fit to reco eff (eff = p0 + p1*pt + p2*pt^2
                            const double fitStart = 200.0; // extended to 1400 GeV
                            const double p0 =     0.955847; // +/- 0.461944    
                            const double p1 = -2.24431e-05; // +/- 0.00128305  
                            const double p2 = -5.68907e-08; // +/- 7.85913e-07
                            double muRecoEff = 0.0;

                            int recoPtBin = muEffReco->GetXaxis()->FindBin(mu1pt);
                            if(recoPtBin >= muEffReco->FindBin(fitStart)) muRecoEff = p0 + p1*mu1pt + p2*mu1pt*mu1pt;
                            else                                          muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            int isoPtBin = muEffIso->GetXaxis()->FindBin(mu1pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            int isoActBin = muEffIso->GetYaxis()->FindBin(cutMuActivity[i]);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff1 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);

                            muRecoEff = 0.0;
                            recoPtBin = muEffReco->GetXaxis()->FindBin(mu2pt);
                            if(recoPtBin >= muEffReco->FindBin(fitStart)) muRecoEff = p0 + p1*mu2pt + p2*mu2pt*mu2pt;
                            else                                          muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            isoPtBin = muEffIso->GetXaxis()->FindBin(mu2pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            isoActBin = muEffIso->GetYaxis()->FindBin(cutMuActivity[j]);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff2 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);
                        }

                        //double tjActR1L1 = 0.0, tjActR1L2 = 0.0;//, jActR1wgm = 0.0, jActR2wgm = 0.0;
                        //for(auto& jet : cleanJetVec)
                        //{
                        //    double dR1 = ROOT::Math::VectorUtil::DeltaR(*jet, cutMuVec[i]);
                        //    double dR2 = ROOT::Math::VectorUtil::DeltaR(*jet, cutMuVec[j]);
                        //    if(dR1 > 0.04 && dR1 < 1.0)
                        //    {
                        //        tjActR1L1 += jet->Pt() / dR1;
                        //    }
                        //    if(dR2 > 0.04 && dR2 < 1.0)
                        //    {
                        //        tjActR1L2 += jet->Pt() / dR2;
                        //    }
                        //}
                        //
                        //if(muEff_jActR1)
                        //{
                        //    muEff1 = muEff_jActR1->GetBinContent(muEff_jActR1->GetXaxis()->FindBin(mu1pt), muEff_jActR1->GetYaxis()->FindBin(tjActR1L1));
                        //    muEff2 = muEff_jActR1->GetBinContent(muEff_jActR1->GetXaxis()->FindBin(mu2pt), muEff_jActR1->GetYaxis()->FindBin(tjActR1L2));
                        //}

                        
                        if((mu1pt > 20 && muEff1 < 1.0e-5) || (mu2pt > 20 && muEff2 < 1.0e-5)) 
                        {
                            std::cout << "SMALL muEff!!! muEff1: " << muEff1 << "\tmuEff2: " << muEff2 << "\t" << mu1pt << "\t" << mu2pt << std::endl;
                            zEff = 1.0e-10;
                        }
                        else zEff = muEff1 * muEff2;

                        //Get mu acceptance
                        /*double muAcc1 = 0.0, muAcc2 = 0.0;

                        muAcc1 = muAcc->GetBinContent(muAcc->GetXaxis()->FindBin(mu1pt), muAcc->GetYaxis()->FindBin(Ht));
                        muAcc2 = muAcc->GetBinContent(muAcc->GetXaxis()->FindBin(mu2pt), muAcc->GetYaxis()->FindBin(Ht));

                        if(muAcc1 < 1.0e-5 || muAcc2 < 1.0e-5) zAcc = 1.0e-10;
                        else                                   zAcc = muAcc1 * muAcc2;*/

                    }
                }
            }
        }

        double genCleanHt = ht;
        for(auto& tlvp : genMuInAcc) if(tlvp->Pt() > 50) genCleanHt -= tlvp->Pt();

        double mu1dRMin = 99.9, mu2dRMin = 99.9;
        for(auto& jet : cleanJetVec)
        {
            double mu1dR = 999.9, mu2dR = 999.9;
            if(cutMuVec.size() >= 1) mu1dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[0]);
            if(cutMuVec.size() >= 2) mu2dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[1]);
            mu1dRMin = std::min(mu1dRMin, mu1dR);
            mu2dRMin = std::min(mu2dRMin, mu2dR);
        }
        
        double mudR = 99.9, genMudR = 99.9, genMudPhi = 99.9, genMudEta = 99.9;

        if(cutMuVec.size() >= 2) mudR    = ROOT::Math::VectorUtil::DeltaR(cutMuVec[0], cutMuVec[1]);
        if(genMu.size() > 2) std::cout << "MORE THAN 2 GEN MUONS: " << genMu.size() << std::endl;
        if(genMu.size() >= 2)    
        {
            genMudR = ROOT::Math::VectorUtil::DeltaR(*(genMu[0]), *(genMu[1]));
            genMudPhi = ROOT::Math::VectorUtil::DeltaPhi(*(genMu[0]), *(genMu[1]));
            genMudEta = genMu[0]->Eta() - genMu[1]->Eta();
        }

        std::vector<double>* jActR1 = new std::vector<double>();
        std::vector<double>* jActR2 = new std::vector<double>();

        for(int i = 0; i < genMu.size(); ++i)
        {
            double tjActR1 = 0.0, tjActR2 = 0.0;//, jActR1wgm = 0.0, jActR2wgm = 0.0;
            //for(auto& jet : cleanJetVec)
            //{
            //    double dR = ROOT::Math::VectorUtil::DeltaR(*jet, *(genMu[i]));
            //    if(dR > 0.1 && dR < 1.0)
            //    {
            //        tjActR1 += jet->Pt() / dR;
            //        tjActR2 += jet->Pt() / (dR * dR);
            //    }
            //}
            for(auto& jet : jetsLVec)
            {
                double dR = ROOT::Math::VectorUtil::DeltaR(jet, *(genMu[i]));
                if(dR > 0.04 && dR < 1.0)
                {
                    tjActR1 += jet.Pt() / dR;
                    tjActR2 += jet.Pt() / (dR * dR);
                }
            }
            jActR1->push_back(tjActR1);
            jActR2->push_back(tjActR2);
        }
        
        //ZPt eff/Acc here
        //if(hZEff) zEff = hZEff->GetBinContent(hZEff->GetXaxis()->FindBin(bestRecoZ.Pt()));//, hZEff->GetYaxis()->FindBin(genCleanHt));//ht));
        //if(hZAcc) zAcc = hZAcc->GetBinContent(hZAcc->GetXaxis()->FindBin(bestRecoZPt));//, hZAcc->GetYaxis()->FindBin(genCleanHt));//ht));

        //functional form [2] - exp([0] + [1]*x)
        double acc_p0 = -2.91374e-01;
        double acc_p1 = -4.42884e-03;
        double acc_p2 =  9.51190e-01;
        if(hZAcc) 
        {
            if(genZPt < 100) zAcc = hZAcc->GetBinContent(hZAcc->GetXaxis()->FindBin(genZPt));
            else             zAcc = acc_p2 - exp(acc_p0 + acc_p1 * genZPt);
        }

        if(pdgIdZDec == 13 && passMuZinvSel && zAcc < 0.05) 
        {
            std::cout << "WARNING: Z acceptance < 0.05, forcing weight to zero! Acc_Z: " << zAcc << std::endl;
            zAcc = 1.0e101;
        }

        if(pdgIdZDec == 13 && passMuZinvSel && zEff < 0.04) 
        {
            std::cout << "WARNING: Z efficiency < 0.05, forcing weight to zero! Eff_Z: " << zEff << "\tZ(Pt): " << bestRecoZPt <<  std::endl;
            zEff = 1.0e101;
        }

        tr.registerDerivedVec("jActR1", jActR1);
        tr.registerDerivedVec("jActR2", jActR2);

        tr.registerDerivedVar("mu1dRMin", mu1dRMin);
        tr.registerDerivedVar("mu2dRMin", mu2dRMin);
        tr.registerDerivedVar("mudR", mudR);
        tr.registerDerivedVar("genMudR", genMudR);
        tr.registerDerivedVar("genMudPhi", genMudPhi);
        tr.registerDerivedVar("genMudEta", genMudEta);
        tr.registerDerivedVar("zEff", zEff);
        tr.registerDerivedVar("zEffWgt", 1.0/zEff);
        tr.registerDerivedVar("zAcc", zAcc);
        tr.registerDerivedVar("zAccWgt", 1.0/zAcc);
    }

    void muInfo(NTupleReader& tr)
    {
        const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
        const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
        const std::vector<TLorentzVector>& muonsLVec    = tr.getVec<TLorentzVector>("muonsLVec");
        const std::vector<double>& muonsRelIso          = tr.getVec<double>("muonsRelIso");
        const std::vector<double>& muonsMiniIso         = tr.getVec<double>("muonsMiniIso");
        const std::vector<int>& W_emuVec                = tr.getVec<int>("W_emuVec");
        const std::vector<double>& muonsCharge          = tr.getVec<double>("muonsCharge");
        const std::vector<TLorentzVector>& jetsLVec     = tr.getVec<TLorentzVector>("jetsLVec");
        const std::vector<double>& recoJetschargedEmEnergyFraction     = tr.getVec<double>("recoJetschargedEmEnergyFraction"); 
        const std::vector<double>& recoJetschargedHadronEnergyFraction = tr.getVec<double>("recoJetschargedHadronEnergyFraction");

        const double& ht                             = tr.getVar<double>("ht");
        const double& met                            = tr.getVar<double>("met");
        const double& metphi                         = tr.getVar<double>("metphi");

        std::vector<const TLorentzVector*>* genMatchIsoMuInAcc = new std::vector<const TLorentzVector*>();
        std::vector<const TLorentzVector*>* genMatchMuInAcc = new std::vector<const TLorentzVector*>();
        std::vector<double>* genMatchMuInAccRes = new std::vector<double>();
        std::vector<const TLorentzVector*>* genMuInAcc = new std::vector<const TLorentzVector*>();
        std::vector<const TLorentzVector*>* genMu = new std::vector<const TLorentzVector*>();
        std::vector<double>* genMatchIsoMuInAccAct = new std::vector<double>();
        std::vector<double>* genMatchMuInAccAct = new std::vector<double>();
        std::vector<double>* genMuInAccAct = new std::vector<double>();
        std::vector<double>* genMuAct = new std::vector<double>();
        std::vector<TLorentzVector>* cutMuVec = new std::vector<TLorentzVector>();
        std::vector<double>* cutMuCharge = new std::vector<double>();
        std::vector<double>* cutMuActivity = new std::vector<double>();

        std::vector<TLorentzVector> cutMuVecRecoOnly;

        int sumCharge = 0;
        int nTriggerMuons = 0;
        for(int i = 0; i < muonsLVec.size(); ++i)
        {
            if(AnaFunctions::passMuon( muonsLVec[i], 0.0, 0.0, AnaConsts::muonsArr)) // emulates muons with pt but no iso requirements.  
            {
                cutMuVecRecoOnly.push_back(muonsLVec[i]);
            }
            if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, AnaConsts::muonsArr))
            {
                if(AnaFunctions::passMuon( muonsLVec[i], muonsRelIso[i], 0.0, AnaConsts::muonsTrigArr)) 
                {
                    if(nTriggerMuons == 0 && muonsLVec[i].Pt() > 17)  nTriggerMuons++;
                    else if(muonsLVec[i].Pt() > 8)  nTriggerMuons++;
                }
                cutMuVec->push_back(muonsLVec[i]);
                cutMuCharge->push_back(muonsCharge[i]);
                cutMuActivity->push_back(AnaFunctions::getMuonActivity(muonsLVec[i], jetsLVec, recoJetschargedHadronEnergyFraction, recoJetschargedEmEnergyFraction, AnaConsts::muonsAct));
                if(muonsCharge[i] > 0) sumCharge++;
                else                   sumCharge--;
            }
        }

        //mu45 non-iso trigger emulation 
        const double effsnom2012ABC[] = {0.928,0.8302,0.8018};
        const double upedge2012ABC[] = { 0.9, 1.2, 2.1};
        bool muTrigMu45 = false;
        for(TLorentzVector& mu : *cutMuVec)
        {
            if(mu.Pt() > 45)
            {
                for(int iBin = 0; iBin < sizeof(effsnom2012ABC)/sizeof(double); ++iBin)
                {
                    if(mu.Eta() < upedge2012ABC[iBin] && tr3 && tr3->Uniform() < effsnom2012ABC[iBin])
                    {
                        muTrigMu45 = true;
                        break;
                    }
                }
                if(muTrigMu45) break;
            }
        }
        
        double genHt = 0.0;

        const double minMuPt = 20.0, highMuPt = 45.0;
        double nuPt1 = -999.9, nuPt2 = -999.9;

        for(int i = 0; i < genDecayPdgIdVec.size() && i < genDecayLVec.size(); ++i)
        {
            if((abs(genDecayPdgIdVec[i]) != 0 &&  abs(genDecayPdgIdVec[i]) < 6) || (abs(genDecayPdgIdVec[i]) > 100 && abs(genDecayPdgIdVec[i]) < 10000)) genHt += genDecayLVec[i].Pt();

            if(genDecayPdgIdVec[i] ==  13) nuPt1 = genDecayLVec[i].Pt();
            if(genDecayPdgIdVec[i] == -13) nuPt2 = genDecayLVec[i].Pt();
                
            if(abs(genDecayPdgIdVec[i]) == 13)
            {
                genMu->push_back(&genDecayLVec[i]);
                genMuAct->push_back(AnaFunctions::getMuonActivity(genDecayLVec[i], jetsLVec, recoJetschargedHadronEnergyFraction, recoJetschargedEmEnergyFraction, AnaConsts::muonsAct));
                if(AnaFunctions::passMuonAccOnly(genDecayLVec[i], AnaConsts::muonsMiniIsoArr) && genDecayLVec[i].Pt() > minMuPt)
                {
                    genMuInAcc->push_back(&genDecayLVec[i]);
                    genMuInAccAct->push_back(genMuAct->back());
                    double dRMin = 999.9;
                    double matchPt = -999.9;
                    for(int j = 0; j < cutMuVecRecoOnly.size(); ++j)
                    {
                        double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutMuVecRecoOnly[j]);
                        if(dR < dRMin)
                        {
                            dRMin = dR;
                            matchPt = cutMuVecRecoOnly[j].Pt();
                        }
                    }
                    if(dRMin < 0.02)
                    {
                        genMatchMuInAcc->push_back(&genDecayLVec[i]);
                        genMatchMuInAccAct->push_back(genMuAct->back());
                        genMatchMuInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                    }
                    
                    dRMin = 999.9;
                    for(int j = 0; j < cutMuVec->size(); ++j)
                    {
                        double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], (*cutMuVec)[j]);
                        if(dR < dRMin)
                        {
                            dRMin = dR;
                        }
                    }
                    if(dRMin < 0.02)
                    {
                        genMatchIsoMuInAcc->push_back(&genDecayLVec[i]);
                        genMatchIsoMuInAccAct->push_back(genMuAct->back());
                    }
                }
            }
        }

        double genZPt = -999.9, genZEta = -999.9, genZmass = -999.9, genZPhi;
        int nZ = 0;
        TLorentzVector genZ;
        for(int j = 0; j <  genDecayPdgIdVec.size(); ++j)
        {
            if(abs(genDecayPdgIdVec[j]) == 23)
            {
                nZ++;
                genZ = genDecayLVec[j];
                genZPt = genDecayLVec[j].Pt();
                genZEta = genDecayLVec[j].Eta();
                genZPhi = genDecayLVec[j].Phi();
                genZmass = genDecayLVec[j].M();
            }
        }
        if(nZ > 1) std::cout << "!!!WARNING MORE THAN 1 Z FOUND!!!" << std::endl;


        int pdgIdZDec = 0;
        if(W_emuVec.size() == 0) pdgIdZDec = 15;
        else if(W_emuVec.size() == 2)
        {
            if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 11) pdgIdZDec = 11;
            else if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 13) pdgIdZDec = 13;
        }

        bool passDiMuTrig  = nTriggerMuons >= 2;
        
        const double zMassMin = 71.0;
        const double zMass    = 91.0;
        const double zMassMax = 111.0;

        double zMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
        TLorentzVector bestRecoZ;
        for(int i = 0; i < cutMuVec->size(); ++i)
        {
            if((*cutMuVec)[i].Pt() < minMuPt) continue;
            for(int j = 0; j < i && j < cutMuVec->size(); ++j)
            {
                if((*cutMuVec)[j].Pt() < minMuPt) continue;
                double zm = ((*cutMuVec)[i] + (*cutMuVec)[j]).M();
                //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                {
                    bestRecoZ = (*cutMuVec)[i] + (*cutMuVec)[j];
                    zMassCurrent = zm;
                }
            }
        }

        //if(genMuInAcc->size() >= 2)
        //{
        //    double mmm = (*genMuInAcc->at(0) + *genMuInAcc->at(1)).M();
        //    if(fabs(genZmass - mmm)/genZmass > 0.2 && mmm > 71 && mmm < 111) std::cout << "genZmass: " << genZmass << "\tgen dimuon mass: " << mmm << "\treco dimuon mass: " << zMassCurrent << "\tgenZPt: " << genZPt << "\tmu1pt: " << (*genMuInAcc->at(0)).Pt() << "\tmu2pt: " << (*genMuInAcc->at(1)).Pt() << std::endl;
        //}



        TLorentzVector metV, metZ;
        metV.SetPtEtaPhiM(met, 0.0, metphi, 0.0);
        metZ.SetPtEtaPhiM(bestRecoZ.Pt(), 0.0, bestRecoZ.Phi(), 0.0);
        TLorentzVector cleanMet = metV + metZ;

        bool passMuZinvSel = (cutMuVec->size() == 2 && sumCharge == 0 && (*cutMuVec)[0].Pt() > highMuPt && (*cutMuVec)[1].Pt() > minMuPt) && (bestRecoZ.M() > zMassMin) && (bestRecoZ.M() < zMassMax);        

        double cutMuPt1 = -999.9;
        double cutMuPt2 = -999.9;
        if(cutMuVec->size() >= 1) cutMuPt1 = cutMuVec->at(0).Pt();
        if(cutMuVec->size() >= 2) cutMuPt2 = cutMuVec->at(1).Pt();

        const unsigned int& run   = tr.getVar<unsigned int>("run");
        const unsigned int& lumi  = tr.getVar<unsigned int>("lumi");
        const unsigned int& event = tr.getVar<unsigned int>("event");
        //if(passMuZinvSel && genZPt > 400 && fabs(nuPt1 - nuPt2) > 300) std::cout << "BONJOUR!!! \t" << genZPt << "\t" << fabs(nuPt1 - nuPt2) << "\t" << run << "\t" << lumi << "\t" << event << std::endl;

        double mindPhiMetJ = 999.9;
        int jc = 0;
        for(const TLorentzVector& jet : jetsLVec)
        {
            if(jc >= 3) break;
            jc++;
            mindPhiMetJ = std::min(mindPhiMetJ, fabs(ROOT::Math::VectorUtil::DeltaPhi(genZ, jet)));
        }
        //if(genZPt > 600) std::cout << "HELLO THERE!!!!" << std::endl;
        //if(genZPt > 600 && mindPhiMetJ < 0.5) std::cout << "BONJOUR!!! \t" << genZPt << "\t" << mindPhiMetJ << "\t" << run << "\t" << lumi << "\t" << event << std::endl;

        double bestRecoZPt = bestRecoZ.Pt();
        double cleanMetPt = cleanMet.Pt();
        //double cleanMet2Pt = cleanMet2.Pt();
        tr.registerDerivedVar("bestRecoZPt", bestRecoZPt);
        tr.registerDerivedVar("bestRecoZM", bestRecoZ.M());
        tr.registerDerivedVar("cleanMetPt", cleanMetPt);
        tr.registerDerivedVar("cleanMetPhi", cleanMet.Phi());
        //tr.registerDerivedVar("cleanMet2Pt", cleanMet2Pt);
        tr.registerDerivedVar("genHt", genHt);
        tr.registerDerivedVar("cutMuPt1", cutMuPt1);
        tr.registerDerivedVar("cutMuPt2", cutMuPt2);
        tr.registerDerivedVar("mindPhiMetJ", mindPhiMetJ);

        tr.registerDerivedVar("ZPtRes", (bestRecoZPt - genZPt)/genZPt);
        tr.registerDerivedVar("ZEtaRes", bestRecoZ.Eta() - genZEta);
        tr.registerDerivedVar("ZPhiRes", bestRecoZ.Phi() - genZPhi);
        tr.registerDerivedVar("ZMRes", (bestRecoZ.M() - genZmass)/genZmass);

        tr.registerDerivedVec("cutMuVec", cutMuVec);
        tr.registerDerivedVec("cutMuActivity", cutMuActivity);
        tr.registerDerivedVec("genMu", genMu);
        tr.registerDerivedVar("ngenMu", static_cast<double>(genMu->size()));
        tr.registerDerivedVec("genMuInAcc", genMuInAcc);
        tr.registerDerivedVec("genMuAct", genMuAct);
        tr.registerDerivedVar("ngenMuInAcc", static_cast<double>(genMuInAcc->size()));
        tr.registerDerivedVec("genMuInAccAct", genMuInAccAct);
        tr.registerDerivedVec("genMatchMuInAcc", genMatchMuInAcc);
        tr.registerDerivedVec("genMatchMuInAccRes", genMatchMuInAccRes);
        tr.registerDerivedVec("genMatchIsoMuInAcc", genMatchIsoMuInAcc);
        tr.registerDerivedVar("ngenMatchMuInAcc", static_cast<double>(genMatchMuInAcc->size()));
        tr.registerDerivedVec("genMatchMuInAccAct", genMatchMuInAccAct);
        tr.registerDerivedVec("genMatchIsoMuInAccAct", genMatchIsoMuInAccAct);
        tr.registerDerivedVar("genZPt", genZPt);
        tr.registerDerivedVar("genZEta", genZEta);
        tr.registerDerivedVar("genZmass", genZmass);
        tr.registerDerivedVar("pdgIdZDec", pdgIdZDec);
        tr.registerDerivedVar("passMuZinvSel", passMuZinvSel);
        tr.registerDerivedVar("passDiMuIsoTrig", passDiMuTrig);
        tr.registerDerivedVar("passSingleMu45", muTrigMu45);
    }

//    void cleanJets(NTupleReader& tr)
//    {
//        const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVec");
//        const std::vector<TLorentzVector>& elesLVec  = tr.getVec<TLorentzVector>("elesLVec");
//        const std::vector<TLorentzVector>& muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
//        const std::vector<double>& elesRelIso        = tr.getVec<double>("elesRelIso");
//        const std::vector<double>& muonsRelIso       = tr.getVec<double>("muonsRelIso");
//        const std::vector<double>& muonsMiniIso      = tr.getVec<double>("muonsMiniIso");
//        const std::vector<double>& recoJetsBtag_0    = tr.getVec<double>("recoJetsBtag_0");
//        const std::vector<int>& muMatchedJetIdx       = tr.getVec<int>("muMatchedJetIdx");
//        const std::vector<int>& eleMatchedJetIdx      = tr.getVec<int>("eleMatchedJetIdx");
//
//        if(elesLVec.size() != elesRelIso.size() || muonsLVec.size() != muonsRelIso.size())
//        {
//            std::cout << "MISMATCH IN VECTOR SIZE!!!!!" << std::endl;
//            return;
//        }
//
//        std::vector<const TLorentzVector*>* cleanJetVec = new std::vector<const TLorentzVector*>();
//        std::vector<TLorentzVector>* zinvJetVec = new std::vector<TLorentzVector>();
//        std::vector<double>* zinvBTag = new std::vector<double>;
//        std::vector<double>* mindR = new std::vector<double>;
//
//        const double jldRMax = 0.15;
//
//        const double HT_jetPtMin = 50;
//        const double HT_jetEtaMax = 2.4;
//        const double MTH_jetPtMin = 30.0;
//
//        double HT = 0.0, HTNoIso = 0.0;
//        TLorentzVector MHT;
//
//        std::vector<bool> keepJet(jetsLVec.size(), true);
//        std::vector<bool> keepJetPFCandMatch(jetsLVec.size(), true);
//
//        if(muonsLVec.size() != muonsMiniIso.size() || muonsLVec.size() != muMatchedJetIdx.size()) std::cout << "Electron vector size missmatch" << std::endl;
//
//        for(int iM = 0; iM < muonsLVec.size() && iM < muonsMiniIso.size() && iM < muMatchedJetIdx.size(); ++iM)
//        {
//            double dRmin = 999.0;
//            int minJMatch = -1;
//
//            if(!AnaFunctions::passMuon(muonsLVec[iM], muonsMiniIso[iM], 0.0, AnaConsts::muonsMiniIsoArr)) continue;
//
//            for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
//            {
//                if(!keepJet[iJet]) continue;
//
//                double dR = ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], muonsLVec[iM]);
//                if(dR < dRmin)
//                {
//                    dRmin = dR;
//                    minJMatch = iJet;
//                }
//            }
//
//            mindR->push_back(dRmin);
//            if(minJMatch >= 0 && dRmin < jldRMax) keepJet[minJMatch] = false;
//
//            //if(muMatchedJetIdx[iM] != minJMatch) std::cout << "Different muon match between dR and PFCand matching found!!!" << std::endl;
//
//            if(muMatchedJetIdx[iM] >= 0) keepJetPFCandMatch[muMatchedJetIdx[iM]] = false;
//            else                         keepJetPFCandMatch[minJMatch] = false;
//
//        }
//
//        if(elesLVec.size() != elesRelIso.size() || elesLVec.size() != eleMatchedJetIdx.size()) std::cout << "Electron vector size mismatch\t" << elesLVec.size() << "\t" << elesRelIso.size() << "\t" << eleMatchedJetIdx.size() << std::endl;
//
//        for(int iE = 0; iE < elesLVec.size() && iE < elesRelIso.size() && iE < eleMatchedJetIdx.size(); ++iE)
//        {
//            double dRmin = 999.0;
//            int minJMatch = -1;
//
//            if(!AnaFunctions::passElectron(elesLVec[iE], elesRelIso[iE], 0.0, AnaConsts::elesArr)) continue;
//
//            for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
//            {
//                if(!keepJet[iJet]) continue;
//
//                double dR = ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], elesLVec[iE]);
//                if(dR < dRmin)
//                {
//                    dRmin = dR;
//                    minJMatch = iJet;
//                }
//            }
//            
//            mindR->push_back(dRmin);
//            if(minJMatch >= 0 && dRmin < jldRMax) keepJet[minJMatch] = false;
//
//            //if(eleMatchedJetIdx[iE] != minJMatch) std::cout << "Different electron match between dR and PFCand matching found!!!" << std::endl;
//
//            if(eleMatchedJetIdx[iE] >= 0) keepJetPFCandMatch[eleMatchedJetIdx[iE]] = false;
//            else                          keepJetPFCandMatch[minJMatch] = false;
//        }
//
//        int jetsKept = 0;
//        for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
//        {
//            if(keepJetPFCandMatch[iJet])
//            {
//                ++jetsKept;
//                cleanJetVec->push_back(&jetsLVec[iJet]);
//                if(AnaFunctions::jetPassCuts(jetsLVec[iJet], AnaConsts::pt30Arr))
//                {
//                    zinvJetVec->push_back(jetsLVec[iJet]);
//                    zinvBTag->push_back(recoJetsBtag_0[iJet]);
//                }
//                if(jetsLVec[iJet].Pt() > HT_jetPtMin && fabs(jetsLVec[iJet].Eta()) < HT_jetEtaMax) HT += jetsLVec[iJet].Pt();
//                if(jetsLVec[iJet].Pt() > MTH_jetPtMin) MHT += jetsLVec[iJet];
//            }
//        }
//
//        tr.registerDerivedVar("nJetsRemoved", static_cast<int>(jetsLVec.size() - jetsKept));
//        tr.registerDerivedVar("cleanHt", HT);
//        tr.registerDerivedVar("cleanMHt", MHT.Pt());
//        tr.registerDerivedVar("cleanMHtPhi", MHT.Phi());
//        tr.registerDerivedVec("cleanJetVec", cleanJetVec);
//        tr.registerDerivedVec("zinvJetVec", zinvJetVec);
//        tr.registerDerivedVec("zinvBTagVec", zinvBTag);
//        tr.registerDerivedVec("minljdR", mindR);
//    }

//    void zinvBaseline(NTupleReader& tr)
//    {
//        const bool& passMuZinvSel = tr.getVar<bool>("passMuZinvSel");
//
//        const double& cleanMetPt  = tr.getVar<double>("cleanMetPt");
//        const double& cleanMetPhi = tr.getVar<double>("cleanMetPhi");
//
//        const std::vector<TLorentzVector>& zinvJetVec  = tr.getVec<TLorentzVector>("cleanJetpt30ArrVec");
//        const std::vector<double>&         zinvBTagVec = tr.getVec<double>("cleanJetpt30ArrBTag");
//
//        TLorentzVector metLVec;
//        metLVec.SetPtEtaPhiM(tr.getVar<double>("cleanMetPt"), 0.0, cleanMetPhi, 0.0);
//
//        // Calculate number of jets and b-tagged jets
//        int cntCSVS = AnaFunctions::countCSVS(zinvJetVec, zinvBTagVec, AnaConsts::cutCSVS, AnaConsts::bTagArr);
//        int cntNJetsPt50Eta24 = AnaFunctions::countJets(zinvJetVec, AnaConsts::pt50Eta24Arr);
//        int cntNJetsPt30Eta24 = AnaFunctions::countJets(zinvJetVec, AnaConsts::pt30Eta24Arr);
//        int cntNJetsPt30      = AnaFunctions::countJets(zinvJetVec, AnaConsts::pt30Arr);
//
//        // Recalculate deltaPhi
//        std::vector<double> *dPhiVec = new std::vector<double>();
//        *dPhiVec = AnaFunctions::calcDPhi(zinvJetVec, cleanMetPhi, 3, AnaConsts::dphiArr);
//
//        // Calculate number of leptons
//        int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
//        int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesArr);
//        int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);
//
//
//        // Pass cuts
//
//        bool passZinvdPhis = true;
//        if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ) passZinvdPhis = false;
//
//        bool passZinvnJets = true;
//        if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ) passZinvnJets = false;
//        if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ) passZinvnJets = false;
//
//        bool passZinvBJets = true;
//        if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ) passZinvBJets = false;
//
//        bool passZinvLeptVeto = true;
//        if( nElectrons != AnaConsts::nElectronsSel ) passZinvLeptVeto = false;
//        if( nIsoTrks != AnaConsts::nIsoTrksSel )     passZinvLeptVeto = false;
//        
//        bool passStandLeptVeto = passZinvLeptVeto;
//        if( nMuons != AnaConsts::nMuonsSel ) passStandLeptVeto = false;
//
//        bool passZinvMET = true;
//        if( tr.getVar<double>("cleanMetPt") < AnaConsts::defaultMETcut ) passZinvMET = false;
//
//        // Calculate top tagger related variables. 
//        // Note that to save speed, only do the calculation after previous base line requirements.
//        int bestTopJetIdx = -1;
//        bool remainPassCSVS = false;
//        int pickedRemainingCombfatJetIdx = -1;
//        double bestTopJetMass = -1;
//        int nTopCandSortedCnt = 0;
//        double MT2 = -1;
//        double mTcomb = -1;
//
//        bool passPreTTag = passZinvLeptVeto && passZinvnJets && passZinvdPhis/* && passZinvBJets */&& passZinvMET && cntNJetsPt30 >= AnaConsts::nJetsSel;
//        bool passZinvTagger = false;
//        if(type3Ptr && type3Ptr2 && passPreTTag)
//        {
//            passZinvTagger = true; //will be set fals below if it fails 
//            topTagger::type3TopTagger * t3Ptr = nullptr;
//            if(passZinvBJets) t3Ptr = type3Ptr;
//            else              t3Ptr = type3Ptr2;
//            t3Ptr->processEvent(zinvJetVec, zinvBTagVec, metLVec);
//            bestTopJetIdx = t3Ptr->bestTopJetIdx;
//            remainPassCSVS = t3Ptr->remainPassCSVS;
//            pickedRemainingCombfatJetIdx = t3Ptr->pickedRemainingCombfatJetIdx;
//            if( bestTopJetIdx != -1 ) bestTopJetMass = t3Ptr->bestTopJetLVec.M();
//
//            nTopCandSortedCnt = t3Ptr->nTopCandSortedCnt;
//            MT2 = t3Ptr->MT2;
//            mTcomb = t3Ptr->mTbJet + 0.5*t3Ptr->mTbestTopJet;
//        }
//
//        if( bestTopJetIdx == -1 )                                                                    passZinvTagger = false;
//        if( ! remainPassCSVS )                                                                       passZinvTagger = false;
//        if( pickedRemainingCombfatJetIdx == -1 && zinvJetVec.size()>=6 )                             passZinvTagger = false;
//        if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ ) ) passZinvTagger = false;
//
//        bool passZinvBaseline      = passZinvLeptVeto && passZinvnJets && passZinvdPhis && passZinvBJets && passZinvMET && passZinvTagger;
//        bool passZinvBaselineNoTag = passZinvLeptVeto && passZinvnJets && passZinvdPhis                  && passZinvMET;
//
//        tr.registerDerivedVec("dPhiZinvVec", dPhiVec);
//
//        tr.registerDerivedVar("bestTopJetIdxZinv", bestTopJetIdx);
//        tr.registerDerivedVar("remainPassCSVSZinv", remainPassCSVS);
//        tr.registerDerivedVar("pickedRemainingCombfatJetIdxZinv", pickedRemainingCombfatJetIdx);
//        tr.registerDerivedVar("bestTopJetMassZinv", bestTopJetMass);
//
//        tr.registerDerivedVar("passZinvMET", passZinvMET);
//        tr.registerDerivedVar("passZinvLeptVeto", passZinvLeptVeto);
//        tr.registerDerivedVar("passStandLeptVeto", passStandLeptVeto);
//        tr.registerDerivedVar("passZinvnJets", passZinvnJets);
//        tr.registerDerivedVar("passZinvdPhis", passZinvdPhis);
//        tr.registerDerivedVar("passZinvBJets", passZinvBJets);
//        tr.registerDerivedVar("passZinvTagger", passZinvTagger);
//        tr.registerDerivedVar("passZinvBaseline", passZinvBaseline);
//        tr.registerDerivedVar("passZinvBaselineNoTag", passZinvBaselineNoTag);
//
//        tr.registerDerivedVar("cntNJetsPt30Eta24Zinv", cntNJetsPt30Eta24);
//        tr.registerDerivedVar("cntCSVSZinv", cntCSVS);
//
//        tr.registerDerivedVar("nTopCandSortedCntZinv", nTopCandSortedCnt);
//        tr.registerDerivedVar("MT2Zinv", MT2);
//        tr.registerDerivedVar("mTcombZinv", mTcomb);
//    }

    void getSearchBin(NTupleReader& tr)
    {
        const int& cntCSVS = tr.getVar<int>("cntCSVSZinv");
        const int& nTopCandSortedCnt = tr.getVar<int>("nTopCandSortedCntZinv");
        const double& cleanMet = tr.getVar<double>("cleanMetPt");
        const double& cleanMetPhi = tr.getVar<double>("cleanMetPhi");
        const double& MT2 = tr.getVar<double>("best_had_brJet_MT2Zinv");
        const std::vector<TLorentzVector>& removedJetsLVec = tr.getVec<TLorentzVector>("removedJetVec");

        int nSearchBin = find_Binning_Index(cntCSVS, nTopCandSortedCnt, MT2, cleanMet);

        //hack, this does not belong here
        TLorentzVector cleanMet2;
        cleanMet2.SetPtEtaPhiM(cleanMet, 0.0, cleanMetPhi, 0.0);
        
        for(auto& jet : removedJetsLVec)
        {
            cleanMet2 += jet;
        }

        std::vector<std::pair<double, double> > * nb0Bins = new std::vector<std::pair<double, double> >();
        const double wnb01 = 2.2322e-01;
        const double wnb02 = 4.3482e-02;
        const double wnb03 = 4.5729e-03;
        if(cntCSVS == 0)
        {
            //nb0Bins->push_back(std::make_pair(find_Binning_Index(0, nTopCandSortedCnt, MT2, cleanMet), 1.0));
            nb0Bins->emplace_back(std::pair<double, double>(find_Binning_Index(1, nTopCandSortedCnt, MT2, cleanMet), wnb01));
            nb0Bins->emplace_back(std::pair<double, double>(find_Binning_Index(2, nTopCandSortedCnt, MT2, cleanMet), wnb02));
            nb0Bins->emplace_back(std::pair<double, double>(find_Binning_Index(3, nTopCandSortedCnt, MT2, cleanMet), wnb03));
        }

        tr.registerDerivedVar("nSearchBin", nSearchBin);
        tr.registerDerivedVar("cleanMet2Pt", double(cleanMet2.Pt()));
        tr.registerDerivedVec("nb0Bins", nb0Bins);
    }

    void printInterestingEvents(NTupleReader& tr)
    {
        const unsigned int& run   = tr.getVar<unsigned int>("run");
        const unsigned int& event = tr.getVar<unsigned int>("event");

        const double& met                            = tr.getVar<double>("met");
        const double& metphi                         = tr.getVar<double>("metphi");

        const int& nMuons_CUT        = tr.getVar<int>("nMuons_CUT");
        const int& nElectrons_CUT    = tr.getVar<int>("nElectrons_CUT");
        const int& cntNJetsPt50Eta24 = tr.getVar<int>("cntNJetsPt50Eta24");
        const int& cntNJetsPt30Eta24 = tr.getVar<int>("cntNJetsPt30Eta24");
        const int& cntNJetsPt30      = tr.getVar<int>("cntNJetsPt30");

        const double& mht    = tr.getVar<double>("mht");
        const double& mhtphi = tr.getVar<double>("mhtphi");
        const double& ht     = tr.getVar<double>("ht");

        //if(met > 1000) std::cout << "run: " << run << "\tevent: " << event << "\tmet: " << met << "\tmetphi: " << metphi << "\tnMuons_CUT: " << nMuons_CUT << "\t nElectrons_CUT: " << nElectrons_CUT << "\tcntNJetsPt30: " << cntNJetsPt30 << "\tcntNJetsPt30Eta24: " << cntNJetsPt30Eta24 << "\tcntNJetsPt50Eta24: " << cntNJetsPt50Eta24 << "\tmht: " << mht << "\tmhtphi: " << mhtphi << "\tht: " << ht << std::endl;
    }
    
    void zinvBaseline(NTupleReader& tr)
    {
        (*blvZinv)(tr);
    }

    void registerFunctions(NTupleReader& tr)
    {
        //Make some global "constants" here
        TH1::AddDirectory(false);
        TFile *f = new TFile("muEffHists.root");
        if(f)
        {
            muEff = static_cast<TH1*>(f->Get("hMuEffPt"));
            muEffReco = static_cast<TH1*>(f->Get("hMuEffPtReco"));
            muEffIso  = static_cast<TH2*>(f->Get("hMuEffPtActIso"));
            muEff_jActR1 = static_cast<TH2*>(f->Get("hZEff_jActR1"));
            muAcc = static_cast<TH2*>(f->Get("hMuAcc"));
            hZEff = static_cast<TH1*>(f->Get("hZEffPt"));
            hZAcc = static_cast<TH1*>(f->Get("hZAccPtSmear"));
            f->Close();
            delete f;
        }

//        type3Ptr2 = new topTagger::type3TopTagger();
//        type3Ptr2->setnJetsSel(AnaConsts::nJetsSel);
//        type3Ptr2->setdobVetoCS(true);
//
        tr3 = new TRandom3();

        blvZinv = new BaselineVessel("Zinv");

        //register functions with NTupleReader
        tr.registerFunction(&muInfo);
        stopFunctions::cjh.setMuonIso("mini");
        stopFunctions::cjh.setRemove(false);
        stopFunctions::cjh.setDisable(false);
        tr.registerFunction(&stopFunctions::cleanJets);
        tr.registerFunction(&generateWeight);
        tr.registerFunction(&zinvBaseline);
        tr.registerFunction(&getSearchBin);
        //tr.registerFunction(&printInterestingEvents);
    }

    void activateBranches(std::set<std::string>& activeBranches)
    {
        for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
        activeBranches.insert("ht");
        activeBranches.insert("run");
        activeBranches.insert("lumi");
        activeBranches.insert("event");
        activeBranches.insert("mht");
        activeBranches.insert("mhtphi");
        activeBranches.insert("genDecayPdgIdVec");
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
        activeBranches.insert("loose_isoTrks_mtw");
        activeBranches.insert("elesMtw");
        activeBranches.insert("loose_isoTrks_iso");
        activeBranches.insert("loose_isoTrksLVec");
        activeBranches.insert("muMatchedJetIdx");
        activeBranches.insert("eleMatchedJetIdx");
        activeBranches.insert("recoJetschargedEmEnergyFraction"); 
        activeBranches.insert("recoJetschargedHadronEnergyFraction");
        activeBranches.insert("elesisEB");
    }
}
