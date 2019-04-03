#ifndef LEPINFO_H 
#define LEPINFO_H 

#include "TypeDefinitions.h"
#include "PhotonTools.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"

#include "TopTagger.h"
#include "TTModule.h"
#include "TopTaggerUtilities.h"
#include "TopTaggerResults.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"

#include "TopTagger/TopTagger/interface/TopObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <set>

namespace plotterFunctions
{
    class LepInfo
    {
    private:
        // use shared_ptr, which will delete and clear dynamically allocated memory
        std::shared_ptr<TRandom3> tr3;
        void lepInfo(NTupleReader& tr)
        {
            std::vector<int> genDecayPdgIdVec;
            std::vector<int> GenPart_statusFlags;
            std::vector<TLorentzVector> genDecayLVec;
            if (tr.checkBranch("GenPartTLV"))
            {
                genDecayPdgIdVec                    = tr.getVec<int>("GenPart_pdgId");
                GenPart_statusFlags                 = tr.getVec<int>("GenPart_statusFlags");
                genDecayLVec                        = tr.getVec<TLorentzVector>("GenPartTLV");
            }
            const auto& muonsLVec                           = tr.getVec<TLorentzVector>("MuonTLV");
            const auto& muonsCharge                         = tr.getVec<int>("Muon_charge");
            const auto& jetsLVec                            = tr.getVec<TLorentzVector>("JetTLV");
            const auto& cutMuVec                            = tr.getVec<TLorentzVector>("cutMuVec"); 
            const auto& cutMuVecRecoOnly                    = tr.getVec<TLorentzVector>("cutMuVecRecoOnly"); 
            const auto& cutMuSummedCharge                   = tr.getVar<int>("cutMuSummedCharge"); 
            const auto& nTriggerMuons                       = tr.getVar<int>("nTriggerMuons"); 
            const auto& cutElecVec                          = tr.getVec<TLorentzVector>("cutElecVec"); 
            const auto& cutElecVecRecoOnly                  = tr.getVec<TLorentzVector>("cutElecVecRecoOnly"); 
            const auto& cutElecSummedCharge                 = tr.getVar<int>("cutElecSummedCharge"); 
            const auto& met                                 = tr.getVar<data_t>("MET_pt");
            const auto& metphi                              = tr.getVar<data_t>("MET_phi");

            bool passMuonVeto = false;
            bool passEleVeto = false;
            
            try
            {
                const auto& passMuonVetoTmp  = tr.getVar<bool>("passMuonVeto");
                const auto& passEleVetoTmp   = tr.getVar<bool>("passEleVeto");
                if(&passMuonVetoTmp != nullptr) passMuonVeto = passMuonVetoTmp;
                if(&passEleVetoTmp != nullptr) passEleVeto = passEleVetoTmp;
            }
            catch(const std::string e)
            {
                //std::cout << "void muInfo(NTupleReader& tr): Caught exception, variable \"" << e << "\" not found" << std::endl;
            }

            auto* genMatchIsoElecInAcc      = new std::vector<TLorentzVector>();
            auto* genMatchElecInAcc         = new std::vector<TLorentzVector>();
            auto* genMatchElecInAccRes      = new std::vector<data_t>();
            auto* genElecInAcc              = new std::vector<TLorentzVector>();
            auto* genElec                   = new std::vector<TLorentzVector>();
            auto* genMatchIsoMuInAcc        = new std::vector<TLorentzVector>();
            auto* genMatchMuInAcc           = new std::vector<const TLorentzVector*>();
            auto* genMatchMuInAccRes        = new std::vector<data_t>();
            auto* genMuInAcc                = new std::vector<TLorentzVector>();
            auto* genMu                     = new std::vector<const TLorentzVector*>();

            //mu45 non-iso trigger emulation
            const double effsnom2012ABC[] = {0.928,0.8302,0.8018};
            const double upedge2012ABC[] = { 0.9, 1.2, 2.1};
            bool muTrigMu45 = false;
            for(const TLorentzVector& mu : cutMuVec)
            {
                if(mu.Pt() > 50)
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

            // muon trigger turn on:     50 GeV
            // electron trigger turn on: 50 GeV
            const double   minMuPt = 20.0,   highMuPt = 50.0;
            const double minElecPt = 20.0, highElecPt = 50.0;

            //Gen info parsing
            const int GENPARTMASK = 0x2100;
            
            if(tr.checkBranch("GenPartTLV") && &genDecayLVec != nullptr)
            {
                for (int i = 0; i < genDecayPdgIdVec.size(); ++i)
                {
                    //int i = W_emuVec[index];
                    int maskedStatusFlag = (GenPart_statusFlags[i] & GENPARTMASK);
                    bool isGoodGenPart = (maskedStatusFlag == GENPARTMASK);
                    if (isGoodGenPart)
                    {
                        //muon efficiency and acceptance
                        if(abs(genDecayPdgIdVec[i]) == 13)
                        {
                            genMu->push_back(&genDecayLVec[i]);
                            if(AnaFunctions::passMuonAccOnly(genDecayLVec[i], AnaConsts::muonsMiniIsoArr) && genDecayLVec[i].Pt() > minMuPt)
                            {
                                genMuInAcc->push_back(genDecayLVec[i]);
                                double dRMin = 999.9;
                                double matchPt = -999.9;
                                for(int j = 0; j < cutMuVecRecoOnly.size(); ++j)
                                {
                                    // difference in angle between gen and reco muons 
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
                                    genMatchMuInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                                }

                                dRMin = 999.9;
                                for(int j = 0; j < cutMuVec.size(); ++j)
                                {
                                    // difference in angle between gen and cut muons 
                                    double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutMuVec[j]);
                                    if(dR < dRMin)
                                    {
                                        dRMin = dR;
                                    }
                                }
                                if(dRMin < 0.02)
                                {
                                    genMatchIsoMuInAcc->push_back(genDecayLVec[i]);
                                }
                            }
                        }

                        //electron efficiency and acceptance
                        if(abs(genDecayPdgIdVec[i]) == 11)
                        {
                            genElec->push_back(genDecayLVec[i]);
                            if(AnaFunctions::passElectronAccOnly(genDecayLVec[i], AnaConsts::elesMiniIsoArr) && genDecayLVec[i].Pt() > minElecPt)
                            {
                                genElecInAcc->push_back(genDecayLVec[i]);
                                double dRMin = 999.9;
                                double matchPt = -999.9;
                                for(int j = 0; j < cutElecVecRecoOnly.size(); ++j)
                                {
                                    // difference in angle between gen and reco electrons
                                    double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutElecVecRecoOnly[j]);
                                    if(dR < dRMin)
                                    {
                                        dRMin = dR;
                                        matchPt = cutElecVecRecoOnly[j].Pt();
                                    }
                                }
                                if(dRMin < 0.02)
                                {
                                    genMatchElecInAcc->push_back(genDecayLVec[i]);
                                    genMatchElecInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                                }

                                dRMin = 999.9;
                                for(int j = 0; j < cutElecVec.size(); ++j)
                                {
                                    // difference in angle between gen and cut electrons
                                    double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutElecVec[j]);
                                    if(dR < dRMin)
                                    {
                                        dRMin = dR;
                                    }
                                }
                                if(dRMin < 0.02)
                                {
                                    genMatchIsoElecInAcc->push_back(genDecayLVec[i]);
                                }
                            }
                        }
                    }
                }
            }

            bool debug = false;

            data_t genZPt = -999.9, genZEta = -999.9, genZmass = -999.9, genZPhi;
            int nZ = 0;
            TLorentzVector genZ;
            int pdgIdZDec = 0;
            if(tr.checkBranch("GenPartTLV") && &genDecayLVec != nullptr)
            {
                for(int j = 0; j <  genDecayPdgIdVec.size(); ++j)
                {
                    //WARNING: make sure to use parentheses; otherwise == happens first and & second!!!
                    int maskedStatusFlag = (GenPart_statusFlags[j] & GENPARTMASK);
                    bool isGoodGenPart = (maskedStatusFlag == GENPARTMASK);
                    if(debug) printf("flag = 0x%x; maskedFlag = 0x%x; mask = 0x%x; pdgid = %d",GenPart_statusFlags[j], maskedStatusFlag, GENPARTMASK, genDecayPdgIdVec[j]);
                    if (isGoodGenPart)
                    {
                        if(abs(genDecayPdgIdVec[j]) == 23)
                        {
                            if (debug) printf(" - behold, a Z boson");
                            nZ++;
                            genZ = genDecayLVec[j];
                            genZPt = genDecayLVec[j].Pt();
                            genZEta = genDecayLVec[j].Eta();
                            genZPhi = genDecayLVec[j].Phi();
                            genZmass = genDecayLVec[j].M();
                        }
                    }
                    if(debug) std::cout << std::endl;
                }
                if (debug)
                {
                    std::cout << "nZ = " << nZ;
                    if(nZ > 1) std::cout << " - WARNING: MORE THAN 1 Z FOUND!!!";
                    std::cout << std::endl;
                }
                else
                {
                    if(nZ > 1) std::cout << " - WARNING: MORE THAN 1 Z FOUND!!! Set debug = true for more info." << std::endl;
                }

            }

            bool passDiMuTrig  = nTriggerMuons >= 2;

            const double zMassMin = 81.0;
            const double zMass    = 91.0;
            const double zMassMax = 101.0;

            double zMuMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
            TLorentzVector bestRecoMuZ;
            TLorentzVector Zrecopt;
            for(int i = 0; i < cutMuVec.size(); ++i)
            {
            Zrecopt =  muonsLVec[0]+muonsLVec[1];//cutMuVec[0] + cutMuVec[1];
                if(cutMuVec[i].Pt() < minMuPt) continue;
                for(int j = 0; j < i && j < cutMuVec.size(); ++j)
                {
                    if(cutMuVec[j].Pt() < minMuPt) continue;
                    double zm = (cutMuVec[i] + cutMuVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zMuMassCurrent - zMass))
                    {
                        bestRecoMuZ = cutMuVec[i] + cutMuVec[j];
                        zMuMassCurrent = zm;
                    }
                }
            }
            double zElecMassCurrent = 1.0e300;
            TLorentzVector bestRecoElecZ;
            for(int i = 0; i < cutElecVec.size(); ++i)
            {
                if(cutElecVec[i].Pt() < minElecPt) continue;
                for(int j = 0; j < i && j < cutElecVec.size(); ++j)
                {
                    if(cutElecVec[j].Pt() < minElecPt) continue;
                    double zm = (cutElecVec[i] + cutElecVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zElecMassCurrent - zMass))
                    {
                        bestRecoElecZ = cutElecVec[i] + cutElecVec[j];
                        zElecMassCurrent = zm;
                    }
                }
            }

            double zElMuMassCurrent = 1.0e300;
            TLorentzVector bestRecoElMuZ;
            for(int i = 0; i < cutMuVec.size(); ++i)
            {
                if(cutMuVec[i].Pt() < minMuPt) continue;
                for(int j = 0; j < cutElecVec.size(); ++j)
                {
                    if(cutElecVec[j].Pt() < minMuPt) continue;
                    double zm = (cutMuVec[i] + cutElecVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zElMuMassCurrent - zMass))
                    {
                        bestRecoElMuZ = cutMuVec[i] + cutElecVec[j];
                        zElMuMassCurrent = zm;
                    }
                }
            }

            TLorentzVector metV, metZ;
            metV.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

            TLorentzVector bestRecoZ = (fabs(bestRecoElecZ.M() - zMass) > fabs(bestRecoMuZ.M() - zMass)) ? (bestRecoMuZ) : (bestRecoElecZ);
            metZ.SetPtEtaPhiM(bestRecoZ.Pt(), 0.0, bestRecoZ.Phi(), 0.0);
            TLorentzVector cleanMet = metV + metZ;
            
            bool passDiMuSel   = passEleVeto  && (cutMuVec.size() == 2   && cutMuSummedCharge == 0   && cutMuVec[0].Pt() > highMuPt     && cutMuVec[1].Pt() > minMuPt);
            bool passDiElecSel = passMuonVeto && (cutElecVec.size() == 2 && cutElecSummedCharge == 0 && cutElecVec[0].Pt() > highElecPt && cutElecVec[1].Pt() > minElecPt);
            bool passElMuSel = (cutMuVec.size() == 1 && cutElecVec.size() == 1 && cutElecSummedCharge == -cutMuSummedCharge && cutMuVec[0].Pt() > highMuPt && cutElecVec[0].Pt() > minMuPt);

            bool passMuZinvSel         = passEleVeto && (cutMuVec.size() == 2   && cutMuSummedCharge == 0   && cutMuVec[0].Pt() > highMuPt     && cutMuVec[1].Pt() > minMuPt)     && (bestRecoMuZ.M() > zMassMin)   && (bestRecoMuZ.M() < zMassMax);
            bool passMuZinvSel_lowpt   = passEleVeto && (cutMuVec.size() == 2   && cutMuSummedCharge == 0   && cutMuVec[0].Pt() > minMuPt      && cutMuVec[1].Pt() > minMuPt)     && (bestRecoMuZ.M() > zMassMin)   && (bestRecoMuZ.M() < zMassMax);
            bool passElecZinvSel       = passMuonVeto && (cutElecVec.size() == 2 && cutElecSummedCharge == 0 && cutElecVec[0].Pt() > highElecPt && cutElecVec[1].Pt() > minElecPt) && (bestRecoElecZ.M() > zMassMin) && (bestRecoElecZ.M() < zMassMax);
            bool passElecZinvSel_lowpt = passMuonVeto && (cutElecVec.size() == 2 && cutElecSummedCharge == 0 && cutElecVec[0].Pt() > minElecPt  && cutElecVec[1].Pt() > minElecPt) && (bestRecoElecZ.M() > zMassMin) && (bestRecoElecZ.M() < zMassMax);
            bool passElMuZinvSel       = (cutMuVec.size() == 1 && cutElecVec.size() == 1 && cutElecSummedCharge == -cutMuSummedCharge && cutMuVec[0].Pt() > highMuPt && cutElecVec[0].Pt() > minMuPt) && (bestRecoElMuZ.M() > zMassMin) && (bestRecoElMuZ.M() < zMassMax);

            double genMuPt     = -999.9;
            double genMuEta    = -999.9;
            double cutMuPt1    = -999.9;
            double cutMuPt2    = -999.9;
            double cutMuEta1   = -999.9;
            double cutMuEta2   = -999.9;
            double cutElecPt1  = -999.9;
            double cutElecPt2  = -999.9;
            double cutElecEta1 = -999.9;
            double cutElecEta2 = -999.9;
            
            // print number of muons
            //printf("num gen mu: %d num mu: %d num cut mu: %d num cut mu reco only: %d\n", genMu->size(), muonsLVec.size(), cutMuVec.size(), cutMuVecRecoOnly.size());

            // print passMuZinvSel conditions
            bool printMuon = false;
            if (printMuon)
            {
                if (cutMuVec.size() > 0)
                {
                    printf("passEleVeto: %d ", passEleVeto);
                    printf("cutMuVec.size() == 2: %d ", cutMuVec.size() == 2);
                    printf("cutMuSummedCharge == 0: %d ", cutMuSummedCharge == 0);
                    printf("cutMuVec[0].Pt() > highMuPt: %d ", cutMuVec[0].Pt() > highMuPt);
                    printf("cutMuVec[1].Pt() > minMuPt: %d ", cutMuVec[1].Pt() > minMuPt);
                    printf("bestRecoMuZ.M() > zMassMin: %d ", bestRecoMuZ.M() > zMassMin);
                    printf("bestRecoMuZ.M() < zMassMax: %d ", bestRecoMuZ.M() < zMassMax);
                    if (passMuZinvSel)
                    {
                        printf(" --- passMuZinvSel");
                    }
                    printf("\n");
                }
                else
                {
                    printf("WARNING: no cut muons found in event.\n");
                }
            }


            if(genMu->size() >= 1)
            {
                genMuPt  = genMu->at(0)->Pt();
                genMuEta = genMu->at(0)->Eta();
            }
            if(cutMuVec.size() >= 1) 
            {
                cutMuPt1  = cutMuVec.at(0).Pt();
                cutMuEta1 = cutMuVec.at(0).Eta();
            }
            if(cutMuVec.size() >= 2) 
            {
                cutMuPt2  = cutMuVec.at(1).Pt();
                cutMuEta2 = cutMuVec.at(1).Eta();
            }
            if(cutElecVec.size() >= 1) 
            {
                cutElecPt1  = cutElecVec.at(0).Pt();
                cutElecEta1 = cutElecVec.at(0).Eta();
            }
            if(cutElecVec.size() >= 2) 
            {
                cutElecPt2  = cutElecVec.at(1).Pt();
                cutElecEta2 = cutElecVec.at(1).Eta();
            }

            double mindPhiMetJ = 999.9;
            int jc = 0;
            for(const TLorentzVector& jet : jetsLVec)
            {
                if(jc >= 3) break;
                jc++;
                mindPhiMetJ = std::min(mindPhiMetJ, fabs(ROOT::Math::VectorUtil::DeltaPhi(genZ, jet)));
            }

            data_t bestRecoZPt  = bestRecoZ.Pt();
            data_t metWithLL    = cleanMet.Pt();
            data_t metphiWithLL = cleanMet.Phi();
            data_t Zrecoptpt    = Zrecopt.Pt();
            tr.registerDerivedVar("bestRecoZPt", bestRecoZPt);
            tr.registerDerivedVar("bestRecoZM", bestRecoZ.M());
            tr.registerDerivedVar("metWithLL", metWithLL);
            tr.registerDerivedVar("metphiWithLL", metphiWithLL);
            tr.registerDerivedVar("cutMuPt1", cutMuPt1);
            tr.registerDerivedVar("cutMuPt2", cutMuPt2);
            tr.registerDerivedVar("cutMuEta1", cutMuEta1);
            tr.registerDerivedVar("cutMuEta2", cutMuEta2);
            tr.registerDerivedVar("cutElecPt1", cutElecPt1);
            tr.registerDerivedVar("cutElecPt2", cutElecPt2);
            tr.registerDerivedVar("cutElecEta1", cutElecEta1);
            tr.registerDerivedVar("cutElecEta2", cutElecEta2);
            tr.registerDerivedVar("mindPhiMetJ", mindPhiMetJ);
            tr.registerDerivedVar("ZPtRes", (bestRecoZPt - genZPt)/genZPt);
            tr.registerDerivedVar("ZEtaRes", bestRecoZ.Eta() - genZEta);
            tr.registerDerivedVar("ZPhiRes", bestRecoZ.Phi() - genZPhi);
            tr.registerDerivedVar("ZMRes", (bestRecoZ.M() - genZmass)/genZmass);
            tr.registerDerivedVec("genMu", genMu);
            const auto& genMu_test = tr.getVec<const TLorentzVector*>("genMu");
            tr.registerDerivedVar("ngenMu", static_cast<data_t>(genMu->size()));
            tr.registerDerivedVec("genMuInAcc", genMuInAcc);
            tr.registerDerivedVar("ngenMuInAcc", static_cast<data_t>(genMuInAcc->size()));
            tr.registerDerivedVec("genMatchMuInAcc", genMatchMuInAcc);
            tr.registerDerivedVec("genMatchMuInAccRes", genMatchMuInAccRes);
            tr.registerDerivedVec("genMatchIsoMuInAcc", genMatchIsoMuInAcc);
            tr.registerDerivedVar("ngenMatchMuInAcc", static_cast<data_t>(genMatchMuInAcc->size()));
            tr.registerDerivedVec("genElec", genElec);
            tr.registerDerivedVar("ngenElec", static_cast<data_t>(genElec->size()));
            tr.registerDerivedVec("genElecInAcc", genElecInAcc);
            tr.registerDerivedVar("ngenElecInAcc", static_cast<data_t>(genElecInAcc->size()));
            tr.registerDerivedVec("genMatchElecInAcc", genMatchElecInAcc);
            tr.registerDerivedVec("genMatchElecInAccRes", genMatchElecInAccRes);
            tr.registerDerivedVec("genMatchIsoElecInAcc", genMatchIsoElecInAcc);
            tr.registerDerivedVar("ngenMatchElecInAcc", static_cast<data_t>(genMatchElecInAcc->size()));
            tr.registerDerivedVar("genZPt", genZPt);
            tr.registerDerivedVar("genZEta", genZEta);
            tr.registerDerivedVar("genZPhi", genZPhi);
            tr.registerDerivedVar("genZmass", genZmass);
            tr.registerDerivedVar("pdgIdZDec", pdgIdZDec);
            tr.registerDerivedVar("passDiMuIsoTrig", passDiMuTrig);
            tr.registerDerivedVar("passSingleMu45", muTrigMu45);
            tr.registerDerivedVar("passDiMuSel", passDiMuSel);
            tr.registerDerivedVar("passDiElecSel", passDiElecSel);
            tr.registerDerivedVar("passElMuSel", passElMuSel);
            tr.registerDerivedVar("passMuZinvSel",         passMuZinvSel);
            tr.registerDerivedVar("passMuZinvSel_lowpt",   passMuZinvSel_lowpt);
            tr.registerDerivedVar("passElecZinvSel",       passElecZinvSel);
            tr.registerDerivedVar("passElecZinvSel_lowpt", passElecZinvSel_lowpt);
            tr.registerDerivedVar("passElMuZinvSel",       passElMuZinvSel);
            tr.registerDerivedVar("Zrecopt", Zrecoptpt);

        }

    public:
        LepInfo() : tr3(new TRandom3())
        {
        }

        void operator()(NTupleReader& tr)
        {
            lepInfo(tr);
        }
    };

}
#endif
