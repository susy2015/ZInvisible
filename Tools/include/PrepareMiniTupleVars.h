#ifndef PREPAREMINITUPLEVARS
#define PREPAREMINITUPLEVARS

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

#include "TopTagger/TopTagger/include/TopObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <vector>
#include <iostream>
#include <string>
#include <set>

namespace plotterFunctions
{
    class PrepareMiniTupleVars
    {
    private:

        static const int BIT_PASSLEPTVETO             = 0x00000001;
        static const int BIT_PASSMUONVETO             = 0x00000002;
        static const int BIT_PASSELEVETO              = 0x00000004;
        static const int BIT_PASSISOTRKVETO           = 0x00000008;
        static const int BIT_PASSNJETS                = 0x00000010;
        static const int BIT_PASSDPHIS                = 0x00000020;
        static const int BIT_PASSBJETS                = 0x00000040;
        static const int BIT_PASSMET                  = 0x00000080;
        static const int BIT_PASSMT2                  = 0x00000100;
        static const int BIT_PASSHT                   = 0x00000200;
        static const int BIT_PASSTAGGER               = 0x00000400;
        static const int BIT_PASSNOISEEVENTFILTER     = 0x00000800;
        static const int BIT_PASSBASELINE             = 0x00001000;
        static const int BIT_PASSBASELINENOTAGMT2     = 0x00002000;
        static const int BIT_PASSBASELINENOTAG        = 0x00004000;
        static const int BIT_PASSLEPTVETOZINV         = 0x00008000;
        static const int BIT_PASSMUONVETOZINV         = 0x00010000;
        static const int BIT_PASSELEVETOZINV          = 0x00020000;
        static const int BIT_PASSISOTRKVETOZINV       = 0x00040000;
        static const int BIT_PASSNJETSZINV            = 0x00080000;
        static const int BIT_PASSDPHISZINV            = 0x00100000;
        static const int BIT_PASSBJETSZINV            = 0x00200000;
        static const int BIT_PASSMETZINV              = 0x00400000;
        static const int BIT_PASSMT2ZINV              = 0x00800000;
        static const int BIT_PASSHTZINV               = 0x01000000;
        static const int BIT_PASSTAGGERZINV           = 0x02000000;
        static const int BIT_PASSNOISEEVENTFILTERZINV = 0x04000000;
        static const int BIT_PASSBASELINEZINV         = 0x08000000;
        static const int BIT_PASSBASELINENOTAGMT2ZINV = 0x10000000;
        static const int BIT_PASSBASELINENOTAGZINV    = 0x20000000;
        static const int BIT_PASSMUZINVSEL            = 0x40000000;
        static const int BIT_PASSELMUZINVSEL          = 0x80000000;

        bool pack_;

        void pack(NTupleReader& tr)
        {
            // standard
            const auto& passLeptVeto =              tr.getVar<bool>("passLeptVeto");
            const auto& passMuonVeto =              tr.getVar<bool>("passMuonVeto");
            const auto& passEleVeto =               tr.getVar<bool>("passEleVeto");
            const auto& passIsoTrkVeto =            tr.getVar<bool>("passIsoTrkVeto");
            const auto& passnJets =                 tr.getVar<bool>("passnJets");
            const auto& passdPhis =                 tr.getVar<bool>("passdPhis");
            const auto& passBJets =                 tr.getVar<bool>("passBJets");
            const auto& passMET =                   tr.getVar<bool>("passMET");
            const auto& passMT2 =                   tr.getVar<bool>("passMT2");
            const auto& passHT =                    tr.getVar<bool>("passHT");
            const auto& passTagger =                tr.getVar<bool>("passTagger");
            const auto& passNoiseEventFilter =      tr.getVar<bool>("passNoiseEventFilter");
            const auto& passBaseline =              tr.getVar<bool>("passBaseline");
            const auto& passBaselineNoTagMT2 =      tr.getVar<bool>("passBaselineNoTagMT2");
            const auto& passBaselineNoTag =         tr.getVar<bool>("passBaselineNoTag");
            // Zinv
            const auto& passLeptVetoZinv =          tr.getVar<bool>("passLeptVetoZinv");
            const auto& passMuonVetoZinv =          tr.getVar<bool>("passMuonVetoZinv");
            const auto& passEleVetoZinv =           tr.getVar<bool>("passEleVetoZinv");
            const auto& passIsoTrkVetoZinv =        tr.getVar<bool>("passIsoTrkVetoZinv");
            const auto& passnJetsZinv =             tr.getVar<bool>("passnJetsZinv");
            const auto& passdPhisZinv =             tr.getVar<bool>("passdPhisZinv");
            const auto& passBJetsZinv =             tr.getVar<bool>("passBJetsZinv");
            const auto& passMETZinv =               tr.getVar<bool>("passMETZinv");
            const auto& passMT2Zinv =               tr.getVar<bool>("passMT2Zinv");
            const auto& passHTZinv =                tr.getVar<bool>("passHTZinv");
            const auto& passTaggerZinv =            tr.getVar<bool>("passTaggerZinv");
            const auto& passNoiseEventFilterZinv =  tr.getVar<bool>("passNoiseEventFilterZinv");
            const auto& passBaselineZinv =          tr.getVar<bool>("passBaselineZinv");
            const auto& passBaselineNoTagMT2Zinv =  tr.getVar<bool>("passBaselineNoTagMT2Zinv");
            const auto& passBaselineNoTagZinv =     tr.getVar<bool>("passBaselineNoTagZinv");
            // ZinvSel
            const auto& passMuZinvSel =             tr.getVar<bool>("passMuZinvSel");
            const auto& passElMuZinvSel =           tr.getVar<bool>("passElMuZinvSel");

            int cuts = 0;

            // standard
            if(passLeptVeto)         cuts |= BIT_PASSLEPTVETO;
            if(passMuonVeto)         cuts |= BIT_PASSMUONVETO;
            if(passEleVeto)          cuts |= BIT_PASSELEVETO;
            if(passIsoTrkVeto)       cuts |= BIT_PASSISOTRKVETO;
            if(passnJets)            cuts |= BIT_PASSNJETS;
            if(passdPhis)            cuts |= BIT_PASSDPHIS;
            if(passBJets)            cuts |= BIT_PASSBJETS;
            if(passMET)              cuts |= BIT_PASSMET;
            if(passMT2)              cuts |= BIT_PASSMT2;
            if(passHT)               cuts |= BIT_PASSHT;
            if(passTagger)           cuts |= BIT_PASSTAGGER;
            if(passNoiseEventFilter) cuts |= BIT_PASSNOISEEVENTFILTER;
            if(passBaseline)         cuts |= BIT_PASSBASELINE;
            if(passBaselineNoTagMT2) cuts |= BIT_PASSBASELINENOTAGMT2;
            if(passBaselineNoTag)    cuts |= BIT_PASSBASELINENOTAG;
            // Zinv
            if(passLeptVetoZinv)         cuts |= BIT_PASSLEPTVETOZINV;
            if(passMuonVetoZinv)         cuts |= BIT_PASSMUONVETOZINV;
            if(passEleVetoZinv)          cuts |= BIT_PASSELEVETOZINV;
            if(passIsoTrkVetoZinv)       cuts |= BIT_PASSISOTRKVETOZINV;
            if(passnJetsZinv)            cuts |= BIT_PASSNJETSZINV;
            if(passdPhisZinv)            cuts |= BIT_PASSDPHISZINV;
            if(passBJetsZinv)            cuts |= BIT_PASSBJETSZINV;
            if(passMETZinv)              cuts |= BIT_PASSMETZINV;
            if(passMT2Zinv)              cuts |= BIT_PASSMT2ZINV;
            if(passHTZinv)               cuts |= BIT_PASSHTZINV;
            if(passTaggerZinv)           cuts |= BIT_PASSTAGGERZINV;
            if(passNoiseEventFilterZinv) cuts |= BIT_PASSNOISEEVENTFILTERZINV;
            if(passBaselineZinv)         cuts |= BIT_PASSBASELINEZINV;
            if(passBaselineNoTagMT2Zinv) cuts |= BIT_PASSBASELINENOTAGMT2ZINV;
            if(passBaselineNoTagZinv)    cuts |= BIT_PASSBASELINENOTAGZINV;
            // ZinvSel
            if(passMuZinvSel)            cuts |= BIT_PASSMUZINVSEL;
            if(passElMuZinvSel)          cuts |= BIT_PASSELMUZINVSEL;

            tr.registerDerivedVar("cuts", cuts);
        }

        void unpack(NTupleReader& tr)
        {
            const int& cuts = tr.getVar<int>("cuts");
            // standard
            tr.registerDerivedVar("passLeptVeto",               static_cast<bool>(cuts & BIT_PASSLEPTVETO));
            tr.registerDerivedVar("passMuonVeto",               static_cast<bool>(cuts & BIT_PASSMUONVETO));
            tr.registerDerivedVar("passEleVeto",                static_cast<bool>(cuts & BIT_PASSELEVETO));
            tr.registerDerivedVar("passIsoTrkVeto",             static_cast<bool>(cuts & BIT_PASSISOTRKVETO));
            tr.registerDerivedVar("passnJets",                  static_cast<bool>(cuts & BIT_PASSNJETS));
            tr.registerDerivedVar("passdPhis",                  static_cast<bool>(cuts & BIT_PASSDPHIS));
            tr.registerDerivedVar("passBJets",                  static_cast<bool>(cuts & BIT_PASSBJETS));
            tr.registerDerivedVar("passMET",                    static_cast<bool>(cuts & BIT_PASSMET));
            tr.registerDerivedVar("passMT2",                    static_cast<bool>(cuts & BIT_PASSMT2));
            tr.registerDerivedVar("passHT",                     static_cast<bool>(cuts & BIT_PASSHT));
            tr.registerDerivedVar("passTagger",                 static_cast<bool>(cuts & BIT_PASSTAGGER));
            tr.registerDerivedVar("passNoiseEventFilter",       static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTER));
            tr.registerDerivedVar("passBaseline",               static_cast<bool>(cuts & BIT_PASSBASELINE));
            tr.registerDerivedVar("passBaselineNoTagMT2",       static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2));
            tr.registerDerivedVar("passBaselineNoTag",          static_cast<bool>(cuts & BIT_PASSBASELINENOTAG));
            // Zinv
            tr.registerDerivedVar("passLeptVetoZinv",           static_cast<bool>(cuts & BIT_PASSLEPTVETOZINV));
            tr.registerDerivedVar("passMuonVetoZinv",           static_cast<bool>(cuts & BIT_PASSMUONVETOZINV));
            tr.registerDerivedVar("passEleVetoZinv",            static_cast<bool>(cuts & BIT_PASSELEVETOZINV));
            tr.registerDerivedVar("passIsoTrkVetoZinv",         static_cast<bool>(cuts & BIT_PASSISOTRKVETOZINV));
            tr.registerDerivedVar("passnJetsZinv",              static_cast<bool>(cuts & BIT_PASSNJETSZINV));
            tr.registerDerivedVar("passdPhisZinv",              static_cast<bool>(cuts & BIT_PASSDPHISZINV));
            tr.registerDerivedVar("passBJetsZinv",              static_cast<bool>(cuts & BIT_PASSBJETSZINV));
            tr.registerDerivedVar("passMETZinv",                static_cast<bool>(cuts & BIT_PASSMETZINV));
            tr.registerDerivedVar("passMT2Zinv",                static_cast<bool>(cuts & BIT_PASSMT2ZINV));
            tr.registerDerivedVar("passHTZinv",                 static_cast<bool>(cuts & BIT_PASSHTZINV));
            tr.registerDerivedVar("passTaggerZinv",             static_cast<bool>(cuts & BIT_PASSTAGGERZINV));
            tr.registerDerivedVar("passNoiseEventFilterZinv",   static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTERZINV));
            tr.registerDerivedVar("passBaselineZinv",           static_cast<bool>(cuts & BIT_PASSBASELINEZINV));
            tr.registerDerivedVar("passBaselineNoTagMT2Zinv",   static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2ZINV));
            tr.registerDerivedVar("passBaselineNoTagZinv",      static_cast<bool>(cuts & BIT_PASSBASELINENOTAGZINV));
            // ZinvSel
            tr.registerDerivedVar("passMuZinvSel",              static_cast<bool>(cuts & BIT_PASSMUZINVSEL));
            tr.registerDerivedVar("passElMuZinvSel",            static_cast<bool>(cuts & BIT_PASSELMUZINVSEL));
        }

    public:
        PrepareMiniTupleVars(bool pack)
        {
            pack_ = pack;
        }

        void operator()(NTupleReader& tr)
        {
            if(pack_) pack(tr);
            else      unpack(tr);
        }
    };
}

#endif
