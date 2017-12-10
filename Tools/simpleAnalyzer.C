#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "BTagCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"

#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"


int main()
{
    TH1::AddDirectory(false);

    TH1* hExample = new TH1D("exmaple","exanple",100,0, 1000);

    AnaSamples::SampleSet        ss(AnaSamples::fileDir, AnaSamples::luminosity);
    AnaSamples::SampleCollection sc(ss);

    try
    {
        //for(auto& fs : sc["DYJetsToLL"])
        for(auto& fs : sc["Signal_T2tt_mStop850_mLSP100"])
        {
            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t);

            std::cout << "Tree: " << fs.treePath << std::endl;
            //std::cout << "sigma*lumi: " << fs.getWeight() << std::endl;

            
            BaselineVessel myBLV(*static_cast<NTupleReader*>(nullptr), "TopTag", "");
            plotterFunctions::PrepareTopVars prepareTopVars;
            plotterFunctions::TriggerInfo triggerInfo(false, true);
            
            BTagCorrector bTagCorrector("allINone_bTagEff.root", "", false);
            ISRCorrector ISRcorrector("allINone_ISRJets.root","","");   
            Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");

            NTupleReader tr(t);
            tr.registerFunction(myBLV);
            tr.registerFunction(prepareTopVars);
            tr.registerFunction(triggerInfo);
            tr.registerFunction(bTagCorrector);
            tr.registerFunction(ISRcorrector);
            tr.registerFunction(pileup);

            double fileWgt = fs.getWeight();

            while(tr.getNextEvent())
            {
                if(tr.getEvtNum() % 100 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

                const double& met    = tr.getVar<double>("met");

                if(tr.getVar<bool>("passNoiseEventFilterTopTag") && tr.getVar<bool>("passSingleLep20") /*and so on*/)
                {
                    hExample->Fill(met, fileWgt);
                }
            }
        }
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }

}
