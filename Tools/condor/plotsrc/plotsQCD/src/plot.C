#include "plot.h"


int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_JetHT-combined.root", "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
//        {"GJets",         "GJets-combined.root",      "hist", kViolet},
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
//        {"DYJets",        "DYJets-combined.root",     "hist", kCyan},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring},
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta},
    };

    //vector containing object we wish to keep in the MC
    std::vector<histInfo> wantedEntries = {
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3}
    };
	
	    //vector containing objects we wish to subtract from the data (background subtraction)
    std::vector<histInfo> unwantedEntries = {
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"tW",            "tW-combined.root",         "hist", kSpring},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    makeplots("QCD",data,bgEntries,wantedEntries,unwantedEntries);
    return 0;
    //vector summarizing signal histograms to include in the plot
//    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
//    };

    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));
    Plotter pltSub(std::move(data), std::move(unwantedEntries), std::move(wantedEntries));
	
    plt.plot("QCD/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/nJets", "Jets/Event", "Events", true);
    plt.plot("QCD/nVertices", "primary vertices/event", "Events", true);
    plt.plot("QCD/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("QCD/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("QCD/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("QCD/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotEff(false, "QCD/genTopMatchPt", "QCD/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/genTopMatchMass", "QCD/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/genTopMatchEta", "QCD/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    pltSub.plotEff(false, "QCD/genTopMatchPt", "QCD/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/genTopMatchMass", "QCD/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/genTopMatchEta", "QCD/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    plt.plotSF("QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "QCD/METTagged",        "QCD/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/HTTagged",         "QCD/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/nJetsTagged",      "QCD/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(true, "QCD/nVerticesTagged",  "QCD/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "QCD/METTagged",        "QCD/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/HTTagged",         "QCD/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/nJetsTagged",      "QCD/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "QCD/nVerticesTagged",  "QCD/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("QCD/METTagged",       "QCD/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("QCD/HTTagged",        "QCD/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("QCD/nJetsTagged",     "QCD/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("QCD/nVerticesTagged", "QCD/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("QCD/METTagged",       "QCD/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("QCD/HTTagged",        "QCD/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("QCD/nJetsTagged",     "QCD/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("QCD/nVerticesTagged", "QCD/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);
	
    pltSub.plot("QCD/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/nJets", "Jets/Event", "Events", true);
    pltSub.plot("QCD/nVertices", "primary vertices/event", "Events", true);
    pltSub.plot("QCD/photon", "Photon p_{T} [GeV]", "Events", true);
    pltSub.plot("QCD/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    pltSub.plot("QCD/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("QCD/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    pltSub.plotSF("QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    pltSub.plotOnlySF("QCD/bestTopPt", "QCD/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("QCD/bestTopMass", "QCD/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("QCD/bestTopEta", "QCD/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    pltSub.plotEff(true, "QCD/METTagged",        "QCD/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/HTTagged",         "QCD/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/nJetsTagged",      "QCD/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(true, "QCD/nVerticesTagged",  "QCD/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(false, "QCD/METTagged",        "QCD/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/HTTagged",         "QCD/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/nJetsTagged",      "QCD/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(false, "QCD/nVerticesTagged",  "QCD/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotSF("QCD/METTagged",       "QCD/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("QCD/HTTagged",        "QCD/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("QCD/nJetsTagged",     "QCD/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    pltSub.plotSF("QCD/nVerticesTagged", "QCD/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    pltSub.plotOnlySF("QCD/METTagged",       "QCD/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("QCD/HTTagged",        "QCD/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("QCD/nJetsTagged",     "QCD/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    pltSub.plotOnlySF("QCD/nVerticesTagged", "QCD/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("QCD/photonTagged",     "QCD/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);	

}
