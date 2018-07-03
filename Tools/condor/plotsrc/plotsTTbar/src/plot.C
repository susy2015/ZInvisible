#include "plot.h"


int main()
{   
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_MET-combined.root", "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
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
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //These are the background that we want to subtract from the data
    std::vector<histInfo> unwantedEntries = {
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring}, 
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3}
    };

    //These are the entires that have the objects we want (hadronically decaying tops)
    std::vector<histInfo> wantedEntries = {
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //vector summarizing signal histograms to include in the plot
//    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
//    };

    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));
    Plotter pltSub(std::move(data), std::move(unwantedEntries), std::move(wantedEntries)); //Plotter to make background subtracted plots

    plt.plot("ttbar/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbar/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/nJets", "Jets/Event", "Events", true);
    plt.plot("ttbar/nVertices", "primary vertices/event", "Events", true);
    plt.plot("ttbar/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("ttbar/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);

    pltSub.plot("ttbar/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/nJets", "Jets/Event", "Events", true);
    pltSub.plot("ttbar/nVertices", "primary vertices/event", "Events", true);
    pltSub.plot("ttbar/photon", "Photon p_{T} [GeV]", "Events", true);
    pltSub.plot("ttbar/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("ttbar/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbar/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbar/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    pltSub.plot("ttbar/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbar/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    plt.plotEff(true, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    pltSub.plotEff(true, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotEff(false, "ttbar/genTopMatchPt", "ttbar/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/genTopMatchMass", "ttbar/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/genTopMatchEta", "ttbar/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    pltSub.plotEff(false, "ttbar/genTopMatchPt", "ttbar/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/genTopMatchMass", "ttbar/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/genTopMatchEta", "ttbar/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    plt.plotSF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    pltSub.plotSF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    pltSub.plotOnlySF("ttbar/bestTopPt", "ttbar/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/bestTopMass", "ttbar/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/bestTopEta", "ttbar/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(true, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(true, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(false, "ttbar/METTagged",        "ttbar/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/HTTagged",         "ttbar/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/nJetsTagged",      "ttbar/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(false, "ttbar/nVerticesTagged",  "ttbar/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    pltSub.plotSF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    pltSub.plotSF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

    pltSub.plotOnlySF("ttbar/METTagged",       "ttbar/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/HTTagged",        "ttbar/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/nJetsTagged",     "ttbar/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    pltSub.plotOnlySF("ttbar/nVerticesTagged", "ttbar/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbar/photonTagged",     "ttbar/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
