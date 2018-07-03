#include "plot.h"

int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_SinglePhoton-combined.root", "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
        {"GJets",         "GJets-combined.root",      "hist", kViolet},
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
        {"DYJets",        "DYJets-combined.root",     "hist", kCyan},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring},
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta},
    };

    //vector summarizing signal histograms to include in the plot
//    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
//    };

    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));

    plt.plot("photon/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/nJets", "Jets/Event", "Events", true);
    plt.plot("photon/nVertices", "primary vertices/event", "Events", true);
    plt.plot("photon/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("photon/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("photon/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("photon/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/bestTopPt", "photon/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/bestTopMass", "photon/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/bestTopEta", "photon/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/bestTopPt", "photon/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/bestTopMass", "photon/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/bestTopEta", "photon/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotEff(false, "photon/genTopMatchPt", "photon/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/genTopMatchMass", "photon/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/genTopMatchEta", "photon/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

//    pltSub.plotEff(false, "photon/genTopMatchPt", "photon/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
//    pltSub.plotEff(false, "photon/genTopMatchMass", "photon/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
//    pltSub.plotEff(false, "photon/genTopMatchEta", "photon/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    plt.plotSF("photon/bestTopPt", "photon/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("photon/bestTopMass", "photon/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("photon/bestTopEta", "photon/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("photon/bestTopPt", "photon/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("photon/bestTopMass", "photon/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("photon/bestTopEta", "photon/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "photon/METTagged",        "photon/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/HTTagged",         "photon/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/nJetsTagged",      "photon/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(true, "photon/nVerticesTagged",  "photon/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "photon/photonTagged",     "photon/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "photon/METTagged",        "photon/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/HTTagged",         "photon/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/nJetsTagged",      "photon/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "photon/nVerticesTagged",  "photon/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "photon/photonTagged",     "photon/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("photon/METTagged",       "photon/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("photon/HTTagged",        "photon/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("photon/nJetsTagged",     "photon/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("photon/nVerticesTagged", "photon/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("photon/photonTagged",     "photon/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("photon/METTagged",       "photon/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("photon/HTTagged",        "photon/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("photon/nJetsTagged",     "photon/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("photon/nVerticesTagged", "photon/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("photon/photonTagged",     "photon/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
