#include "plot.h"

int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_SingleMuon-combined.root", "PEX0",       kBlack};

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

    plt.plot("dilepton/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/nJets", "Jets/Event", "Events", true);
    plt.plot("dilepton/nVertices", "primary vertices/event", "Events", true);
    plt.plot("dilepton/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("dilepton/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("dilepton/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("dilepton/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/bestTopPt", "dilepton/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/bestTopMass", "dilepton/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/bestTopEta", "dilepton/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/bestTopPt", "dilepton/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/bestTopMass", "dilepton/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/bestTopEta", "dilepton/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotSF("dilepton/bestTopPt", "dilepton/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("dilepton/bestTopMass", "dilepton/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("dilepton/bestTopEta", "dilepton/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("dilepton/bestTopPt", "dilepton/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("dilepton/bestTopMass", "dilepton/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("dilepton/bestTopEta", "dilepton/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "dilepton/METTagged",        "dilepton/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/HTTagged",         "dilepton/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/nJetsTagged",      "dilepton/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(true, "dilepton/nVerticesTagged",  "dilepton/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "dilepton/photonTagged",     "dilepton/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "dilepton/METTagged",        "dilepton/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/HTTagged",         "dilepton/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/nJetsTagged",      "dilepton/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "dilepton/nVerticesTagged",  "dilepton/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "dilepton/photonTagged",     "dilepton/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("dilepton/METTagged",       "dilepton/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("dilepton/HTTagged",        "dilepton/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("dilepton/nJetsTagged",     "dilepton/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("dilepton/nVerticesTagged", "dilepton/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("dilepton/photonTagged",     "dilepton/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("dilepton/METTagged",       "dilepton/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("dilepton/HTTagged",        "dilepton/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("dilepton/nJetsTagged",     "dilepton/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("dilepton/nVerticesTagged", "dilepton/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("dilepton/photonTagged",     "dilepton/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
