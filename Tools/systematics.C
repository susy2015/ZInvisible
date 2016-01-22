#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"
#include "Systematic.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
    using namespace std;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"plot",             no_argument, 0, 'p'},
        {"dataSets",   required_argument, 0, 'D'},
	{"luminosity", required_argument, 0, 'L'}
    };

    string dataSets = "ZJetsToNuNu";
    double lumi = AnaSamples::luminosity;

    while((opt = getopt_long(argc, argv, "pD:L:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'p':
            break;

        case 'D':
            dataSets = optarg;
            break;

	case 'L':
            lumi = atof(optarg);
            break;
        }
    }

    AnaSamples::SampleSet        ss(AnaSamples::fileDir, lumi);
    AnaSamples::SampleCollection sc(ss);

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {sc["DYJetsToLL"]};
        fileMap["ZJetsToNuNu"] = {sc["ZJetsToNuNu"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_400to600"]};
        fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
        fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
    }
    else
    {
        if(ss[dataSets] != ss.null())
        {
            fileMap[dataSets] = {ss[dataSets]};
            for(const auto& colls : ss[dataSets].getCollections())
            {
                fileMap[colls] = {ss[dataSets]};
            }
        }
        else if(sc[dataSets] != sc.null())
        {
            fileMap[dataSets] = {sc[dataSets]};
            int i = 0;
            for(const auto& fs : sc[dataSets])
            {
                fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
            }
        }
    }

    for(auto& vfs : fileMap)
    {
        for(auto& fs : vfs.second)
        {
            size_t start = fs.filePath.rfind('/');
            size_t stop  = fs.filePath.rfind('.');
            if(int(stop) - (int(start) + 1) > 0)
            {
                fs.treePath = fs.filePath.substr(start + 1, stop - start - 1);
                fs.filePath = "stupidFileList.txt";
                fs.readFileList();
            }
        }
    }

    RegisterFunctions *rf = new RegisterFunctionsSyst;

    vector<Plotter::HistSummary> vh;

    Plotter::DatasetSummary dsDY_mm( "DY#rightarrow#mu#mu", fileMap["DYJetsToLL"],  "", "");
    Plotter::DatasetSummary dsDY_ee( "DY#rightarrowee",     fileMap["DYJetsToLL"],  "", "");
    Plotter::DatasetSummary dsZ_nunu("Z#rightarrow#nu#nu",  fileMap["ZJetsToNuNu"], "", "");

    Plotter::DataCollection dcDY_mm_cleamMET( "single", "cleanMetPt", {dsDY_mm});
    Plotter::DataCollection dcDY_ee_cleamMET( "single", "cleanMetPt", {dsDY_ee});
    Plotter::DataCollection dcZ_nunu_cleamMET("single", "cleanMetPt", {dsZ_nunu});

    vh.push_back(PHS("cleanMet",     {dcDY_mm_cleamMET, dcDY_ee_cleamMET, dcZ_nunu_cleamMET}, {1, 2}, "",            100, 0, 1500,  true,  true,  "MET [GeV]",  "Norm Events"));
    vh.push_back(PHS("cleanMet_cut", {dcDY_mm_cleamMET, dcDY_ee_cleamMET, dcZ_nunu_cleamMET}, {1, 2}, "passMETZinv", 100, 0, 1500,  true,  true,  "MET [GeV]",  "Norm Events"));

    TF1 *fit = new TF1("fit","pol1");
    fit->SetParameters(0.98, 0.001);
    Systematic test("systWgtTest", "cleanMetPt", fit);
    test.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(test, std::placeholders::_1));

    //Met shape syst 
    TF1 *MET_fit = new TF1("MET_fit","pol1");
    //Chi2                      =      10.0185
    //NDf                       =            8
    //p0                        =      1.08688   +/-   0.0689733   
    //p1                        =  -0.00079217   +/-   0.000223161 
    MET_fit->SetParameters(1.08688, -0.00079217);
    Systematic METSyst("systWgtMET", "cleanMetPt", MET_fit);
    METSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(METSyst, std::placeholders::_1));

    //MT2 shape syst 
    TF1 *MT2_fit = new TF1("MT2_fit","pol1");
    //Chi2                      =      3.44991
    //NDf                       =            5
    //p0                        =      1.11303   +/-   0.0809109   
    //p1                        = -0.000670589   +/-   0.000259671 
    MT2_fit->SetParameters(1.11303, -0.000670589);
    Systematic MT2Syst("systWgtMT2", "best_had_brJet_MT2Zinv", MT2_fit);
    MT2Syst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(MT2Syst, std::placeholders::_1));

    //NT shape uncertainty
    double bins[] = {0, 1, 2, 8};
    TH1 *NT_hist = new TH1D("NT_hist","NT_hist", 3, bins);
    NT_hist->SetBinContent(1, 1.00401505518);
    NT_hist->SetBinContent(2, 0.981637373833);
    NT_hist->SetBinContent(3, 0.617594525252);

    Systematic NTSyst("systWgtNT", "nTopCandSortedCntZinv", NT_hist);
    NTSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NTSyst, std::placeholders::_1));

    //NB shape uncertainty
    double binsNB[] = {0, 1, 2, 3, 8};
    TH1 *NB_hist = new TH1D("NB_hist","NB_hist", 3, binsNB);
    NB_hist->SetBinContent(1, 0.986938983347);
    NB_hist->SetBinContent(2, 0.981886953393);
    NB_hist->SetBinContent(3, 0.610958135550);
    NB_hist->SetBinContent(4, 0.467353090212);

    Systematic NBSyst("systWgtNB", "cntCSVSZinv", NB_hist);
    NBSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NBSyst, std::placeholders::_1));

    Plotter::DataCollection dcDY_nunu_nJetWgtSyst( "single", {{"njSystWeightedSB", dsZ_nunu}, {"njSystUnweightedSB", dsZ_nunu}});
    vh.emplace_back(PHS("systNJetWgtStat", {dcDY_nunu_nJetWgtSyst}, {2, 1}, "", 45, 0, 45, true, true, "Search Bin", "Events"));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, true, "systematics.root");
    plotter.setLumi(lumi);
    plotter.setPlotDir("Syst_plots");
    plotter.setDoHists(true);
    plotter.setDoTuple(false);
    plotter.setPrintInterval(10000);
    plotter.setRegisterFunction(rf);
    plotter.read();
    plotter.saveHists();
    plotter.plot();
}
