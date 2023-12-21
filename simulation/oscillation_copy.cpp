#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>

// datafiles
#define INPUT_TREE "./dat/theoretical/noref/g_main.root"

// physical constants
#define pi M_PI
#define m 1.6749e-27
#define h 6.62607015e-34
#define DATASIZE 1e6
constexpr double g = 9.8;

// alignment constants
constexpr double Lambda_min = 2e-10;
constexpr double Lambda_max = 12e-10;
constexpr double total_length = 1.;
constexpr double daq_freq = 62.5e6;
constexpr int daq_downsizing = 256;

int oscillation()
{
    // PREPARATIONS

    // alignment calculations & settings
    double dt = 1 / (daq_freq / daq_downsizing);
    double dLambda = h / total_length / m * dt;
    int L_max = (int)((Lambda_max - Lambda_min) / dLambda);
    L_max = 10000;
    std::cout << "///// ALIGNMENT CALCULATIONS /////" << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "dLambda = " << dLambda << std::endl;
    std::cout << "L_max = " << L_max << std::endl;

    // open the datafile
    std::cout << "///// FILE OPEN /////" << std::endl;
    std::ifstream inputFileO(INPUT_FILE_O);
    std::ifstream inputFileH(INPUT_FILE_H);
    if (!inputFileO.is_open() || !inputFileH.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }
    else
    {
        std::cout << "FILE OPENED" << std::endl;
    }

    // Initialize the histograms
    TH1F *histO = new TH1F("histO", "", L_max, Lambda_min, Lambda_max);
    TH1F *histH = new TH1F("histH", "", L_max, Lambda_min, Lambda_max);

    // DATA INPUT
    double lambda;
    int count = 0;
    std::cout << "///// DATA INPUT /////" << std::endl;
    while (inputFileO >> lambda)
    {
        histO->Fill(lambda);
        count++;
        if (count > DATASIZE)break;
    }
    count = 0;
    std::cout << "O file has all input" << std::endl;
    while (inputFileH >> lambda)
    {
        histH->Fill(lambda);
        count++;
        if (count > DATASIZE)break;
    }
    std::cout << "H file has all input" << std::endl;

    // DATA CALCULATIONS
    std::cout << "///// DATA CALCULATIONS /////" << std::endl;
    // H-O/H+O
    TH1F *h1 = (TH1F *)histO->Clone();
    TH1F *h2 = (TH1F *)histH->Clone();
    h1->Add(histH, -1);
    h2->Add(histO, 1);
    h1->Divide(h2);
    h1->Rebin(5);

    // FITTINGS
    std::cout << "///// FITTINGS /////" << std::endl;
    TF1 *f_propto = new TF1("f_propto", "([0] + [1] * x + [2] * x * x) * cos([3] * x + [4])");
    f_propto->SetParNames("k_{+}", "#theta_{0}");
    f_propto->SetParameters(1, 0, 0, 5.8e+11, 0);
    f_propto->SetParLimits(3, 1e+10, 1e+12);
    // f_propto->SetParLimits(1, -pi, pi);
    // TF1 *f = new TF1("f", "[0] * cos([1] / x + [2] * x + [3]) + [4]");
    // f->SetParNames("A", "k_{-1}", "k_{+1}", "#theta_{0}", "BG");
    // f->SetParameters(1, 0, 2 * pi * 1e+11, 0, 0);
    // f->SetParLimits(0, 1., 1.);
    // f->SetParLimits(1, 0., 0.);
    // f->SetParLimits(2, 1e+10, 1e+12);
    // f->SetParLimits(3, -pi, pi);
    // f->SetParLimits(4, 0, 0);
    h1->Fit("f_propto", "", "", 7e-10, Lambda_max);

    // GRAPH SETTINGS
    h1->Draw();

    // Clean up allocated memory
    delete histO;
    delete histH;
    // delete h1;
    // delete h2;

    std::cout << "///// COMPLETED /////" << std::endl;

    return 0;
}
