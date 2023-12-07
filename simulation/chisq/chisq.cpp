#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>

// datafiles
#define INPUT_FILE_O "../dat/outO_g_mix.dat"
#define INPUT_FILE_H "../dat/outH_g_mix.dat"
#define OUTPUT_FILE "../dat/chi2.dat"
#define OUTPUT_TREE "../dat/chi2_check.root"

// physical constants
#define pi 3.141592
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
#define DATASIZE 1e6

// alignment constants
#define Lambda_min 2e-10
#define Lambda_max 9e-10
#define TOTAL_LENGTH 1
#define DAQ_FREQ 62.5e+6 / 256

constexpr double bgO = 0;
constexpr double bgH = 0;
constexpr double offset = 0;

constexpr int daq_downsizing = 16;

constexpr double dt = daq_downsizing / DAQ_FREQ;
constexpr double dLMD = h / TOTAL_LENGTH / m * dt;
// constexpr int L_max = (int)((Lambda_max - Lambda_min) / dLMD);
constexpr int L_max = 10000;    // WARNING: L_max is changed

int chisq()
{
    // PREPARATIONS
    std::cout << "///// PREPARATIONS /////" << std::endl;

    // open the datafile
    std::cout << "// files open ///" << std::endl;
    std::ifstream inputFileO(INPUT_FILE_O);
    std::ifstream inputFileH(INPUT_FILE_H);
    std::ofstream output(OUTPUT_FILE);
    if (!inputFileO.is_open() || !inputFileH.is_open() || !output.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }
    else
    {
        std::cout << "Both files are opened." << std::endl;
    }

    // initialize the histograms
    TH1F *histO = new TH1F("histO", "", L_max, Lambda_min, Lambda_max);
    TH1F *histH = new TH1F("histH", "", L_max, Lambda_min, Lambda_max);

    // data input
    std::cout << "/// data input ///" << std::endl;
    double lambda;
    int count = 0;
    while (inputFileO >> lambda)
    {
        histO->Fill(lambda);
        count++;
        if(count>DATASIZE)break;
    }
    count = 0;
    std::cout << "O file has all input." << std::endl;
    while (inputFileH >> lambda)
    {
        histH->Fill(lambda);
        count++;
        if (count > DATASIZE)break;
    }
    std::cout << "H file has all input." << std::endl;


    // DATA CALCULATIONS
    std::cout << "///// DATA CALCULATIONS /////" << std::endl;

    TH1F *h1 = (TH1F *)histO->Clone();
    TH1F *h2 = (TH1F *)histH->Clone();
    h1->Add(histH, -1);
    h2->Add(histO, 1);
    h1->Divide(h2);

    // FITTINGS
    std::cout << "///// FITTINGS /////" << std::endl;
    // k (initial condition of wavelength) loop
    double k, chiSquare;
    int point = 0;
    TGraph *gr = new TGraph();
    for (k = 1e11; k < 1e12; k += 1e9)
    {
        TF1 *f_propto = new TF1("f_propto", "cos([0] * x + [1])");
        f_propto->SetParNames("k_{+}", "#theta_{0}");
        f_propto->SetParameters(k, 0);
        h1->Fit("f_propto", "", "", 7e-10, Lambda_max);

        chiSquare = f_propto->GetChisquare();
        output << k << " " << f_propto->GetParameter(0) << " " << chiSquare << std::endl;
        gr->SetPoint(point, k, chiSquare);
        point++;

        delete f_propto;
    }

    // Set the graph attributes
    gr->SetTitle("Chi-Square Values vs Initial Conditions of k");
    gr->GetXaxis()->SetTitle("Initial Condition of k");
    gr->GetYaxis()->SetTitle("Chi-Square Value");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kBlue);

    // Draw the graph
    gr->Draw("AP");

    delete gr;
    return 0;
}
