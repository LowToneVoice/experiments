#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>

// datafiles
#define INPUT_FILE_O "./dat/montecarlo/ref/lambda/g_main_O.dat"
#define INPUT_FILE_H "./dat/montecarlo/ref/lambda/g_main_H.dat"
#define OUTPUT_FILE "./oscil_graph/montecarlo/ref/lambda/g_main.pdf"
#define OUTPUT_FILE_ZOOM "./oscil_graph/montecarlo/ref/lambda/g_main_zoom.pdf"

// PHYS & MATH CONSTANTS
#define h 6.62607015e-34
#define pi M_PI
#define J_per_eV 1.6022e-19
#define m 1.675e-27
constexpr double g = 9.8;
constexpr double hbar = h / (2 * pi);

// alignment constants
constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 12e-10;
constexpr int daq_freq = 62.5e6;

// alignment variables
constexpr double lambda_min_used = 8e-10;
constexpr double lambda_max_used = 8.5e-10;
constexpr double theta_min = .2 * pi / 180;
constexpr double theta_max = 1.5 * pi / 180;
constexpr double total_length = 1.;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;

constexpr int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);

int oscillation_lambda()
{
    // open the datafile
    std::ifstream inputFileO(INPUT_FILE_O);
    std::ifstream inputFileH(INPUT_FILE_H);
    if (!inputFileO.is_open() || !inputFileH.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }

    // make histogram
    TH1F *histO = new TH1F("Obeam", "Oscil", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH = new TH1F("Hbeam", "Oscil", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histO_zoom = new TH1F("Obeam", "Oscil", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1F *histH_zoom = new TH1F("Hbeam", "Oscil", N_loop_lambda, lambda_min_used, lambda_max_used);

    double lambda, theta;

    while (inputFileO >> lambda)
    {
        histO->Fill(lambda);
        histO_zoom->Fill(lambda);
    }
    while (inputFileH >> lambda)
    {
        histH->Fill(lambda);
        histH_zoom->Fill(lambda);
    }

    // data calculations
    TH1F *h1 = (TH1F *)histO->Clone();
    TH1F *h2 = (TH1F *)histH->Clone();
    h1->Add(histH, -1);
    h2->Add(histO, 1);
    h1->Divide(h2);
    TH1F *h1_zoom = (TH1F *)histO_zoom->Clone();
    TH1F *h2_zoom = (TH1F *)histH_zoom->Clone();
    h1_zoom->Add(histH_zoom, -1);
    h2_zoom->Add(histO_zoom, 1);
    h1_zoom->Divide(h2_zoom);

    // fittings
    TF1 *f = new TF1("f", "([0] + [1] * x + [2]) * cos([3] * x + [4]) + [5]");
    f->SetParNames("A0", "A1", "A2", "k_{+}", "#theta_{0}", "BG");
    f->SetParameters(1, 0, 0, 5.8e11, 0, 0);

    // draw histogram
    TCanvas *c = new TCanvas("c", "oscillation", 700, 500);
    h1->Fit("f", "", "", lambda_min, lambda_max);
    h1->Draw("c");
    c->Print(OUTPUT_FILE);

    TCanvas *c1 = new TCanvas("c1", "oscillation", 700, 500);
    h1_zoom->Fit("f", "", "", lambda_min_used, lambda_max_used);
    h1_zoom->Draw("c1");
    c1->Print(OUTPUT_FILE_ZOOM);

    delete histO;
    delete histH;
    delete histO_zoom;
    delete histH_zoom;

    return 0;
}
