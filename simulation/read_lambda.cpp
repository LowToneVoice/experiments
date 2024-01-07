#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>

// FILES
#define INPUT_TREE "./dat/montecarlo/ref/lambda/g_mix.root"
#define INPUT_O "./dat/montecarlo/ref/lambda/g_mix_O.dat"
#define INPUT_H "./dat/montecarlo/ref/lambda/g_mix_H.dat"
#define OUTPUT "./beam_count/montecarlo/ref/lambda/g_mix.pdf"
#define OUTPUT_ZOOM "./beam_count/montecarlo/ref/lambda/g_mix_zoom.pdf"

#define MAX_DATA_SIZE 1e+6

// CHANNELS
#define H_TDC 0
#define O_TDC 1
#define H_ADC 2
#define O_ADC 3

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

int read_lambda(double theta_output = 1.05 * pi / 180)
{
    // open the datafile
    std::ifstream inputFileO(INPUT_O);
    std::ifstream inputFileH(INPUT_H);
    if (!inputFileO.is_open() || !inputFileH.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }

    // make histogram
    TH1F *histO = new TH1F("Obeam", "Obeam", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH = new TH1F("Hbeam", "Hbeam", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histO_zoom = new TH1F("Obeam", "Obeam", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1F *histH_zoom = new TH1F("Hbeam", "Hbeam", N_loop_lambda, lambda_min_used, lambda_max_used);

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
    // TGraph *gr = new TGraph(L_max, lambdas, values);
    // gr->Draw("L");

    // TF1 *fHG = new TF1("fH", "[ 0 ] * cos([ 1 ] * x) + [ 2 ]");
    // fHG->SetParameters(1e+3, 5e+11, 1e+3);
    // fHG->SetParNames("A", "k_{2}", "Const");
    // histH->Fit("fH", "", "", Lambda_min, Lambda_max);
    // // TF1 *fOA = new TF1("fA", "5000 * cos([ 0 ] / x) + 5000");
    // // fOA->SetParameter(0, 3e-8);
    // // fOA->SetParName(0, "k_{1}");
    // // histO->Fit("fA", "", "", Lambda_min, Lambda_max);
    // TF1 *f = new TF1("f", "500 * cos([ 0 ] / x + [ 1 ] * x) + 500");
    // f->SetParameters(3e-8, 7e+13);
    // f->SetParNames("k_{1}", "k_{2}");
    // histH->Fit("fA", "", "", Lambda_min, Lambda_max);

    // draw histogram
    TCanvas *c = new TCanvas("c", "O (blue) & H (red) beam count", 700, 500);
    // histO->Rebin(5);
    // histH->Rebin(5);
    histO->SetTitle("O (blue) & H (red) beam count");
    histH->SetTitle("O (blue) & H (red) beam count");
    histH->SetLineColor(2);
    histO->Draw("c");
    histH->Draw("SAME");
    c->Print(OUTPUT);

    TCanvas *c1 = new TCanvas("c1", "O (blue) & H (red) beam count", 700, 500);
    // histO_zoom->Rebin(5);
    // histH_zoom->Rebin(5);
    histO_zoom->SetTitle("O (blue) & H (red) beam count");
    histH_zoom->SetTitle("O (blue) & H (red) beam count");
    histH_zoom->SetLineColor(2);
    histO_zoom->Draw("c1");
    histH_zoom->Draw("SAME");
    c1->Print(OUTPUT_ZOOM);

    inputFileO.close();
    inputFileH.close();

    return 0;
}
