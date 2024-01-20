#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>

// data labels
std::string PHASE_CONTRIB = "mix";
std::string MAIN_SUB = "mix";
std::string ANGLE_DELTA_DEG = "30";
std::string TIME_MIN = "10";
std::string ANGLE_FROM_PARALLEL_DEG = "2e-3";

// fitting range and initial conditions
constexpr double fit_range[] = {7, 12};
constexpr double fit_par_height = 300;
constexpr double fit_par_position = 9;
constexpr double fit_par_width = 1;

// datafiles
std::string INPUT_TREE = "./dat/montecarlo/ref/lambda/" + PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + ".root";
std::string OUTPUT_FILE = "./oscil_graph/montecarlo/ref/lambda/" + PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + ".pdf";
std::string OUTPUT_FILE_ZOOM = "./oscil_graph/montecarlo/ref/lambda/zoom/" + PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + ".pdf";
std::string OUTPUT_FILE_FOURIER = "./oscil_graph/montecarlo/ref/lambda/" + PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + "_fourier.pdf";
std::string OUTPUT_FILE_FOURIER_WIDE = "./oscil_graph/montecarlo/ref/lambda/" + PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + "_fourier_wide.pdf";

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
constexpr double lambda_min_used = 6.9e-10;
constexpr double lambda_max_used = 8.4e-10;
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
    // open the data file
    TFile *file = TFile::Open(INPUT_TREE);
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }
    TTree *tree;
    file->GetObject("tree", tree);
    if(!tree)
    {
        std::cerr << "There is no tree available." << std::endl;
        return 1;
    }

    // conditions
    TCut CutH = ("channel==" + std::to_string(H_TDC)).c_str();
    TCut CutO = ("channel==" + std::to_string(O_TDC)).c_str();

    // make histogram
    TH1F *histO = new TH1F("histO", "Oscil", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH = new TH1F("histH", "Oscil", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH_zoom = new TH1F("histH_zoom", "Oscil", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1F *histO_zoom = new TH1F("histO_zoom", "Oscil", N_loop_lambda, lambda_min_used, lambda_max_used);
    tree->Draw("lambda>>histH", CutH);
    tree->Draw("lambda>>histO", CutO);
    tree->Draw("lambda>>histH_zoom", CutH);
    tree->Draw("lambda>>histO_zoom", CutO);

    double lambda, theta;

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
    TH1 *hFourier_wide = h1_zoom->FFT(0, "MAG");
    TH1 *hFourier = h1_zoom->FFT(0, "MAG");
    hFourier->GetXaxis()->SetRangeUser(0, 50);

    // fittings
    TF1 *f = new TF1("gaus", "[0]*exp(-0.5*((x-[1])/[2])^2)", fit_range[0], fit_range[1]);
    f->SetParameters(fit_par_height, fit_par_position, fit_par_width);
    hFourier->Fit("gaus","","", fit_range[0], fit_range[1]);
    Double_t p0 = f->GetParameter(1);
    Double_t p0e = f->GetParError(1);
    Double_t chi2 = f->GetChisquare();
    Int_t Ndof = f->GetNDF();

    // draw histogram
    TCanvas *c = new TCanvas("c", "oscillation", 700, 500);
    // h1->Fit("f", "", "", 7e-10, 10e-10);
    h1->Draw("c");
    c->Print(OUTPUT_FILE);

    TCanvas *c1 = new TCanvas("c1", "oscillation", 700, 500);
    // h1_zoom->Fit("f", "", "", lambda_min_used, lambda_max_used);
    h1_zoom->Draw("c1");
    c1->Print(OUTPUT_FILE_ZOOM);

    TCanvas *c2 = new TCanvas("c2", "Fourier components", 700, 500);
    hFourier->Draw();
    c2->Print(OUTPUT_FILE_FOURIER);

    delete histO;
    delete histH;
    delete histO_zoom;
    delete histH_zoom;
    if (hFourier)
    {
        TCanvas *c2 = new TCanvas("c2", "Fourier components", 700, 500);
        hFourier->Draw();
        c2->Print(OUTPUT_FILE_FOURIER);
        delete hFourier; // Don't forget to delete the histogram when you're done with it
        double resultK = -p0 * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
        std::cout << "k = " << resultK << std::endl;
    }
    else
    {
        std::cerr << "Failed to perform FFT." << std::endl;
    }


    return 0;
}
