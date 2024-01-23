#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>

// data labels
std::string PHASE_CONTRIB = "g";
std::string MAIN_SUB = "mix";
std::string ANGLE_DELTA_DEG = "90";
std::string TIME_MIN = "20";
std::string ANGLE_FROM_PARALLEL_DEG = "1e-2";

// fitting range and initial conditions
constexpr double fit_width = 2e11;
// constexpr double fit_range[] = {7, 12};
constexpr double fit_par_height = 30;
// constexpr double fit_par_position = 9;
// constexpr double fit_par_width = 1;

// datafiles
std::string FORMAT = PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG;
std::string INPUT_TREE = "./dat/montecarlo/root/" + FORMAT + ".root";
std::string OUTPUT_FILE = "./oscil_graph/montecarlo/normal/" + FORMAT + ".pdf";
std::string OUTPUT_FILE_ZOOM = "./oscil_graph/montecarlo/zoom/" + FORMAT + ".pdf";
std::string OUTPUT_FILE_FOURIER = "./oscil_graph/montecarlo/fourier/" + FORMAT + ".pdf";
std::string OUTPUT_FILE_FOURIER_WIDE = "./oscil_graph/montecarlo/fourier_wide/" + FORMAT + ".pdf";

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
constexpr double gap = 189e-6;

// alignment variables
constexpr double lambda_min_used = 6.9e-10;
constexpr double lambda_max_used = 8.4e-10;
constexpr double total_length = 1.;
constexpr double theta = 1.05 * pi / 180;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;
constexpr double mirror_distance = 150e-3;

constexpr int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);

int oscillation_lambda()
{
    // open the data file
    TFile *file = TFile::Open(INPUT_TREE.c_str());
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
    TH1F *histO = new TH1F("histO", ";(I_H-I_O)/(I_H+I_O);lambda [m]", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH = new TH1F("histH", ";(I_H-I_O)/(I_H+I_O);lambda [m]", N_loop_lambda, lambda_min, lambda_max);
    TH1F *histH_zoom = new TH1F("histH_zoom", ";(I_H-I_O)/(I_H+I_O);lambda [m]", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1F *histO_zoom = new TH1F("histO_zoom", ";(I_H-I_O)/(I_H+I_O);lambda [m]", N_loop_lambda, lambda_min_used, lambda_max_used);
    tree->Draw("lambda>>histH", CutH);
    tree->Draw("lambda>>histO", CutO);
    tree->Draw("lambda>>histH_zoom", CutH);
    tree->Draw("lambda>>histO_zoom", CutO);

    double lambda;

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

    // transform FFT from n space to K space
    TH1 *hFourier_n = h1_zoom->FFT(0, "MAG");
    int Nbins = hFourier_n->GetXaxis()->GetNbins();
    double k_max = Nbins * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
    TH1D *hFourier = new TH1D("hFourierK", "Fourier transformation of (I_H-I_O)/(I_H+I_O); wave number [/m]; amp", Nbins, 0, k_max);
    for (int bin = 0; bin < Nbins; bin++)
    {
        hFourier->SetBinContent(bin, hFourier_n->GetBinContent(bin));
    }

    // fittings
    const double peak_position = 2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * sin(std::stod(ANGLE_DELTA_DEG) * pi / 180);
    const double fit_range[] = {peak_position - fit_width, peak_position + fit_width};
    TF1 *f = new TF1("gaus", "[0] * exp(-0.5 * ((x-[1]) / [2])^2)", fit_range[0], fit_range[1]);
    f->SetParameters(fit_par_height, peak_position, fit_width);
    hFourier->Fit("gaus", "", "", fit_range[0], fit_range[1]);
    Double_t p0 = f->GetParameter(1);
    Double_t p0e = f->GetParError(1);
    Double_t chi2 = f->GetChisquare();
    Int_t Ndof = f->GetNDF();

    TH1D *hFourier_wide = (TH1D *)hFourier->Clone();
    hFourier->GetXaxis()->SetRangeUser(0, 1e12);

    // draw histogram
    TCanvas *c1 = new TCanvas("c1", "oscillation", 700, 500);
    TCanvas *c2 = new TCanvas("c2", "oscillation zoom", 700, 500);
    TCanvas *c3 = new TCanvas("c3", "Fourier", 700, 500);
    TCanvas *c4 = new TCanvas("c4", "Fourier wide", 700, 500);
    c1->cd();
    h1->Draw();
    c1->Print(OUTPUT_FILE.c_str());
    c2->cd();
    h1_zoom->Draw();
    c2->Print(OUTPUT_FILE_ZOOM.c_str());
    c3->cd();
    hFourier->Draw();
    c3->Print(OUTPUT_FILE_FOURIER.c_str());
    c4->cd();
    hFourier_wide->Draw();
    c4->Print(OUTPUT_FILE_FOURIER_WIDE.c_str());

    // result of peak position
    double resultK = -p0 * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
    double resultK_error = p0e * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
    std::cout << "k = " << -p0 << " +- " << p0e << std::endl;

    return 0;
}
