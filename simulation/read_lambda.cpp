#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>

// FILES
#define INPUT_TREE "./dat/montecarlo/ref/lambda/g_mix_90deg_N8e6.root"
#define OUTPUT "./beam_count/montecarlo/ref/lambda/g_mix_90deg_N8e6.pdf"
#define OUTPUT_ZOOM "./beam_count/montecarlo/ref/lambda/g_mix_90deg_N8e6_zoom.pdf"

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

int read_lambda()
{
    // open the datafile
    TFile *file = TFile::Open(INPUT_TREE);
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }
    TTree *tree;
    file->GetObject("tree", tree);
    if(!tree)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }

    TCanvas *c1 = new TCanvas("c1", "histograms on lambda", 600, 400);
    TCanvas *c2 = new TCanvas("c2", "histograms on lambda", 600, 400);
    TCut CutH = ("channel==" + std::to_string(H_TDC)).c_str();
    TCut CutO = ("channel==" + std::to_string(O_TDC)).c_str();

    c1->cd();
    TH1D *histH = new TH1D("histH", "lambda-count of H beam", N_loop_lambda, lambda_min, lambda_max);
    TH1D *histO = new TH1D("histO", "lambda-count of O beam", N_loop_lambda, lambda_min, lambda_max);
    histO->SetTitle("O (blue) & H (red) beam count");
    histH->SetTitle("O (blue) & H (red) beam count");
    histH->SetLineColor(2);
    tree->Draw("lambda>>histH", CutH);
    tree->Draw("lambda>>histO", CutO, "same");
    c1->Print(OUTPUT);

    c2->cd();
    TH1D *histH_zoom = new TH1D("histH_zoom", "lambda-count of H beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1D *histO_zoom = new TH1D("histO_zoom", "lambda-count of O beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    histO_zoom->SetTitle("O (blue) & H (red) beam count");
    histH_zoom->SetTitle("O (blue) & H (red) beam count");
    histH_zoom->SetLineColor(2);
    tree->Draw("lambda>>histH_zoom", CutH);
    tree->Draw("lambda>>histO_zoom", CutO, "same");
    c2->Print(OUTPUT_ZOOM);

    return 0;
}
