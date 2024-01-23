#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>

using namespace std;

// data labels
string PHASE_CONTRIB = "mix";
string MAIN_SUB = "mix";
string ANGLE_DELTA_DEG = "30";
string TIME_MIN = "60";
string ANGLE_FROM_PARALLEL_DEG = "1e-1";
string LMD_USED_MIN = "7e-10";
string LMD_USED_MAX = "10e-10";

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
constexpr double total_length = 1.;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
constexpr double d_lambda = h / total_length / m * dt;

int read_lambda(
    string phase_contrib_input = PHASE_CONTRIB,
    string main_sub_input = MAIN_SUB,
    string time_min_input = TIME_MIN,
    string angle_delta_deg_input = ANGLE_DELTA_DEG,
    string angle_from_parallel_deg_input = ANGLE_FROM_PARALLEL_DEG,
    string lmd_used_min_input = LMD_USED_MIN,
    string lmd_used_max_input = LMD_USED_MAX
)
{
    // data file names
    string FILE_FORMAT = phase_contrib_input + "_" + main_sub_input + "_" + angle_delta_deg_input + "deg_" + time_min_input + "min_ALPHA" + angle_from_parallel_deg_input + "_lmd" + lmd_used_min_input + "to" + lmd_used_max_input;
    string INPUT_TREE = "./dat/montecarlo/root/" + FILE_FORMAT + ".root";
    string OUTPUT = "./beam_count/montecarlo/normal/" + FILE_FORMAT + ".pdf";
    string OUTPUT_ZOOM = "./beam_count/montecarlo/zoom/" + FILE_FORMAT + ".pdf";

    // input values
    const double lambda_min_used = stod(lmd_used_min_input);
    const double lambda_max_used = stod(lmd_used_max_input);
    const int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);

    // open the datafile
    TFile *file = TFile::Open(INPUT_TREE.c_str());
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
    histO->SetTitle("O (blue) & H (red); lambda [m]; beam count");
    histH->SetTitle("O (blue) & H (red) beam count; lambda [m]; beam count");
    histH->SetLineColor(2);
    tree->Draw("lambda>>histH", CutH);
    tree->Draw("lambda>>histO", CutO, "same");
    c1->Print(OUTPUT.c_str());

    c2->cd();
    TH1D *histH_zoom = new TH1D("histH_zoom", "lambda-count of H beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1D *histO_zoom = new TH1D("histO_zoom", "lambda-count of O beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    histO_zoom->SetTitle("O (blue) & H (red); lambda [m]; beam count");
    histH_zoom->SetTitle("O (blue) & H (red); lambda [m]; beam count");
    histH_zoom->SetLineColor(2);
    tree->Draw("lambda>>histH_zoom", CutH);
    tree->Draw("lambda>>histO_zoom", CutO, "same");
    c2->Print(OUTPUT_ZOOM.c_str());

    return 0;
}
