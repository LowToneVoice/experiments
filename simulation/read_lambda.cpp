#include "simulate.h"

int read_lambda_perSec_withCut(
    string phase_contrib_input = PHASE_CONTRIB,
    string main_sub_input = MAIN_SUB,
    string time_min_input = TIME_MIN,
    string angle_delta_deg_input = ANGLE_DELTA_DEG,
    string angle_from_parallel_deg_input = ANGLE_FROM_PARALLEL_DEG,
    string lmd_used_min_input = LMD_USED_MIN,
    string lmd_used_max_input = LMD_USED_MAX,
    string file_extension = FILE_EXTENSION)
{
    // data file names
    string file_format = FILE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);
    string INPUT_TREE = "./dat/montecarlo/root/" + file_format + ".root";
    string OUTPUT = "./beam_count/montecarlo/normal/" + file_format + "." + file_extension;
    string OUTPUT_ZOOM = "./beam_count/montecarlo/zoom/" + file_format + "." + file_extension;
    string title_format = TITLE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);

    // input values
    const double lambda_min_used = stod(lmd_used_min_input);
    const double lambda_max_used = stod(lmd_used_max_input);
    const int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);
    const double time_sec = stod(time_min_input) * 60;

    // open the datafile
    TFile *file = TFile::Open(INPUT_TREE.c_str());
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }
    TTree *tree;
    file->GetObject("tree", tree);
    if (!tree)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }

    TCanvas *c1 = new TCanvas("c1", "histograms on lambda", 600, 400);
    TCanvas *c2 = new TCanvas("c2", "histograms on lambda", 600, 400);
    TCut CutH = ("channel==" + std::to_string(H_TDC)).c_str();
    TCut CutO = ("channel==" + std::to_string(O_TDC)).c_str();
    TString title = Form("O (blue) & H (red) %s; lambda [m]; beam count", title_format.c_str());

    c1->cd();
    TH1D *histH = new TH1D("histH", "lambda-count of H beam", N_loop_lambda, lambda_min, lambda_max);
    TH1D *histO = new TH1D("histO", "lambda-count of O beam", N_loop_lambda, lambda_min, lambda_max);

    // divide by time
    for (int bin = 0; bin < histH->GetXaxis()->GetNbins(); bin++)
    {
        histH->SetBinContent(bin, histH->GetBinContent(bin) / time_sec);
    }
    for (int bin = 0; bin < histO->GetXaxis()->GetNbins(); bin++)
    {
        histO->SetBinContent(bin, histO->GetBinContent(bin) / time_sec);
    }

    histO->SetTitle(title);
    histH->SetTitle(title);
    histH->SetLineColor(2);
    tree->Draw("lambda>>histH", CutH);
    tree->Draw("lambda>>histO", CutO, "same");
    c1->Print(OUTPUT.c_str());

    c2->cd();
    TH1D *histH_zoom = new TH1D("histH_zoom", "lambda-count of H beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    TH1D *histO_zoom = new TH1D("histO_zoom", "lambda-count of O beam", N_loop_lambda, lambda_min_used, lambda_max_used);
    histO_zoom->SetTitle(title);
    histH_zoom->SetTitle(title);
    histH_zoom->SetLineColor(2);
    tree->Draw("lambda>>histH_zoom", CutH);
    tree->Draw("lambda>>histO_zoom", CutO, "same");
    c2->Print(OUTPUT_ZOOM.c_str());

    return 0;
}
