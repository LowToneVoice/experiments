#include "simulate.h"

/*
    FROM read_lambda.cpp
*/

int read_lambda(
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


/*
    FROM oscillation_lambda.cpp
*/

int oscillation_lambda(
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
    string OUTPUT_FILE = "./oscil_graph/montecarlo/normal/" + file_format + "." + file_extension;
    string OUTPUT_FILE_ZOOM = "./oscil_graph/montecarlo/zoom/" + file_format + "." + file_extension;
    string OUTPUT_FILE_FOURIER = "./oscil_graph/montecarlo/fourier/" + file_format + "." + file_extension;
    string OUTPUT_FILE_FOURIER_WIDE = "./oscil_graph/montecarlo/fourier_wide/" + file_format + "." + file_extension;
    string title_format = TITLE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);

    // input values
    const double lambda_min_used = stod(lmd_used_min_input);
    const double lambda_max_used = stod(lmd_used_max_input);
    const double time_min = stod(time_min_input);
    const double angle_delta_deg = stod(angle_delta_deg_input);
    const double angle_from_parallel_deg = stod(angle_from_parallel_deg_input);

    const int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);

    // open the data file
    TFile *file = TFile::Open(INPUT_TREE.c_str());
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }
    TTree *tree;
    file->GetObject("tree", tree);
    if (!tree)
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
    double count[] = {0, 0};

    // data calculations
    TH1F *h_oscil = (TH1F *)histO->Clone();
    TH1F *h_oscil_zoom = (TH1F *)histO_zoom->Clone();
    for (int bin = 0; bin < histO->GetXaxis()->GetNbins() && bin < histH->GetXaxis()->GetNbins(); bin++)
    {
        count[0] = histO->GetBinContent(bin);
        count[1] = histH->GetBinContent(bin);
        if (count[0] + count[1] == 0)
        {
            std::cerr << "Warning: Division by zero in the loop. bin = " << bin << "of normal" << std::endl;
            h_oscil->SetBinContent(bin, 0);
            h_oscil->SetBinError(bin, 100);
        }
        else
        {
            h_oscil->SetBinContent(bin, (count[0] - count[1]) / (count[0] + count[1]));
            h_oscil->SetBinError(bin, 2 * sqrt(count[0] * count[1] / pow(count[0] + count[1], 3)));
        }
    }
    for (int bin = 0; bin < histO_zoom->GetXaxis()->GetNbins() && bin < histH_zoom->GetXaxis()->GetNbins(); bin++)
    {
        count[0] = histO_zoom->GetBinContent(bin);
        count[1] = histH_zoom->GetBinContent(bin);
        if (count[0] + count[1] == 0)
        {
            std::cerr << "Warning: Division by zero in the loop. bin = " << bin << "of zoom" << std::endl;
        }
        else
        {
            h_oscil->SetBinContent(bin, (count[0] - count[1]) / (count[0] + count[1]));
            h_oscil->SetBinError(bin, 2 * sqrt(count[0] * count[1] / pow(count[0] + count[1], 3)));
        }
    }

    // transform FFT from n space to K space
    TH1 *hFourier_n = h_oscil_zoom->FFT(0, "MAG");
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
    h_oscil->SetTitle(Form("%s; #lambda [m]; (I_H - I_O) / (I_H + I_O)", title_format.c_str()));
    h_oscil->GetYaxis()->SetRangeUser(-2, 2);
    h_oscil->Draw("eh");
    c1->Print(OUTPUT_FILE.c_str());
    c2->cd();
    h_oscil_zoom->SetTitle(Form("%s; #lambda [m]; (I_H - I_O) / (I_H + I_O)", title_format.c_str()));
    h_oscil->GetYaxis()->SetRangeUser(-2, 2);
    h_oscil_zoom->Draw("eh");
    c2->Print(OUTPUT_FILE_ZOOM.c_str());
    c3->cd();
    hFourier->SetTitle(Form("%s; wave number [/m]; Fourier amp", title_format.c_str()));
    hFourier->Draw("eh");
    c3->Print(OUTPUT_FILE_FOURIER.c_str());
    c4->cd();
    hFourier_wide->SetTitle(Form("%s; wave number [/m]; Fourier amp", title_format.c_str()));
    hFourier_wide->Draw("eh");
    c4->Print(OUTPUT_FILE_FOURIER_WIDE.c_str());

    // result of peak position
    double resultK = -p0 * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
    double resultK_error = p0e * 2 * TMath::Pi() / (lambda_max_used - lambda_min_used);
    std::cout << "k = " << -p0 << " +- " << p0e << std::endl;

    return 0;
}
