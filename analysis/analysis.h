#include "config.h"

// from count2d.cpp
int count2d(
    string input_tree,
    string output_img_file)
{
    // open the datafile
    TFile *file = TFile::Open(input_tree.c_str());
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file" << std::endl;
        return 1;
    }

    TTree *T;
    file->GetObject("T", T);
    if (!T)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }

    TCanvas *c1 = new TCanvas("c1", "Beam count at xy position", 600, 400);
    T->Draw("y:x", "f==4&&x>=0&&x<=1&&y>=00&&y<=1", "colz");

    // The following is only in order to name the graph.
    TH1F *histogram = (TH1F *)gPad->GetPrimitive("htemp");
    if (histogram)
        histogram->SetTitle("Beam count at xy position");
    else
    {
        std::cerr << "There is no histogram" << std::endl;
    }

    c1->Print(output_img_file.c_str());

    return 0;
}

// from lambda_roi_hist.cpp
int lambda_roi_hist(
    double time_sec,
    double x_min_H,
    double x_max_H,
    double y_min_H,
    double y_max_H,
    double x_min_O,
    double x_max_O,
    double y_min_O,
    double y_max_O,
    string input_tree,
    string output_tree,
    string output_img_tof,
    string output_img_lmd)
{
    Double_t lambda;
    Int_t channel;
    double tof;
    int count;

    // open the datafiles
    TFile *input = TFile::Open(input_tree.c_str());
    if (!input || input->IsZombie())
    {
        cerr << "Failed to open the input." << endl;
        return 1;
    }
    TTree *T;
    input->GetObject("T", T);
    if (!T)
    {
        cerr << "Failed to retrieve the tree from the input." << endl;
        return 1;
    }

    // output file and tree
    TFile output(output_tree.c_str(), "recreate");
    TTree tree("tree", "tree");
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lambda);

    // canvases and regions
    TCanvas *ctof = new TCanvas("ctof", "count on tof H", 600, 400);
    TCanvas *ctof_perSec = new TCanvas("ctof_perSec", "count on tof H per second", 600, 400);
    TCanvas *clmd = new TCanvas("clmd", "count on wavelength", 600, 400);
    TCut xCutH = ("x>=" + to_string(x_min_H) + "&&x<=" + to_string(x_max_H)).c_str();
    TCut yCutH = ("y>=" + to_string(y_min_H) + "&&y<=" + to_string(y_max_H)).c_str();
    TCut xCutO = ("x>=" + to_string(x_min_O) + "&&x<=" + to_string(x_max_O)).c_str();
    TCut yCutO = ("y>=" + to_string(y_min_O) + "&&y<=" + to_string(y_max_O)).c_str();

    // CALCULATION ON CHANNEL H
    // tof-count draw
    channel = H_TDC;
    TH1D *histHtof = new TH1D("histHtof", "tof_count (H: blue, O: red); tof (µs); count", 500, 0, 50000);
    TH1D *histHtof_perSec = new TH1D("histHtof_perSec", "tof count per second (H: blue, O: red); tof (µs); count (/s)", 500, 0, 50000);
    ctof->cd();
    T->Draw("tof>>histHtof", "f==4" && xCutH && yCutH);
    ctof_perSec->cd();
    T->Draw("tof>>histHtof_perSec", "f==4" && xCutH && yCutH, "same eh");
    // divide by time
    for (int bin = 1; bin <= histHtof_perSec->GetXaxis()->GetNbins(); bin++)
    {
        histHtof_perSec->SetBinContent(bin, histHtof_perSec->GetBinContent(bin) / time_sec);
        histHtof_perSec->SetBinError(bin, histHtof_perSec->GetBinError(bin) / time_sec);
    }
    // filling tree
    for (int i = 0; i < histHtof->GetXaxis()->GetNbins(); i++)
    {
        count = histHtof->GetBinContent(i);
        tof = histHtof->GetXaxis()->GetBinCenter(i);
        lambda = tof * 1e-6 * factor2tof; // tof is written in µs
        for (int j = 0; j < count; j++)
        {
            tree.Fill();
        }
    }

    // CALCULATION ON CHANNEL O
    channel = O_TDC;
    TH1D *histOtof = new TH1D("histOtof", "tof_count (H: blue, O: red); tof (µs); count", 500, 0, 40000);
    TH1D *histOtof_perSec = new TH1D("histOtof_perSec", "tof count per second (H: blue, O: red); tof (µs); count (/s)", 500, 0, 40000);
    ctof->cd();
    T->Draw("tof>>histOtof", "f==4" && xCutO && yCutO);
    ctof_perSec->cd();
    T->Draw("tof>>histOtof_perSec", "f==4" && xCutO && yCutO);
    // divide by time
    for (int bin = 1; bin <= histOtof_perSec->GetXaxis()->GetNbins(); bin++)
    {
        histOtof_perSec->SetBinContent(bin, histOtof_perSec->GetBinContent(bin) / time_sec);
        histOtof_perSec->SetBinError(bin, histOtof_perSec->GetBinError(bin) / time_sec);
    }
    // filling tree
    for (int i = 0; i < histOtof->GetXaxis()->GetNbins(); i++)
    {
        count = histOtof->GetBinContent(i);
        tof = histOtof->GetXaxis()->GetBinCenter(i);
        lambda = tof * 1e-6 * factor2tof; // tof is written in µs
        for (int j = 0; j < count; j++)
        {
            tree.Fill();
        }
    }

    // LAMBDA COUNT
    clmd->cd();
    TH1D *histHlmd = new TH1D("histHlmd", "lambda-count (H: blue, O: red): lambda [m]; count [/s]", 500, 0, lambda_max_displayed);
    TH1D *histOlmd = new TH1D("histOlmd", "lambda-count (H: blue, O: red): lambda [m]; count [/s]", 500, 0, lambda_max_displayed);
    histOlmd->SetLineColor(2);
    tree.Draw("lambda>>histHlmd", "channel==0", "eh");
    tree.Draw("lambda>>histOlmd", "channel==0", "eh same");

    ctof_perSec->Print(output_img_tof.c_str());
    clmd->Print(output_img_lmd.c_str());

    // ENDING
    tree.Write();
    input->Close();
    output.Close();
    return 0;
}

// from simulation/motecarlo.h
int read_lambda_perSec(
    double time_sec,
    string input_root,
    string output_img,
    string output_zoom_img,
    string histTitle = "NoTitle")
{
    // open the datafile
    TFile *file = TFile::Open(input_root.c_str());
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

    // build canvases
    TCanvas *c1 = new TCanvas("c1", "histograms on lambda", 600, 400);
    TCanvas *c2 = new TCanvas("c2", "histograms on lambda", 600, 400);
    TCut CutH = ("channel==" + std::to_string(H_TDC)).c_str();
    TCut CutO = ("channel==" + std::to_string(O_TDC)).c_str();
    TString title = Form("O (blue) & H (red) %s; lambda [m]; beam count [/s]", histTitle.c_str());

    // normal histogram
    c1->cd();
    TH1D *histH = new TH1D("histH", "lambda-count of H beam", nBins, lambda_min, lambda_max);
    TH1D *histO = new TH1D("histO", "lambda-count of O beam", nBins, lambda_min, lambda_max);
    histO->SetTitle(title);
    histH->SetTitle(title);
    histH->SetLineColor(2);
    tree->Draw("lambda>>histH", CutH, "eh");
    tree->Draw("lambda>>histO", CutO, "same eh");

    c2->cd();
    TH1D *histH_zoom = new TH1D("histH_zoom", "lambda-count of H beam", nBins_zoom, lambda_min_focus, lambda_max_focus);
    TH1D *histO_zoom = new TH1D("histO_zoom", "lambda-count of O beam", nBins_zoom, lambda_min_focus, lambda_max_focus);
    histO_zoom->SetTitle(title);
    histH_zoom->SetTitle(title);
    histH_zoom->SetLineColor(2);
    tree->Draw("lambda>>histH_zoom", CutH, "eh");
    tree->Draw("lambda>>histO_zoom", CutO, "same eh");

    // divide by time
    for (int bin = 1; bin <= histH->GetXaxis()->GetNbins(); bin++)
    {
        histH->SetBinContent(bin, histH->GetBinContent(bin) / time_sec);
        histH->SetBinError(bin, histH->GetBinError(bin) / time_sec);
    }
    for (int bin = 1; bin <= histO->GetXaxis()->GetNbins(); bin++)
    {
        histO->SetBinContent(bin, histO->GetBinContent(bin) / time_sec);
        histO->SetBinError(bin, histO->GetBinError(bin) / time_sec);
    }
    for (int bin = 1; bin <= histH_zoom->GetXaxis()->GetNbins(); bin++)
    {
        histH_zoom->SetBinContent(bin, histH_zoom->GetBinContent(bin) / time_sec);
        histH_zoom->SetBinError(bin, histH_zoom->GetBinError(bin) / time_sec);
    }
    for (int bin = 1; bin <= histO_zoom->GetXaxis()->GetNbins(); bin++)
    {
        histO_zoom->SetBinContent(bin, histO_zoom->GetBinContent(bin) / time_sec);
        histO_zoom->SetBinError(bin, histO_zoom->GetBinError(bin) / time_sec);
    }

    c1->Print(output_img.c_str());
    c2->Print(output_zoom_img.c_str());

    return 0;
}

int oscillation_lambda(
    string input_root,
    string output_img_normal,
    string output_img_zoom,
    string output_img_fourier,
    string output_img_fourier_wide,
    double delta_deg = 30,
    double theta_deg = 1.,
    string histTitle = "NoTitle")
{
    // open the data file
    TFile *file = TFile::Open(input_root.c_str());
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
    TH1F *histO = new TH1F("histO", ";(I_H-I_O)/(I_H+I_O);lambda [m]", nBins, lambda_min, lambda_max);
    TH1F *histH = new TH1F("histH", ";(I_H-I_O)/(I_H+I_O);lambda [m]", nBins, lambda_min, lambda_max);
    TH1F *histH_zoom = new TH1F("histH_zoom", ";(I_H-I_O)/(I_H+I_O);lambda [m]", nBins, lambda_min_focus, lambda_max_focus);
    TH1F *histO_zoom = new TH1F("histO_zoom", ";(I_H-I_O)/(I_H+I_O);lambda [m]", nBins, lambda_min_focus, lambda_max_focus);
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
    int nBinsFourier = hFourier_n->GetXaxis()->GetNbins();
    double k_max = nBinsFourier * 2 * TMath::Pi() / (lambda_max_focus - lambda_min_focus);
    TH1D *hFourier = new TH1D("hFourierK", "Fourier transformation of (I_H-I_O)/(I_H+I_O); wave number [/m]; amp", nBinsFourier, 0, k_max);
    for (int bin = 0; bin < nBinsFourier; bin++)
    {
        hFourier->SetBinContent(bin, hFourier_n->GetBinContent(bin));
    }

    // fittings
    double peak_position = 2. * M_PI * g * m * m / (h * h) * 2. * gap * mirror_distance / tan(2. * theta_deg * M_PI / 180) * sin(delta_deg * M_PI / 180);
    // oscillation
    const double fit_range_oscil[] = {6.9, 8.5}; // TODO: find a good value
    TF1 *f_oscil = new TF1("oscil", "[0] * cos(x * [1] + x / [2] + [3]) + [4] * (x - [5])^2 + [6]", fit_range_oscil[0], fit_range_oscil[1]);
    f_oscil->SetParameters(0.5, peak_position, 1e10, 0, 1e-10, 7.5e-10, 0);
    h_oscil->Fit("oscil", "", "", fit_range_oscil[0], fit_range_oscil[1]);
    Double_t k1oscil = f_oscil->GetParameter(1);
    Double_t k2oscil = f_oscil->GetParameter(2);
    // peak
    const double fit_range_peak[] = {peak_position - fit_width_peak, peak_position + fit_width_peak};
    TF1 *f_peak = new TF1("gaus", "[0] * exp(-0.5 * ((x-[1]) / [2])^2)", fit_range_peak[0], fit_range_peak[1]);
    f_peak->SetParameters(fit_peak_height, peak_position, fit_width_peak);
    hFourier->Fit("gaus", "", "", fit_range_peak[0], fit_range_peak[1]);
    Double_t p0_peak = f_peak->GetParameter(1);
    Double_t p0e_peak = f_peak->GetParError(1);
    Double_t chi2 = f_peak->GetChisquare();
    Int_t Ndof = f_peak->GetNDF();

    TH1D *hFourier_wide = (TH1D *)hFourier->Clone();
    // hFourier->GetXaxis()->SetRangeUser(0, 1e12);

    // draw histogram
    TCanvas *c1 = new TCanvas("c1", "oscillation", 700, 500);
    TCanvas *c2 = new TCanvas("c2", "oscillation zoom", 700, 500);
    TCanvas *c3 = new TCanvas("c3", "Fourier", 700, 500);
    TCanvas *c4 = new TCanvas("c4", "Fourier wide", 700, 500);
    c1->cd();
    h_oscil->SetTitle(Form("%s; #lambda [m]; (I_H - I_O) / (I_H + I_O)", histTitle.c_str()));
    h_oscil->GetYaxis()->SetRangeUser(-2, 2);
    h_oscil->Draw("eh");
    c1->Print(output_img_normal.c_str());
    c2->cd();
    h_oscil_zoom->SetTitle(Form("%s; #lambda [m]; (I_H - I_O) / (I_H + I_O)", histTitle.c_str()));
    h_oscil->GetYaxis()->SetRangeUser(-2, 2);
    h_oscil_zoom->Draw("eh");
    c2->Print(output_img_zoom.c_str());
    c3->cd();
    hFourier->SetTitle(Form("%s; wave number [/m]; Fourier amp", histTitle.c_str()));
    hFourier->Draw("eh");
    c3->Print(output_img_fourier.c_str());
    c4->cd();
    hFourier_wide->SetTitle(Form("%s; wave number [/m]; Fourier amp", histTitle.c_str()));
    hFourier_wide->Draw("eh");
    c4->Print(output_img_fourier_wide.c_str());

    // result of peak position
    double resultK = -p0_peak * 2 * TMath::Pi() / (lambda_max_focus - lambda_min_focus);
    double resultK_error = p0e_peak * 2 * TMath::Pi() / (lambda_max_focus - lambda_min_focus);
    std::cout << "k = " << -p0_peak << " +- " << p0e_peak << std::endl;

    return 0;
}
