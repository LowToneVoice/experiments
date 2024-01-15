#include<iostream>
#include<TFile.h>
#include<TH1.h>
#include<TH1F.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TAxis.h>
#include<TArrayD.h>
#include<TROOT.h>

// datafiles
#define INPUT_TREE "./BL05/20231224225335_list_copy.root"
#define OUTPUT_TREE "./dat/experiments/20231224225335.root"
#define OUTPUT_FILE "./beam_count/test/20231224225335.pdf"
#define OUTPUT_FILE_TOF "./beam_count/test/20231224225335_tof.pdf"

#define h 6.62607015e-34
#define m 1.6749e-27

#define H_CHANNEL 0
#define O_CHANNEL 1

constexpr double x_range_H[] = {.41, .44};
constexpr double y_range_H[] = {.34, .48};
constexpr double x_range_O[] = {.52, .54};
constexpr double y_range_O[] = {.32, .48};

constexpr double total_length = 1.;
constexpr double d_total_length = .005;

constexpr double factor_tof2lambda = h / total_length / m;
constexpr double factor_tof2lambda_max = h / (total_length - d_total_length) / m;

int lambda_roi_histogram()
{
    Double_t lambda;
    Double_t lambdaUpperLim;
    Int_t channel;
    double tof;
    int count;

    // open the datafile
    TFile *input = TFile::Open(INPUT_TREE);
    if (!input || input->IsZombie())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }
    TTree *T;
    input->GetObject("T", T);
    if(!T)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }

    // output file and tree
    TFile output(OUTPUT_TREE, "recreate");
    TTree tree("tree", "tree");
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lambda);
    tree.Branch("lambdaUpperLim", &lambdaUpperLim);

    TCanvas *c1 = new TCanvas("c1", "count on tof of H channel", 600, 400);
    TCanvas *c2 = new TCanvas("c2", "count on tof of O channel", 600, 400);
    TCanvas *c3 = new TCanvas("c3", "count on lambda of both channel", 600, 400);

    // trans of channel H
    channel = H_CHANNEL;
    c1->cd();
    TCut xCutH = ("x>=" + std::to_string(x_range_H[0]) + "&&x<=" + std::to_string(x_range_H[1])).c_str();
    TCut yCutH = ("y>=" + std::to_string(y_range_H[0]) + "&&y<=" + std::to_string(y_range_H[1])).c_str();
    T->Draw("tof>>h_tof(500,0,40000)", "f==4"&&xCutH&&yCutH);
    TH1D *histH = (TH1D *)gROOT->FindObject("h_tof");
    for (int i = 0; i < histH->GetXaxis()->GetNbins(); i++)
    {
        count = histH->GetBinContent(i);
        tof = histH->GetXaxis()->GetBinCenter(i);
        lambda = tof * factor_tof2lambda * 1e-6;    // tof is written in µs
        lambdaUpperLim = tof * factor_tof2lambda_max * 1e-6;
        for (int j = 0; j < count; j++)
            tree.Fill();
    }
    c1->Print(OUTPUT_FILE_TOF);

    // trans of channel O
    channel = O_CHANNEL;
    c2->cd();
    TCut xCutO = ("x>=" + std::to_string(x_range_O[0]) + "&&x<=" + std::to_string(x_range_O[1])).c_str();
    TCut yCutO = ("y>=" + std::to_string(y_range_O[0]) + "&&y<=" + std::to_string(y_range_O[1])).c_str();
    T->Draw("tof>>o_tof(500,0,40000)", "f==4" && xCutO && yCutO);
    TH1D *histO = (TH1D *)gROOT->FindObject("o_tof");
    std::cout << histO->GetXaxis()->GetBinCenter(50) << std::endl;
    for (int i = 0; i < histO->GetXaxis()->GetNbins(); i++)
    {
        count = histO->GetBinContent(i);
        tof = histO->GetXaxis()->GetBinCenter(i);
        lambda = tof * factor_tof2lambda * 1e-6;
        lambdaUpperLim = tof * lambdaUpperLim * 1e-6;
        for (int j = 0; j < count; j++)
            tree.Fill();
    }

    c3->cd();
    tree.Draw("lambda>>h_lambda(500,0,1e-8)", "channel==0", "");
    tree.Draw("lambda>>o_lambda(500,0,1e-8)", "channel==1", "SAME");
    TH1D *hist = (TH1D *)gROOT->FindObject("o_lambda");
    hist->SetLineColor(kRed);
    c3->Print(OUTPUT_FILE);

    tree.Write();
    input->Close();
    output.Close();

    return 0;
}

int main()
{
    return lambda_roi_histogram();
}
