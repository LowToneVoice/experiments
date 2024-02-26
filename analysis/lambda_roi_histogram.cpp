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
#define OUTPUT_TREE_UPPER_LIM "./dat/experiments/20231224225335_upperLim.root"
#define OUTPUT_FILE "./beam_count/test/20231224225335.pdf"
#define OUTPUT_FILE_TOF "./beam_count/test/20231224225335_tof.pdf"

constexpr double x_range_H[] = {.41, .44};
constexpr double y_range_H[] = {.34, .48};
constexpr double x_range_O[] = {.52, .54};
constexpr double y_range_O[] = {.32, .48};

#define h 6.62607015e-34
#define m 1.6749e-27

#define H_CHANNEL 0
#define O_CHANNEL 1

constexpr double lambda_max_displayed = 1.5e-8;
constexpr double total_length = 1.;
constexpr double d_total_length = .005;

constexpr double factor_tof2lambda = h / total_length / m;
constexpr double factor_tof2lambda_max = h / (total_length - d_total_length) / m;

int lambda_roi_histogram()
{
    Double_t lambda;
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

    // canvases and regions
    TCanvas *ctof = new TCanvas("ctof", "count on tof of H channel", 600, 400);
    TCanvas *cBothLambda = new TCanvas("cBothLambda", "count on lambda of both channel", 600, 400);
    TCut xCutH = ("x>=" + std::to_string(x_range_H[0]) + "&&x<=" + std::to_string(x_range_H[1])).c_str();
    TCut yCutH = ("y>=" + std::to_string(y_range_H[0]) + "&&y<=" + std::to_string(y_range_H[1])).c_str();
    TCut xCutO = ("x>=" + std::to_string(x_range_O[0]) + "&&x<=" + std::to_string(x_range_O[1])).c_str();
    TCut yCutO = ("y>=" + std::to_string(y_range_O[0]) + "&&y<=" + std::to_string(y_range_O[1])).c_str();


    // NORMAL
    ctof->cd();

    // trans of channel H
    channel = H_CHANNEL;
    TH1D *histHtof = new TH1D("histHtof", "tof-count (H: blue, O: red)", 500, 0, 40000);
    T->Draw("tof>>histHtof", "f==4" && xCutH && yCutH);
    // filling histogram
    for (int i = 0; i < histHtof->GetXaxis()->GetNbins(); i++)
    {
        count = histHtof->GetBinContent(i);
        tof = histHtof->GetXaxis()->GetBinCenter(i);
        lambda = tof * factor_tof2lambda * 1e-6;    // tof is written in µs
        for (int j = 0; j < count; j++)
            tree.Fill();
    }

    // trans of channel O
    channel = O_CHANNEL;
    TH1D *histOtof = new TH1D("histOtof", "tof-count (H: blue, O: red)", 500, 0, 40000);
    histOtof->SetLineColor(2);
    T->Draw("tof>>histOtof", "f==4" && xCutO && yCutO, "same");
    // filling histogram
    for (int i = 0; i < histOtof->GetXaxis()->GetNbins(); i++)
    {
        count = histOtof->GetBinContent(i);
        tof = histOtof->GetXaxis()->GetBinCenter(i);
        lambda = tof * factor_tof2lambda * 1e-6; // tof is written in µs
        for (int j = 0; j < count; j++)
            tree.Fill();
    }

    ctof->Print(OUTPUT_FILE_TOF);

    // both histograms of lambda-count
    cBothLambda->cd();
    TH1D *histH = new TH1D("histH", "lambda-count (H: blue, O: red)", 500, 0, lambda_max_displayed);
    TH1D *histO = new TH1D("histO", "lambda-count (H: blue, O: red)", 500, 0, lambda_max_displayed);
    histO->SetLineColor(2);
    tree.Draw("lambda>>histH", "channel==0", "");
    tree.Draw("lambda>>histO", "channel==1", "SAME");
    cBothLambda->Print(OUTPUT_FILE);

    tree.Write();
    output.Close();


    // // UPPER LIMIT
    // Double_t lambdaUpperLim;
    // Int_t channelUpperLim;
    // TFile outputUpperLim(OUTPUT_TREE_UPPER_LIM, "recreate");
    // TTree treeUpperLim("tree", "tree");
    // treeUpperLim.Branch("channel", &channelUpperLim);
    // treeUpperLim.Branch("lambda", &lambdaUpperLim);

    // // channel H
    // channelUpperLim = H_CHANNEL;
    // for (int i = 0; i < histHtof->GetXaxis()->GetNbins(); i++)
    // {
    //     count = histHtof->GetBinContent(i);
    //     tof = histHtof->GetXaxis()->GetBinCenter(i);
    //     lambdaUpperLim = tof * factor_tof2lambda_max * 1e-6; // tof is written in µs
    //     for (int j = 0; j < count; j++)
    //         treeUpperLim.Fill();
    // }

    // // trans of channel O
    // channelUpperLim = O_CHANNEL;
    // for (int i = 0; i < histOtof->GetXaxis()->GetNbins(); i++)
    // {
    //     count = histOtof->GetBinContent(i);
    //     tof = histOtof->GetXaxis()->GetBinCenter(i);
    //     lambdaUpperLim = tof * factor_tof2lambda_max * 1e-6; // tof is written in µs
    //     for (int j = 0; j < count; j++)
    //         treeUpperLim.Fill();
    // }

    // treeUpperLim.Write();
    // outputUpperLim.Close();

    input->Close();

    return 0;
}

int main()
{
    return lambda_roi_histogram();
}
