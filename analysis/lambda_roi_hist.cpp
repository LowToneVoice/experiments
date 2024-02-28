#include "config.h"

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
    double factor2tof = h / total_length / m;

    // open the datafiles
    TFile *input = TFile::Open(input_tree);
    if (!input||input->IsZombie())
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
    TFile output(output_tree, "recreate");
    TTree tree("tree", "tree");
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lambda);

    // canvases and regions
    TCanvas *ctof = new TCanvas("ctof", "count on tof", 600, 400);
    TCanvas *ctof_perSec = new TCanvas("ctof", "count on tof", 600, 400);
    TCanvas *clmd = new TCanvas("clmd", "count on wavelength", 600, 400);
    TCut xCutH = ("x>=" + to_string(x_min_H) + "&&x<=" + to_string(x_max_H)).c_str();
    TCut yCutH = ("y>=" + to_string(y_min_H) + "&&y<=" + to_string(y_max_H)).c_str();
    TCut xCutO = ("x>=" + to_string(x_min_O) + "&&x<=" + to_string(x_max_O)).c_str();
    TCut yCutO = ("y>=" + to_string(y_min_O) + "&&y<=" + to_string(y_max_O)).c_str();


    // CALCULATION ON CHANNEL H
    // tof-count draw
    channel = H_TDC;
    TH1D *histHtof = new TH1D("histHtof", "tof_count (H: blue, O: red); tof (µs); count", 500, 0, 40000);
    TH1D *histHtof_perSec = new TH1D("histHtof_perSec", "tof count per second (H: blue, O: red); tof (µs); count (/s)", 500, 0, 40000);
    ctof->cd();
    T->Draw("tof>>histHtof", "f==4" && xCutH && yCutH);
    ctof_perSec->cd();
    T->Draw("tof>>histHtof_perSec", "f==4" && xCutH && yCutH, "same eh");
    // divide by time
    for (int bin = 1; bin <= histHtof_perSec->GetXaxis()->GetNbins(); bin++)
    {
        histHtof_perSec->SetBinContent(bin, histH_perSec->GetBinContent(bin) / time_sec);
        histHtof_perSec->SetBinError(bin, histHtof_perSec->GetBinError(bin) / time_sec);
    }
    // filling tree
    for (int i = 0; i < histHtof->GetXaxis()->GetNbins(); i++)
    {
        count = histHtof->GetBinContent(i);
        tof = histHtof->GetXaxis()->GetBinCenter(i);
        lambda = tof * 1e-6 * factor2tof;   // tof is written in µs
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
        histOtof_perSec->SetBinContent(bin, histO_perSec->GetBinContent(bin) / time_sec);
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

    ctof_perSec->Print(output_img_tof);
    clmd->Print(output_img_lmd);


    // ENDING
    tree.Write();
    input->Close();
    output.Close();
    return 0;
}
