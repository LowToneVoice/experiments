#include<iostream>
#include<TFile.h>
#include<TH1F.h>
#include<TTree.h>
#include<TCanvas.h>

// datafiles
#define INPUT_TREE "./BL05/20231224225335_list_copy.root"
#define OUTPUT_FILE "./beam_count/test/20231224225335.pdf"

#define h 6.62607015e-34
#define m 1.6749e-27

constexpr double total_length = 1.;
constexpr double d_total_length = .005;

constexpr double factor_tof2lambda = h / total_length / m;

int lambda_roi_histogram()
{
    // open the datafile
    TFile *file = TFile::Open(INPUT_TREE);
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }

    TTree *T;
    file->GetObject("T", T);
    if(!T)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }


    TCanvas *c1 = new TCanvas("c1", "count on wave length", 600, 400);
    // Draw the histogram with the rescaled x-axis using the added branch
    T->Draw("tof * tof_factor>>h_tof(500,0,40000)", "x>=0.41&&x<=0.44&&y>=0.34&&y<=0.48", "");

    return 0;
}

int main()
{
    return lambda_roi_histogram();
}
