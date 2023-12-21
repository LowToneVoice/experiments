#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>

// datafiles
#define INPUT_TREE "./dat/montecarlo/noref/g_main.root"

// physical constants
#define pi M_PI
#define m 1.6749e-27
#define h 6.62607015e-34
#define DATASIZE 1e6
constexpr double g = 9.8;

// alignment constants
constexpr double Lambda_min = 2e-10;
constexpr double Lambda_max = 12e-10;
constexpr double total_length = 1.;
constexpr double daq_freq = 62.5e6;
constexpr int daq_downsizing = 256;

int oscillation()
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
        file->Close();
        delete tree;
        return 1;
    }

    Int_t channel;
    Double_t lambda, theta;
    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("lambda", &lambda);
    tree->SetBranchAddress("theta", &theta);

    tree->Draw("lambda", "theta==1.05 && channel==0");

    return 0;
}
