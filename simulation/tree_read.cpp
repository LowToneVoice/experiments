#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

// physical constants
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34

// alignment constants
#define Lambda_min 2e-10
#define Lambda_max 9e-10
#define TOTAL_LENGTH 1
#define DAQ_FREQ 62.5e+6 / 256

#define DATASIZE 1e+15

void tree_read()
{
    // PREPARATION
    // Alignment calculations & settings
    double dt = 1 / DAQ_FREQ;
    double dLambda = h / TOTAL_LENGTH / m * dt;
    int L_max = (int)((Lambda_max - Lambda_min) / dLambda);
    // L_max = 10000;
    std::cout << "///// ALIGNMENT CALCULATIONS /////" << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "dLambda = " << dLambda << std::endl;
    std::cout << "L_max = " << L_max << std::endl;

    // Open datafile
    TFile *file = TFile::Open("tree.root");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Unable to open the file 'tree.root'" << std::endl;
        return;
    }
    TTree *tree;
    file->GetObject("tree", tree);
    if (!tree)
    {
        std::cerr << "Error: Unable to retrieve the tree from the file" << std::endl;
        file->Close();
        delete file;
        return;
    }

    Int_t channel;
    Double_t lambda;

    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("lambda", &lambda);
    Long64_t totalEntries = tree->GetEntries();

    std::cout << "Entries in the tree: " << totalEntries << std::endl;
    // std::cout << "---------------------------------" << std::endl;
    // std::cout << "  Channel   |   Lambda" << std::endl;
    // std::cout << "---------------------------------" << std::endl;

    // for (Long64_t i = 0; i < totalEntries && i < DATASIZE; ++i)
    // {
    //     tree->GetEntry(i);
    //     std::cout << "    " << channel << "       |    " << lambda << std::endl;
    // }

    tree->Draw("lambda", "channel==0");

    file->Close();
    delete file;
}
