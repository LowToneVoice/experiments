#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#define INPUT_FILE_O "dat/outO_mix_mix.dat"
#define INPUT_FILE_H "dat/outH_mix_mix.dat"

// physical constants
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34

// alignment constants
#define Lambda_min 2e-10
#define Lambda_max 9e-10
#define TOTAL_LENGTH 1
#define DAQ_FREQ 62.5e+6 / 256

#define DATASIZE 1e+18

void tree_make()
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

    // open the datafile
    std::cout << "///// FILE OPEN /////" << std::endl;
    std::ifstream inputFileO(INPUT_FILE_O);
    std::ifstream inputFileH(INPUT_FILE_H);
    if (!inputFileO.is_open() || !inputFileH.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return;
    }
    else
    {
        std::cout << "FILE OPENED" << std::endl;
    }

    // tree making
    std::cout << "///// DATA INPUT /////" << std::endl;
    Double_t lambda;
    Int_t channel;
    int count = 0;

    TFile *outputFileH = new TFile("tree.root", "recreate");
    if (!outputFileH || !outputFileH->IsOpen())
    {
        std::cerr << "Failed to create or open output file!" << std::endl;
        return;
    }
    TTree *treeH = new TTree("treeH", "treeH");
    if (!treeH || !treeO)
    {
        std::cerr << "Failed to create TTree!" << std::endl;
        outputFileH->Close();
        outputFileO->Close();
        delete outputFileH;
        delete outputFileO;
        return;
    }
    TH1F *histH = new TH1F("histH", "", L_max, Lambda_min, Lambda_max);
    histH->SetXTitle("lambda [m]");
    histH->SetYTitle("(I_H-I_O)/(I_H+I_O)");
    tree->Branch("lambda", &lambda);
    tree->Branch("channel", &channel);

    // H input
    std::cout << "H input start" << std::endl;
    channel = 0;
    while (inputFileH >> lambda)
    {
        histH->Fill(lambda);
        count++;
        if (count > DATASIZE)
        {
            std::cout << "File H has reached the max datasize." << std::endl;
            break;
        }
    }
    inputFileH.close();
    std::cout << "Input H has ended" << std::endl;

    // O files
    TFile *outputFileO = new TFile("tree.root", "recreate");
    TTree *treeO = new TTree("treeO", "treeO");
    TH1F *histO = new TH1F("histO", "", L_max, Lambda_min, Lambda_max);
    histO->SetXTitle("lambda [m]");
    histO->SetYTitle("(I_H-I_O)/(I_H+I_O)");

    // O input
    std::cout << "O input start" << std::endl;
    channel = 1;
    count = 0;
    while (inputFileO >> lambda)
    {
        histO->Fill(lambda);
        count++;
        if (count > DATASIZE)
        {
            std::cout << "File O has reached the max datasize." << std::endl;
            break;
        }
    }
    inputFileO.close();
    std::cout << "Input O has ended" << std::endl;

    // Save tree
    tree->Write();

    delete tree;
    delete outputFile;

    return;
}
