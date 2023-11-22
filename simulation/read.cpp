#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TH1F.h>

#define INPUT_FILE_O "dat/outO_g_mix.dat"
#define INPUT_FILE_H "dat/outH_g_mix.dat"
#define MAX_DATA_SIZE 1e+6

// physical constants
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34

// alignment constants
#define Lambda_min 2e-10
#define Lambda_max 9e-10
#define TOTAL_LENGTH 1
#define DAQ_FREQ 62.5e+6 / 256

void read()
{
    // alignment calculations & settings
    double dt = 1 / DAQ_FREQ;
    double dLambda = h / TOTAL_LENGTH / m * dt;
    int L_max = (int)((Lambda_max - Lambda_min) / dLambda);
    L_max = 10000;
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
        return 1;
    }
    else
    {
        std::cout << "FILE OPENED" << std::endl;
    }

    // make histogram
    TH1F *histO = new TH1F("Obeam", "Obeam", L_max, Lambda_min, Lambda_max);
    TH1F *histH = new TH1F("Hbeam", "Hbeam", L_max, Lambda_min, Lambda_max);

    double lambda;
    // double lambdas[L_max] = {};
    // double values[L_max] = {};

    // for (int L = 0; L < L_max; L++)
    // {
    //     lambdas[L] = Lambda_min + (Lambda_max - Lambda_min) * L / L_max;
    // }

    // while (inputFileO >> lambda)
    // {
    //     histO->Fill(lambda);
    //     // values[int((lambda - Lambda_min) / (Lambda_max - Lambda_min) * L_max)] += 1.;
    // }
    int count = 0;
    while (inputFileH >> lambda)
    {
        histH->Fill(lambda);
        count++;
        if (count > MAX_DATA_SIZE)break;
        // values[int((lambda - Lambda_min) / (Lambda_max - Lambda_min) * L_max)] += 1.;
    }
    // TGraph *gr = new TGraph(L_max, lambdas, values);
    // gr->Draw("L");

    TF1 *fHG = new TF1("fH", "[ 0 ] * cos([ 1 ] * x) + [ 2 ]");
    fHG->SetParameters(1e+3, 5e+11, 1e+3);
    fHG->SetParNames("A", "k_{2}", "Const");
    histH->Fit("fH", "", "", Lambda_min, Lambda_max);
    // // TF1 *fOA = new TF1("fA", "5000 * cos([ 0 ] / x) + 5000");
    // // fOA->SetParameter(0, 3e-8);
    // // fOA->SetParName(0, "k_{1}");
    // // histO->Fit("fA", "", "", Lambda_min, Lambda_max);
    // TF1 *f = new TF1("f", "500 * cos([ 0 ] / x + [ 1 ] * x) + 500");
    // f->SetParameters(3e-8, 7e+13);
    // f->SetParNames("k_{1}", "k_{2}");
    // histH->Fit("fA", "", "", Lambda_min, Lambda_max);

    // draw histogram
    // histO->Draw();
    histH->Draw();

    inputFileO.close();
    inputFileH.close();
}
