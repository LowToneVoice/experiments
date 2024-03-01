#include<config.h>

int oscil_cut(
    double time_noCut,
    double time_Cut,
    string input_noCut,
    string input_Cut,
    string output_normal,
    string output_zoom,
    string output_fourier,
    string output_fourier_wide,
    string histTitle = "NoTitle"
)
{
    // open the data file
    TFile *file_noCut = TFile::Open(input_noCut.c_str());
    TFile *file_Cut = TFile::Open(input_Cut.c_str());
    if(!file_noCut||!file_Cut||file_noCut->IsZombie()||file_Cut->IsZombie())
    {
        cerr << "Failed to open the file." << endl;
        return 0;
    }
    TTree *tree_noCut, *tree_Cut;
    file->GetObject("tree_noCut", tree_noCut);
    file->GetObject("tree_Cut", tree_Cut);
    if (!tree_noCut || !tree_Cut)
    {
        cerr << "Failed to retrieve the tree." << endl;
        return 0;
    }

    // conditions
    TCut CutHnoCut = ("channel==" + to_string(H_TDC)).c_str();
    TCut CutHCutC = ("channel==" + to_string(H_TDC_Ccut)).c_str();
    TCut CutHCutD = ("channel==" + to_string(H_TDC_Dcut)).c_str();
    TCut CutOnoCut = ("channel==" + to_string(O_TDC)).c_str();
    TCut CutHCutC = ("channel==" + to_string(O_TDC_Ccut)).c_str();
    TCut CutHCutD = ("channel==" + to_string(O_TDC_Dcut)).c_str();

    // make histogram
    TH1F *histO = new TH1F("histO", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
    TH1F *histOcutC = new TH1F("histOcutC", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
    TH1F *histOcutD = new TH1F("histOcutD", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
    TH1F *histH = new TH1F("histH", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
    TH1F *histHcutC = new TH1F("histHcutC", ";(I_H-I_O)/(I_H+I_O);lambda [m]", nBins, lambda_min, lambda_max);
    TH1F *histHcutD = new TH1F("histHcutD", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
    TH1F *hist_num = new TH1F("hist_num", ";lambda [m];(I_H-I_O)/(I_H+I_O)", nBins, lambda_min, lambda_max);
}
