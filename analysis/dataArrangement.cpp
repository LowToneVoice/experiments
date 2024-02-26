#include "../simulation/simulate.h"

int dataArrangement(
    string inputNoCutNoExt,
    string inputCutNoExt,
    double delta_deg,
    double time_sec
)
{
    // files
    TFile *file_noCut = TFile::Open(inputNoCutNoExt);
    TFile *file_Cut = TFile::Open(inputCutNoExt);
    TTree *tree, *treeCut, *treeNoCut;
    file_noCut->GetObject("tree", treeNoCut);
    file_Cut->GetObject("tree", treeCut);

    Int_t channel;
    Double_t lambda;
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lambda);
}
