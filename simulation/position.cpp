#include<iostream>
#include<TFile.h>
#include<TH1F.h>
#include<TTree.h>
#include<TCanvas.h>

// datafiles
#define INPUT_TREE "./BL05/20231224225335_list_copy.root"
#define OUTPUT_FILE "./beam_count/position/20231224225335.pdf"

int position()
{
    // open the datafile
    TFile *file = TFile::Open(INPUT_TREE);
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file" << std::endl;
        return 1;
    }

    TTree *T;
    file->GetObject("T", T);
    if(!T)
    {
        std::cerr << "Failed to retrieve the tree from the file." << std::endl;
        return 1;
    }

    TCanvas *c1 = new TCanvas("c1", "Beam count at xy position", 600, 400);
    T->Draw("y:x", "f==4&&x>=0&&x<=1&&y>=00&&y<=1", "colz");

    // The following is only in order to name the graph.
    TH1F *histogram = (TH1F *)gPad->GetPrimitive("htemp");
    if (histogram)
        histogram->SetTitle("Beam count at xy position");

    c1->Print(OUTPUT_FILE);

    return 0;
}

int main()
{
    return position();
}
