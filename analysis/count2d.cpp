#include"config.h"

#define input_default "../simulation/BL05/20231224225335_list_copy.root"
#define output_default "./beam_count/position/20231224225335.pdf"

int count2d(
    string input_tree = input_default,
    string output_img_file = output_default
)
{
    // open the datafile
    TFile *file = TFile::Open(input_tree.c_str());
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the input file" << std::endl;
        return 1;
    }

    TTree *T;
    file->GetObject("T", T);
    if (!T)
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
    else
    {
        std::cerr << "There is no histogram" << std::endl;
    }

    c1->Print(output_img_file.c_str());

    return 0;
}

