#include<TH1D.h>
#include<TVirtualFFT.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TMath.h>

int fft_copy()
{
    TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 800, 600);
    TPad *c1_1 = new TPad("c1_1", "c1_1", .01, .67, .49, .99);
    TPad *c1_2 = new TPad("c1_2", "c1_2", .51, .67, .99, .99);
    TPad *c1_3 = new TPad("c1_3", "c1_3", .01, .34, .49, .65);
    TPad *c1_4 = new TPad("c1_4", "c1_4", .51, .34, .99, .65);
    TPad *c1_5 = new TPad("c1_5", "c1_5", .01, .01, .49, .32);
    TPad *c1_6 = new TPad("c1_6", "c1_6", .51, .01, .99, .32);
    c1_1->Draw();
    c1_2->Draw();
    c1_3->Draw();
    c1_4->Draw();
    c1_5->Draw();
    c1_6->Draw();

    c1_1->cd();
    TH1::AddDirectory(kFALSE);

    // A function to sample
    TF1 *fsin = new TF1("fsin", "sin(x)+sin(2*x)+sin(0.5*x)+1", 0, 4 * TMath::Pi());
    fsin->Draw();

    Int_t n = 30;
    TH1D *hsin = new TH1D("hsin", "hsin", n + 1, 0, 4 * TMath::Pi());
    Double_t x;

    // Fill the histogram with function values
    for (Int_t i = 0; i <= n; i++)
    {
        x = (Double_t(i) / n) * (4 * TMath::Pi());
        hsin->SetBinContent(i + 1, fsin->Eval(x));
    }
    hsin->Draw("same");

    // Compute the transform and look at the magnitude of the output
    c1_2->cd();
    TH1 *hm = nullptr;
    TVirtualFFT::SetTransform(nullptr);
    hm = hsin->FFT(hm, "MAG");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->SetStats(kFALSE);
    // NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function (in this case 4*Pi);
    // y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
    // hm->GetXaxis()->SetLimits(0, 1.0 / (4 * TMath::Pi()));
    hm->Draw();

    return 1;
}
