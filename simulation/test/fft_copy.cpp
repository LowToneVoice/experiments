#include<TH1D.h>
#include<TVirtualFFT.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TMath.h>

int fft_copy()
{
    // Histograms
    // =========
    // prepare the canvas for drawing
    TCanvas *myc = new TCanvas("myc", "Fast Fourier Transform", 800, 600);
    myc->SetFillColor(45);
    TPad *c1_1 = new TPad("c1_1", "c1_1", 0.01, 0.67, 0.49, 0.99);
    TPad *c1_2 = new TPad("c1_2", "c1_2", 0.51, 0.67, 0.99, 0.99);
    TPad *c1_3 = new TPad("c1_3", "c1_3", 0.01, 0.34, 0.49, 0.65);
    TPad *c1_4 = new TPad("c1_4", "c1_4", 0.51, 0.34, 0.99, 0.65);
    TPad *c1_5 = new TPad("c1_5", "c1_5", 0.01, 0.01, 0.49, 0.32);
    TPad *c1_6 = new TPad("c1_6", "c1_6", 0.51, 0.01, 0.99, 0.32);
    c1_1->Draw();
    c1_2->Draw();
    c1_3->Draw();
    c1_4->Draw();
    c1_5->Draw();
    c1_6->Draw();
    c1_1->SetFillColor(30);
    c1_1->SetFrameFillColor(42);
    c1_2->SetFillColor(30);
    c1_2->SetFrameFillColor(42);
    c1_3->SetFillColor(30);
    c1_3->SetFrameFillColor(42);
    c1_4->SetFillColor(30);
    c1_4->SetFrameFillColor(42);
    c1_5->SetFillColor(30);
    c1_5->SetFrameFillColor(42);
    c1_6->SetFillColor(30);
    c1_6->SetFrameFillColor(42);

    c1_1->cd();
    TH1::AddDirectory(kFALSE);

    // A function to sample
    TF1 *fsin = new TF1("fsin", "sin(x)+sin(2*x)+sin(0.5*x)+1", 0, 4 * TMath::Pi());
    fsin->Draw();

    Int_t n = 25;
    TH1D *hsin = new TH1D("hsin", "hsin", n + 1, 0, 4 * TMath::Pi());
    Double_t x;

    // Fill thew histogram with function values
    for (Int_t i = 0; i <= n; i++)
    {
        x = (Double_t(i) / n) * (4 * TMath::Pi());
        hsin->SetBinContent(i + 1, fsin->Eval(x));
    }
    hsin->Draw();
    fsin->GetXaxis()->SetLabelSize(.05);
    fsin->GetYaxis()->SetLabelSize(.05);

    c1_2->cd();
    TH1 *hm = nullptr;
    TVirtualFFT::SetTransform(nullptr);
    hm = hsin->FFT(hm, "MAG");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->SetStats(kFALSE);

    // Print the frequencies corresponding to each bin
    // Print the frequencies corresponding to each bin within the expected range
    Double_t binWidth = 1.0 / (4 * TMath::Pi() * n);
    Double_t nyquistFreq = 0.5 / binWidth;

    for (Int_t i = 0; i <= n; i++)
    {
        Double_t frequency = i * binWidth;

        // Skip frequencies beyond the Nyquist frequency
        if (frequency > nyquistFreq)
        {
            break;
        }

        Double_t magnitude = hm->GetBinContent(i + 1);
        std::cout << "Frequency: " << frequency << ", Magnitude: " << magnitude << std::endl;
    }

    hm->Draw();

    return 1;
}
