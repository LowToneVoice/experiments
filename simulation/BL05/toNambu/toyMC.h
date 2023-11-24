#include <sstream>
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"    
#include "TF1.h" 
#include "TF2.h" 
#include "TMath.h"
#include "TLegend.h"
#include "TPaveText.h"                                                         
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TCut.h"
#include "TGaxis.h"
#include "TSystem.h"

#define sq(x)  ((x)*(x))

#include <TH3.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TLine.h>

const Bool_t use_set_max_entry_loop = false;
const Int_t  max_entry_loop = 1E5;

// Physical constants
const Double_t kMn = 1.6749275e-27;  // kg
const Double_t kConversion1 = 395.6; // lambda = tof*1.e-6*Conversion1/Distance;
const Double_t kConversion2 = 8.0*kMn*TMath::Qe()*1e-9/pow(TMath::H()/2.0/TMath::Pi(), 2.0)*1e-18;
// Conversion of potential energy to q^2

// Etalons
const Double_t dist_etalons = 142.0; //[mm]
const Double_t ref_lambda   = 0.88; // mirror reflective wavelength [nm]

// RPMT
const Double_t dist_detector = 17.75;
const Int_t    nbins_x  = 640;
const Double_t min_x    = 0.;
const Double_t max_x    = 128.;
const Int_t    nbins_y  = 640;
const Double_t min_y    = 0.;
const Double_t max_y    = 128.;
const Int_t    nbins_t  = 200;
// const Int_t    nbins_t  = 400;
const Double_t min_t    = 0.;
#if MRCUT
const Double_t max_t    = 120.; //for nuke pulse 
const Double_t flap_t   = 0.; //for nuke pulse
TCut cut_tof    = "MRflag>0";
#else
const Double_t max_t    = 40.;
const Double_t flap_t   = 10.;
TCut cut_tof    = "";
#endif

// Calculation conditions
const Int_t    LLD            = 500.;
//  const Double_t HLD  = 7400.;
const Double_t rpmt_range     = 128.0;

Double_t roi_xmin    = 60.0;
Double_t roi_xmax    = 100.0;
Double_t roi_ymin    = 56.0;
Double_t roi_ymax    = 75.0;
Double_t direct_xmin = 69.0;
Double_t direct_xmax = 73.0;
// theta=-1.05 condition
Double_t O_xmin      = 75.0;
Double_t O_xmax      = 77.5;
Double_t H_xmin      = 90.5; // 0617 overnight
Double_t H_xmax      = 93.0; // 0617 overnight

//Cuts
TCut cut_direct = Form("x*128.>%f&&x*128.<%f", direct_xmin, direct_xmax);
TCut cut_H      = Form("x*128.>%f&&x*128.<%f", H_xmin, H_xmax);
TCut cut_O      = Form("x*128.>%f&&x*128.<%f", O_xmin, O_xmax);
TCut cut_roi    = Form("x*128.>%f&&x*128.<%f&&y*128.>%f&&y*128.<%f",
			 roi_xmin, roi_xmax, roi_ymin, roi_ymax);

TCanvas* GetSquarePlot (TString name, TString title,
                        int w, float l, float r, float b, float t) {
  // https://root-forum.cern.ch/t/draw-square-th2/29135/6
  // This function creates a canvas making sure the plot is square.
  // The first parameter w is the canvas width and l, r, b and t are the
  // left, right, bottom and top magins. The canvas height is computed with these
  // parameters in order to have a square plot.
  int h = ((1.-(l+r))*w)/(1.-(b+t));
  TCanvas *c = new TCanvas(name,title,w,h);
  c->SetLeftMargin(l),
    c->SetRightMargin(r),
    c->SetBottomMargin(b),
    c->SetTopMargin(t);
  // c->Draw();
  return c;
}
Double_t GetUCoord(TString axis, TString minmax) {
  Double_t axisminmax;
  if      (axis == "x" && minmax == "min") axisminmax = gPad->GetUxmin();
  else if (axis == "x" && minmax == "max") axisminmax = gPad->GetUxmax();
  else if (axis == "y" && minmax == "min") axisminmax = gPad->GetUymin();
  else if (axis == "y" && minmax == "max") axisminmax = gPad->GetUymax();
  else {
    std::cerr << "axis=" << axis
              << " minmax=" << minmax << " is invalid" << std::endl;
    exit(1);
  }
  if      (axis == "x" &&  gPad->GetLogx()) return pow(10., axisminmax);
  else if (axis == "x" && !gPad->GetLogx()) return axisminmax;
  else if (axis == "y" &&  gPad->GetLogy()) return pow(10., axisminmax);
  else if (axis == "y" && !gPad->GetLogy()) return axisminmax;
  else {
    std::cerr << "axis gPad->GetLogx or y error" << std::endl;
    exit(1);
  }
}
void InitROOT() {
  // ======== initial setting ===============
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);
  TH1::StatOverflows(kTRUE);
  gStyle->SetOptStat(1001111);//to delete number pallete
  gStyle->SetOptFit(1110);
  //  gStyle->SetOptStat(0);
}

TGraphErrors *ConvertHistToGraph(TH1D* hist) {
  TGraphErrors* gr = new TGraphErrors(hist);
  //remove t0
  gr->RemovePoint(0);
  Int_t num = gr->GetN();
  for(Int_t i=0; i<num; i++){
    Double_t x = gr->GetX()[i];
    Double_t y = gr->GetY()[i];
    Double_t ey = gr->GetErrorY(i);
    gr->SetPointError(i,0.,ey);
    if(flap_t == 0.) continue;
    else if(x<flap_t) {      
      gr->SetPoint(i,x+max_t,y);
    }
  }
  return gr;
}

Double_t GetAoverApulsBerror(Double_t A, Double_t B, Double_t EA, Double_t EB){
  Double_t e1 = TMath::Sqrt( EA*EA/(A*A) + EB*EB/(B*B) );
  Double_t e2 = A*B/(A+B)/(A+B);
  return e1*e2;
}

TH1D* GetH1overH1pulseH2(TH1D* h1, TH1D* h2){
  TH1D *hsum = (TH1D*)h1->Clone("hsum");
  TH1D *hr   = (TH1D*)h1->Clone("hr");
  hsum->Add(h2);
  hr->Divide(hsum);
  for (Int_t i = 0; i < hr->GetNbinsX()+1; i++) {    
    Double_t A = h1->GetBinContent(i);
    Double_t B = h2->GetBinContent(i);
    Double_t EA = h1->GetBinError(i);
    Double_t EB = h2->GetBinError(i);
    if(B==0.)continue;
    Double_t err_tot = GetAoverApulsBerror(A,B,EA,EB);
    hr->SetBinError(i,err_tot);
  }
  hr->SetName(Form("%s_%s",h1->GetName(),h2->GetName()));
  return hr;
}

TH1D* GetH1minusH2overH1pulseH2(TH1D* h1, TH1D* h2){
  TH1D *hsum = (TH1D*)h1->Clone("hsum");
  TH1D *hr   = (TH1D*)h1->Clone("hr");
  hsum->Add(h2);
  hr->Add(h2);
  hr->Divide(hsum);
  for (Int_t i = 0; i < hr->GetNbinsX()+1; i++) {    
    Double_t A = h1->GetBinContent(i);
    Double_t B = h2->GetBinContent(i);
    Double_t EA = h1->GetBinError(i);
    Double_t EB = h2->GetBinError(i);
    if(B==0.)continue;
    Double_t err_tot = GetAoverApulsBerror(A,B,EA,EB);
    hr->SetBinError(i,err_tot);
  }
  hr->SetName(Form("%s_%s",h1->GetName(),h2->GetName()));
  return hr;
}

