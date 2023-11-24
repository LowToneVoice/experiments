#include "toyMC.h"
#include <TSpline.h>

TH1D* GetFluxHist(){
  TFile *file = new TFile("./flux.root");
  if( !file->IsOpen() ) exit(0);
  TH1D *h = dynamic_cast<TH1D*>(file->Get("h_mr0"));
  /*
  TSpline3* sp = new TSpline3(h);
  sp->Print();
  TCanvas *c0 = new TCanvas("c0","c0",800,600);
  h->Draw("eh");
  sp->Draw("same");
  sp->SetLineColor(2);
  */
  return h;
}

TF1* InitFunc(Char_t* funcname=(Char_t*)"func",Double_t fit_begin=30., Double_t fit_end=50.){
  Double_t conv_ms = dist_detector/kConversion1*1.e3;
  Double_t tof_center = ref_lambda*conv_ms;

  TF1 *func = new TF1(funcname,
		      Form("[0]*(1-[3]/(1-[3])*(TMath::Sin(2*TMath::Pi()*([1]*%0.2f/x-[2]*x))-1))",conv_ms),
		      //		      Form("[0]*(1-[3]/(1-[3])*(TMath::Sin([1]*([2]-%0.2f/x))-1)",conv_ms),
		      fit_begin,fit_end);
  func->Print();
  //  TF1 *func = new TF1(funcname,Form("[0]*(1+[3]*TMath::Sin([1]*([2]-%0.2f/x)))+[4]*(x-%0.2f)",conv_ms,tof_center),fit_begin,fit_end);
  //  TF1 *func = new TF1(funcname,Form("[0]*(1+[3]*TMath::Sin([1]*([2]-%0.2f/(x+[5]))))+[4]*((x)-%0.2f)",conv_ms,tof_center),fit_begin,fit_end);
  func->Print();
  func->SetNpx(1000);
  func->SetParName(0, "Const.");
  func->SetParName(1, "2D#delta#theta [nm]");
  func->SetParName(2, "P_{1} [nm^{-1}]");
  func->SetParName(3, "Visibility");
  func->SetParName(4, "Slope [ms^{-1}]");
  
  //  func->SetParLimits(0,0.,1.);
  //  func->SetParLimits(1,-300.,300.);
  //  func->SetParLimits(3,0.,1.);
  //  func->SetParLimits(4,-0.06,0.0);

  func->SetParameter(0,0.4);
  func->SetParameter(1,150.);
  func->SetParameter(2,0.01);
  func->SetParameter(3,0.3);
  func->SetParameter(4,-0.003);
  
  return func;
}

TF1* FitFunc(Char_t* funcname=(Char_t*)"fitfunc",Double_t fit_begin=0.2, Double_t fit_end=1.){
  Double_t conv_ms = dist_detector/kConversion1*1.e3;
  //  Double_t tof_center = ref_lambda*conv_ms;
  Double_t tof_center = ref_lambda;
  
  TF1 *func = new TF1(funcname,
		      Form("[2]*TMath::Sin(2*TMath::Pi()*([0]*%0.2f/x-[1]*x))",conv_ms),
		      //		      "[2]*TMath::Sin([0]*x+[1]/x)",
		      fit_begin,fit_end);
  func->Print();
  func->SetNpx(1e6);
  func->SetParName(0, "2D#delta#theta [nm]");
  func->SetParName(1, "Phase [nm^{-1}]");
  func->SetParName(2, "Visibility");
  return func;
}

void WriteText(TF1* func, Char_t* comment){
  TString filename = "output.dat";
  ofstream outFile(filename.Data(), ios::app);
  outFile << comment <<"\t";
  outFile << func->GetParameter(0) <<"\t";
  outFile << func->GetParError(0) <<"\t";
  outFile << func->GetParameter(1) <<"\t";
  outFile << func->GetParError(1) <<"\t";
  outFile << func->GetParameter(2) <<"\t";
  outFile << func->GetParError(2) <<endl;
  outFile.close();
  return;
}

//Int_t nbin=1000;
//Double_t fit_begin   = 40.;
//Double_t fit_end     = 50.;
void toyMC(Int_t nbin = 1000, Double_t fit_begin = 12., Double_t fit_end = 50.){
  InitROOT();
  TGaxis::SetMaxDigits(3);
  gSystem->Exec("mkdir -p fig");

  TCanvas *c0 = new TCanvas();
  TH1D* hflux1 = GetFluxHist();
  TH1D* hflux2 = GetFluxHist();
  hflux1->Draw("eh");
  
  TCanvas *c1 = new TCanvas();
    //Fitting
  Double_t begin   =  10.;
  Double_t end     =  50.;
  Double_t fwidth  = end - begin; 
  Double_t Const   =  1.0;
  Double_t dD      = 3.;
  Double_t Phase   = 1.e-5;
  Double_t Vis     = 0.6;
  TF1* func1 = InitFunc((Char_t*)"func1",begin,end);
  func1->FixParameter(0,Const);
  func1->SetParameter(1,dD);
  func1->SetParameter(2,Phase);
  func1->SetParameter(3,Vis);
  TF1* func2 = InitFunc((Char_t*)"func2",begin,end);
  func2->FixParameter(0,Const);
  func2->SetParameter(1,-dD);
  func2->SetParameter(2,-Phase);
  func2->SetParameter(3,Vis);

  func2->SetLineColor(4);
  func1->Draw("");
  func2->Draw("same");
  c1->SaveAs("./fig/c1.pdf");

  TCanvas *c2 = new TCanvas();
  hflux1->Multiply(func1);
  hflux2->Multiply(func2);
  hflux1->SetLineColor(2);
  hflux2->SetLineColor(4);
  hflux1->SetLineStyle(2);
  hflux2->SetLineStyle(3);

  hflux1->Draw("h");
  hflux2->Draw("hsame");
  c2->SaveAs("./fig/c2.pdf");
  
  TCanvas *c3 = new TCanvas();
  //  Double_t hmin=0., hmax=50.;
  //  TH1D* ht1 = new TH1D("ht1","tof;tof [ms]; count[arb.]",nbin,hmin,hmax);
  //  TH1D* ht2 = new TH1D("ht2","tof;tof [ms]; count[arb.]",nbin,hmin,hmax);

  Double_t hmin=0., hmax=1.;
  TH1D* ht1 = new TH1D("ht1",";#lambda_{n} [nm]; count[1/pm/1000s]",nbin,hmin,hmax);
  TH1D* ht2 = new TH1D("ht2",";#lambda_{n} [nm]; count[1/pm/1000s]",nbin,hmin,hmax);

  Double_t hwidth = (hmax-hmin)/nbin;
  Double_t conv_ms = dist_detector/kConversion1*1.e3;
  //  const Int_t num = 1e5; //total count is 2 times
  //  const Int_t num = 7.5e6; // 1000s s equivarent with total flux is 1.5e4 count/s (Slit Full open/Collimator 0.1x10 mm)
  const Int_t num = 15e6; // 1000s s equivarent with total flux is 3e4 count/s (Slit Full open/Collimator 0.2x10 mm)
  for(Int_t i=0; i<num; i++){
    //    ht1->Fill( func1->GetRandom());
    //    ht2->Fill( func2->GetRandom());
    ht1->Fill( hflux1->GetRandom()/conv_ms);
    ht2->Fill( hflux2->GetRandom()/conv_ms);
  }
  ht1->SetLineColor(2);
  ht1->Draw("eh");
  ht2->Draw("ehsames");

  ht1->GetXaxis()->SetTitleSize(0.05);
  ht1->GetYaxis()->SetTitleSize(0.05);
  ht1->GetXaxis()->SetLabelSize(0.05);
  ht1->GetYaxis()->SetLabelSize(0.05);
  ht1->GetYaxis()->SetTitleOffset(0.9);
  ht2->SetLineStyle(2);
  ht1->SetStats(0);
  ht2->SetStats(0);
  TLegend* leg = new TLegend(0.6, 0.8, 0.95, 0.95,"");
  leg->AddEntry(ht1, "H beam", "pl");
  leg->AddEntry(ht2, "O beam", "pl");
  leg->Draw();

  c3->SaveAs("./fig/c3.pdf");

  TCanvas *c4 = new TCanvas();
  TH1D* ht3 = GetH1overH1pulseH2(ht1, ht2);
  //  TH1D* ht3 = GetH1minusH2overH1pulseH2(ht1, ht2);
  ht3->Print();
  ht3->SetLineColor(1);
  TH1D* htn = (TH1D*)ht3->Clone();
  htn->SetName("n");
  htn->Sumw2();
  htn->Scale(-2.);
  TF1* f0 = new TF1("f0","1",begin/conv_ms-1.e-3,1.);
  htn->Add(f0,1);
  htn->Draw("eh");
  htn->SetStats(1);
  htn->SetTitle("");
  htn->GetYaxis()->SetTitle("n  [(I_{h} - I_{o}) / (I_{h} + I_{o})]");
  htn->GetYaxis()->SetRangeUser(-1.2,+1.2);
  htn->GetXaxis()->SetTitleSize(0.05);
  htn->GetYaxis()->SetTitleSize(0.05);
  htn->GetYaxis()->SetTitleOffset(0.9);
  htn->GetXaxis()->SetLabelSize(0.05);
  htn->GetYaxis()->SetLabelSize(0.05);

  TF1* fitfunc = FitFunc((Char_t*)"fitfunc",0.225,1.0);
  fitfunc->SetParameter(0,dD/conv_ms);
  fitfunc->SetParameter(1,Phase);
  fitfunc->SetParameter(2,Vis);

  fitfunc->FixParameter(0,dD/conv_ms);
  fitfunc->SetParameter(1,Phase);
  fitfunc->FixParameter(2,Vis);
  
  htn->Fit("fitfunc","R");
  c4->SaveAs("./fig/c4.pdf");

  //  WriteText(fitfunc, (Char_t*)Form(" %0.2f %0.0f %0.0f %d #binwidth = %0.2f ms, Fitting %0.0f - %0.0f ms, n=%d",
  //				   hwidth,fit_begin,fit_end,num,hwidth,fit_begin,fit_end,num));

  WriteText(fitfunc, (Char_t*)Form("%0.2f %0.0f %0.0f %d",hwidth,fit_begin,fit_end,num));

  return;
}


