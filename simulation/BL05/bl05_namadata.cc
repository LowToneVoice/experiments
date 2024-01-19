#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"

constexpr double mass = 1.675e-27; // 中性子の質量
constexpr double hbar = 1.055e-34;
constexpr double L_tof = 18022e-03;
constexpr double TEISU = 2. * 3.14159 * hbar / (mass * L_tof);

void bl05_namadata(){
    double MYPI = 3.1415;
    double teisu = 2.*MYPI*hbar/(mass*L_tof) *1.e-06;//tofを波長にする係数

    TFile *file = new TFile("./bl05/20210208160817_list.root");//ファイルの読み込み
    TTree *tree = (TTree*) file -> Get("T");

    TH1D *h1 = new TH1D("h1", "h1", 1000, 0, 45000*teisu);//空のヒストグラム

    tree -> Draw(Form("%e*tof>>h1", teisu), "f==4");//波長分布の書き込み

    double maxkp = tree -> GetMaximum("kp");
    h1 -> Scale(25./maxkp);//縦軸をcpsに変換

    h1 -> SetTitle("BL05 wavelength [cps]");
    h1 -> SetXTitle("lambda[m]");
    h1 -> SetYTitle("count per second");
    h1 -> Draw();

    h1 -> SaveAs("./bl05/BL05_wavelength.root");

    std::cout << 45000*teisu <<std::endl;
}
