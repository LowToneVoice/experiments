// from 2021 P2

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

// SI単位系
constexpr double V_sio2 = 90.9e-09 * 1.6022e-19;

constexpr double V_per45 = 297.4e-09 * 1.6022e-19;
constexpr double V_per45_d = 104.4e-09 * 1.6022e-19;

constexpr double V_ge = 94.1e-09 * 1.6022e-19;

constexpr double D_per45 = 154.8e-10 / 0.93; // 厚さが1/0.93倍になっていたことに気付く
constexpr double D_ge = 108.4e-10 / 0.93;

constexpr double V_ni = 243.4e-09 * 1.6022e-19;
constexpr double V_ti = 0.0;

constexpr double D_ni = 133.5e-10 / 0.93;
constexpr double D_ti = 98.3e-10 / 0.93;

constexpr int Nlayer = 8;

constexpr double hbar = 1.055e-34;
constexpr double m = 1.675e-27;
constexpr double g = 9.8;
constexpr double PI = 3.14159265;

constexpr double theta_min = 0.6;
constexpr double theta_max = 1.5;
constexpr double d_theta = 0.001;
constexpr int Nroop = (int)((theta_max - theta_min) / d_theta);

/** 各層の反射行列のようなものの計算 **/
std::vector<std::vector<double>> matrix(double k0, double theta, double V, double D)
{
    std::vector<std::vector<double>> M(2, std::vector<double>(2, 0.));

    double kx_0 = k0 * std::sin(theta);
    double kx_in_square = kx_0 * kx_0 - 2. * m * V / (hbar * hbar);

    if (kx_in_square < 0.0)
    {
        // nが虚数の時
        double kx_in = sqrt(-kx_in_square);
        double nx = kx_in / kx_0;
        double delta = kx_0 * D;

        M.at(0).at(0) = std::cosh(nx * delta);
        M.at(0).at(1) = std::sinh(nx * delta) / nx;
        M.at(1).at(0) = std::sinh(nx * delta) * nx;
        M.at(1).at(1) = std::cosh(nx * delta);
    }
    else
    {
        // nが実数の時
        double kx_in = sqrt(kx_in_square);
        double nx = kx_in / kx_0;
        double delta = kx_0 * D;

        M.at(0).at(0) = std::cos(nx * delta);
        M.at(0).at(1) = std::sin(nx * delta) / nx;
        M.at(1).at(0) = -std::sin(nx * delta) * nx;
        M.at(1).at(1) = std::cos(nx * delta);
    }

    return M;
}

/** 反射率の計算 **/
double reflect(double k0, double rad, double Vg, double A, double B, double C, double D)
{
    double kx_0 = k0 * std::sin(rad);
    double kx_g_square = kx_0 * kx_0 - 2. * m * Vg / (hbar * hbar);

    double R;

    if (kx_g_square < 0.0)
    {
        // ngが虚数の時
        double kx_g = sqrt(-kx_g_square);
        double n0 = 1.;
        double ng = kx_g / kx_0;

        std::complex<double> ue(-ng * A - C, -n0 * ng * B - n0 * D);
        std::complex<double> sita(ng * A + C, -n0 * ng * B - n0 * D);
        std::complex<double> r = ue / sita;
        R = std::norm(r);
    }
    else
    {
        // ngが実数の時
        double kx_g = sqrt(kx_g_square);
        double n0 = 1.;
        double ng = kx_g / kx_0;

        std::complex<double> ue(-n0 * ng * B - C, ng * A - n0 * D);
        std::complex<double> sita(-n0 * ng * B + C, -ng * A - n0 * D);
        std::complex<double> r = ue / sita;
        R = std::norm(r);
    }

    return R;
}

/** 行列の積の計算　eigenってライブラリ使えばよかった **/
std::vector<std::vector<double>> DOT(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N)
{
    std::vector<std::vector<double>> product(2, std::vector<double>(2, 0.));
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int l = 0; l < 2; l++)
            {
                product.at(i).at(j) += M.at(i).at(l) * N.at(l).at(j);
            }
        }
    }
    return product;
}

void reffortheta()
{

    double lambda = 8.0e-10;      // 波長
    double k0 = 2. * PI / lambda; // 波数、2πがつくことに注意
    double theta[Nroop];          // 入射角
    double R[3][Nroop];           // 反射率

    for (int j = 0; j < Nroop; j++)
    {

        theta[j] = theta_min + d_theta * j; // deg
        double rad = theta[j] * PI / 180.;  // rad

        /** 磁気ミラーの計算**/
        std::vector<std::vector<double>> M_per45 = matrix(k0, rad, V_per45, D_per45); // permalloy45の行列
        std::vector<std::vector<double>> M_ge = matrix(k0, rad, V_ge, D_ge);          // Geの行列

        std::vector<std::vector<double>> layermatrix_per_ge(2, std::vector<double>(2, 0.)); // 二個分の行列
        std::vector<std::vector<double>> refmatrix_per_ge(2, std::vector<double>(2, 0.));   // 多層膜の行列

        layermatrix_per_ge = DOT(M_ge, M_per45);

        // std::cout << layermatrix.at(0).at(0) <<" "<< layermatrix.at(0).at(1) <<" "<< layermatrix.at(1).at(0) <<" "<< layermatrix.at(1).at(1) << std::endl;

        // まずは単位行列
        refmatrix_per_ge.at(0).at(0) = 1.;
        refmatrix_per_ge.at(0).at(1) = 0.;
        refmatrix_per_ge.at(1).at(0) = 0.;
        refmatrix_per_ge.at(1).at(1) = 1.;

        for (int l = 0; l < Nlayer; l++)
        {
            refmatrix_per_ge = DOT(layermatrix_per_ge, refmatrix_per_ge);
        }

        double A = refmatrix_per_ge.at(0).at(0);
        double B = refmatrix_per_ge.at(0).at(1);
        double C = refmatrix_per_ge.at(1).at(0);
        double D = refmatrix_per_ge.at(1).at(1);

        R[0][j] = reflect(k0, rad, V_sio2, A, B, C, D);

        // std::cout << R[1][j] << std::endl;

        /**全反射ミラーの計算**/
        std::vector<std::vector<double>> M_ni = matrix(k0, rad, V_ni, D_ni); // permalloy45の行列
        std::vector<std::vector<double>> M_ti = matrix(k0, rad, V_ti, D_ti); // Geの行列

        std::vector<std::vector<double>> layermatrix_ni_ti(2, std::vector<double>(2, 0.)); // 二個分の行列
        std::vector<std::vector<double>> refmatrix_ni_ti(2, std::vector<double>(2, 0.));   // 多層膜の行列

        layermatrix_ni_ti = DOT(M_ti, M_ni);

        // std::cout << layermatrix.at(0).at(0) <<" "<< layermatrix.at(0).at(1) <<" "<< layermatrix.at(1).at(0) <<" "<< layermatrix.at(1).at(1) << std::endl;

        // まずは単位行列
        refmatrix_ni_ti.at(0).at(0) = 1.;
        refmatrix_ni_ti.at(0).at(1) = 0.;
        refmatrix_ni_ti.at(1).at(0) = 0.;
        refmatrix_ni_ti.at(1).at(1) = 1.;

        for (int l = 0; l < Nlayer; l++)
        {
            refmatrix_ni_ti = DOT(layermatrix_ni_ti, refmatrix_ni_ti);
        }

        double E = refmatrix_ni_ti.at(0).at(0);
        double F = refmatrix_ni_ti.at(0).at(1);
        double G = refmatrix_ni_ti.at(1).at(0);
        double H = refmatrix_ni_ti.at(1).at(1);

        R[1][j] = reflect(k0, rad, V_sio2, E, F, G, H);

        // 磁気ミラーのdownに対しての反射率
        std::vector<std::vector<double>> M_per45_d = matrix(k0, rad, V_per45_d, D_per45); // permalloy45の行列
        std::vector<std::vector<double>> M_ge_d = matrix(k0, rad, V_ge, D_ge);            // Geの行列

        std::vector<std::vector<double>> layermatrix_per_ge_d(2, std::vector<double>(2, 0.)); // 二個分の行列
        std::vector<std::vector<double>> refmatrix_per_ge_d(2, std::vector<double>(2, 0.));   // 多層膜の行列

        layermatrix_per_ge_d = DOT(M_ge, M_per45_d);

        // std::cout << layermatrix.at(0).at(0) <<" "<< layermatrix.at(0).at(1) <<" "<< layermatrix.at(1).at(0) <<" "<< layermatrix.at(1).at(1) << std::endl;

        // まずは単位行列
        refmatrix_per_ge_d.at(0).at(0) = 1.;
        refmatrix_per_ge_d.at(0).at(1) = 0.;
        refmatrix_per_ge_d.at(1).at(0) = 0.;
        refmatrix_per_ge_d.at(1).at(1) = 1.;

        for (int l = 0; l < Nlayer; l++)
        {
            refmatrix_per_ge_d = DOT(layermatrix_per_ge_d, refmatrix_per_ge_d);
        }

        double I = refmatrix_per_ge_d.at(0).at(0);
        double J = refmatrix_per_ge_d.at(0).at(1);
        double K = refmatrix_per_ge_d.at(1).at(0);
        double L = refmatrix_per_ge_d.at(1).at(1);

        R[2][j] = reflect(k0, rad, V_sio2, I, J, K, L);
    }

    TCanvas *c_reffortheta1 = new TCanvas("c_reffortheta1", "c_reffortheta1");
    TMultiGraph *mg1 = new TMultiGraph();
    TGraph *graph1 = new TGraph(Nroop, theta, R[0]);
    TGraph *graph1_d = new TGraph(Nroop, theta, R[2]);
    mg1->SetTitle(Form("Permalloy45/Ge at %f angstrom;incident angle [degree];reflectivity", lambda * 1.e10));
    graph1->SetLineColor(1);
    graph1_d->SetLineColor(2);
    mg1->Add(graph1);
    mg1->Add(graph1_d);
    mg1->Draw("AL");

    TCanvas *c_reffortheta2 = new TCanvas("c_reffortheta2", "c_reffortheta2");
    TGraph *graph2 = new TGraph(Nroop, theta, R[1]);
    graph2->SetTitle(Form("Ni/Ti at %f angstrom", lambda * 1.e10));
    graph2->GetXaxis()->SetTitle("incident angle [degree]");
    graph2->GetYaxis()->SetTitle("reflectivity");
    graph2->Draw("AL");
}
