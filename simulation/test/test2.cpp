// from 2021 P2

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMultiGraph.h"

// PHYS & MATH CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define pi M_PI
#define PI M_PI
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
#define J_per_eV 1.6022e-19
constexpr std::complex<double> i = std::complex<double>(0, 1);
constexpr double hbar = h / 2 / pi;
constexpr double V_sio2 = 90.9e-09 * J_per_eV;
constexpr double V_per45 = 297.4e-09 * J_per_eV;
constexpr double V_per45_d = 104.4e-09 * J_per_eV;
constexpr double V_ge = 94.1e-09 * J_per_eV;
constexpr double V_ni = 224.e-09 * J_per_eV;
constexpr double V_ti = -40.e-9 * J_per_eV;
constexpr double rho_Ni = 9.1e28;
constexpr double bc_Ni = 3.4e-15;

constexpr double delta_D = 1;
// 厚さが1/0.93倍になっていたことに気付く
// constexpr double delta_D = 0.93;

constexpr double D_per45 = 154.8e-10 / delta_D;
constexpr double D_ge = 108.4e-10 / delta_D;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;

constexpr int N_bilayer = 8;

constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 12e-10;
constexpr double d_lambda = 0.01e-10;

constexpr int N_loop = (int)((lambda_max - lambda_min) / d_lambda);

/** 各層の反射行列のようなものの計算 **/
Eigen::Matrix2d mat_layer(double k0, double theta, double V, double D)
{
    Eigen::Matrix2d M;
    double kx_0 = k0 * std::sin(theta);
    double kx_in_square = kx_0 * kx_0 - 2. * m * V / (hbar * hbar);

    // n in real
    if (kx_in_square > 0.)
    {
        double kx_in = sqrt(kx_in_square);
        double nx = kx_in / kx_0;
        double delta = kx_0 * D;

        M(0, 0) = std::cos(nx * delta);
        M(0, 1) = std::sin(nx * delta) / nx;
        M(1, 0) = -std::sin(nx * delta) * nx;
        M(1, 1) = std::cos(nx * delta);
    }
    // n in imaginary
    else
    {
        double kx_in = sqrt(-kx_in_square);
        double nx = kx_in / kx_0;
        double delta = kx_0 * D;

        M(0, 0) = std::cosh(nx * delta);
        M(0, 1) = std::sinh(nx * delta) / nx;
        M(1, 0) = std::sinh(nx * delta) * nx;
        M(1, 1) = std::cosh(nx * delta);
    }

    return M;
}
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

/**屈折率のRe,Imの判定　いつか使えるかも**/
int n_hantei(double k0, double rad, double V)
{
    double kx_0 = k0 * std::sin(rad);
    double kx_square = kx_0 * kx_0 - 2. * m * V / (hbar * hbar);
    if (kx_square < 0.0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

/** 反射率の計算 **/
/* 真空から入射して Vg の物質へ抜けるときの反射率 */
std::complex<double> reflect(double k0, double rad, double Vg /* always V_SiO2 */, double A, double B, double C, double D)
{
    double kx_0 = k0 * std::sin(rad);
    double kx_g_square = kx_0 * kx_0 - 2. * m * Vg / (hbar * hbar);

    if (kx_g_square < 0.0)
    {
        // ngが虚数の時
        double kx_g = sqrt(-kx_g_square);
        double n0 = 1.;
        double ng = kx_g / kx_0;

        std::complex<double> ue(-ng * A - C, -n0 * ng * B - n0 * D);
        std::complex<double> sita(ng * A + C, -n0 * ng * B - n0 * D);
        std::complex<double> r = ue / sita;
        return r;
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
        return r;
    }
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

void test2()
{

    double lambda[N_loop]; // 波長
    double theta = 1.05;         // 入射角
    double rad = theta * PI / 180.;
    double R[3][N_loop]; // 反射率
    double hantei[2][N_loop] = {0.};
    std::complex<double> r;

    for (int j = 0; j < N_loop; j++)
    {

        lambda[j] = lambda_min + d_lambda * j;
        double k0 = 2. * PI / lambda[j];

        hantei[0][j] = ((double)n_hantei(k0, rad, V_per45) + (double)n_hantei(k0, rad, V_ge)) / 2.;
        hantei[1][j] = ((double)n_hantei(k0, rad, V_ni) + (double)n_hantei(k0, rad, V_ti)) / 2.;

        // std::cout << hantei[0][j] << " " <<hantei[1][j] << std::endl;
        /** 磁気ミラーの計算**/
        Eigen::Matrix2d M_per45 = mat_layer(k0, rad, V_per45, D_per45); // permalloy45の行列
        Eigen::Matrix2d M_ge = mat_layer(k0, rad, V_ge, D_ge);          // Geの行列

        Eigen::Matrix2d layermatrix_per_ge; // 二個分の行列
        Eigen::Matrix2d refmatrix_per_ge;   // 多層膜の行列

        layermatrix_per_ge = M_ge * M_per45;


        // std::cout << layermatrix.at(0).at(0) <<" "<< layermatrix.at(0).at(1) <<" "<< layermatrix.at(1).at(0) <<" "<< layermatrix.at(1).at(1) << std::endl;

        // まずは単位行列
        refmatrix_per_ge(0,0) = 1.;
        refmatrix_per_ge(0,1) = 0.;
        refmatrix_per_ge(1,0) = 0.;
        refmatrix_per_ge(1,1) = 1.;

        for (int l = 0; l < N_bilayer; l++)
        {
            refmatrix_per_ge = layermatrix_per_ge * refmatrix_per_ge;
        }

        double A = refmatrix_per_ge(0,0);
        double B = refmatrix_per_ge(0,1);
        double C = refmatrix_per_ge(1,0);
        double D = refmatrix_per_ge(1,1);

        r = reflect(k0, rad, V_sio2, A, B, C, D);
        R[0][j] = std::norm(r);

        // std::cout << R[1][j] << std::endl;

        /**全反射ミラーの計算**/
        Eigen::Matrix2d M_ni = mat_layer(k0, rad, V_ni, D_ni);  // Niの行列
        Eigen::Matrix2d M_ti = matrix(k0, rad, V_ti, D_ti); // Tiの行列

        Eigen::Matrix2d layermatrix_ni_ti; // 二個分の行列
        Eigen::Matrix2d refmatrix_ni_ti;   // 多層膜の行列

        layermatrix_ni_ti = M_ti * M_ni;

        // std::cout << layermatrix.at(0).at(0) <<" "<< layermatrix.at(0).at(1) <<" "<< layermatrix.at(1).at(0) <<" "<< layermatrix.at(1).at(1) << std::endl;

        // まずは単位行列
        refmatrix_ni_ti(0,0) = 1.;
        refmatrix_ni_ti(0,1) = 0.;
        refmatrix_ni_ti(1,0) = 0.;
        refmatrix_ni_ti(1,1) = 1.;

        for (int l = 0; l < N_bilayer; l++)
        {
            refmatrix_ni_ti = layermatrix_ni_ti * refmatrix_ni_ti;
        }

        double E = refmatrix_ni_ti(0,0);
        double F = refmatrix_ni_ti(0,1);
        double G = refmatrix_ni_ti(1,0);
        double H = refmatrix_ni_ti(1,1);

        r = reflect(k0, rad, V_sio2, E, F, G, H);
        R[1][j] = std::norm(r);

        // 磁気ミラーのdownに対しての反射率
        std::vector<std::vector<double>>
            M_per45_d = matrix(k0, rad, V_per45_d, D_per45);                              // permalloy45の行列
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

        for (int l = 0; l < N_bilayer; l++)
        {
            refmatrix_per_ge_d = DOT(layermatrix_per_ge_d, refmatrix_per_ge_d);
        }

        double I = refmatrix_per_ge_d.at(0).at(0);
        double J = refmatrix_per_ge_d.at(0).at(1);
        double K = refmatrix_per_ge_d.at(1).at(0);
        double L = refmatrix_per_ge_d.at(1).at(1);

        r = reflect(k0, rad, V_sio2, I, J, K, L);
        R[2][j] = std::norm(r);
    }

    TCanvas *c_refforlambda1 = new TCanvas("c_refforlambda1", "c_refforlambda1");
    TMultiGraph *mg1 = new TMultiGraph();
    TGraph *graph1 = new TGraph(N_loop, lambda, R[0]);
    TGraph *graph1_d = new TGraph(N_loop, lambda, R[2]);
    mg1->SetTitle(Form("Permalloy45/Ge at %f degree;lambda [m];reflectivity", theta));
    graph1->SetLineColor(1);
    graph1_d->SetLineColor(2);

    mg1->Add(graph1);
    mg1->Add(graph1_d);
    mg1->Draw("AL");

    TCanvas *c_refforlambda2 = new TCanvas("c_refforlambda2", "c_refforlambda2");

    TMultiGraph *mg2 = new TMultiGraph();
    TGraph *graph2 = new TGraph(N_loop, lambda, R[1]);
    TGraph *nhantei2 = new TGraph(N_loop, lambda, hantei[1]);
    graph2->SetTitle(Form("Ni/Ti at %f degree;lambda [m];reflectivity", theta));
    graph2->SetLineColor(1);
    graph2->Draw("AL");
}
