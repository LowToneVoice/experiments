#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#define J_per_eV 1.6022e-19
#define hbar 1.055e-34
#define m 1.675e-27
#define g 9.8
#define pi 3.14159265

constexpr double delta_D = 1; // ratio of layers thickness
// constexpr double delta_D = .93;

// optical potentials
constexpr double V_sio2 = 90.9e-9 * J_per_eV;
constexpr double V_ni = 224.e-9 * J_per_eV;
constexpr double V_ti = -40.e-9 * J_per_eV;

// constants about the multilayer mirror
constexpr double D_ni = 13.35e-9 / delta_D;
constexpr double D_ti = 9.83e-9 / delta_D;
constexpr int N_bilayer = 8;

// constants about alignments
constexpr double lambda_min = 2.e-10;
constexpr double lambda_max = 9.e-10;
constexpr double d_lambda = .01e-10;

// loop No.
constexpr int N_loop = (int)((lambda_max - lambda_min) / d_lambda);

// MATRIX OF EACH LAYER
std::vector<std::vector<double>> matrix(double k0, double theta, double V, double D)
{
    std::vector<std::vector<double>> M(2, std::vector<double>(2, 0.));
    double kx_0 = k0 * std::sin(theta);
    double kx_in_square = kx_0 * kx_0 - 2. * m * V / (hbar * hbar);

    // in case n is imaginary
    if (kx_in_square < 0)
    {
        double kx_in = sqrt(-kx_in_square);
        double nx = kx_in / kx_0;
        double zeta = kx_0 * D;

        M.at(0).at(0) = std::cosh(nx * zeta);
        M.at(0).at(1) = std::sinh(nx * zeta) / nx;
        M.at(1).at(0) = -std::sinh(nx * zeta) * nx;
        M.at(1).at(1) = std::cosh(nx * zeta);
    }
    else
    {
        double kx_in = sqrt(kx_in_square);
        double nx = kx_in / kx_0;
        double zeta = kx_0 * D;

        M.at(0).at(0) = std::cos(nx * zeta);
        M.at(0).at(1) = std::sin(nx * zeta) / nx;
        M.at(1).at(0) = -std::sin(nx * zeta) * nx;
        M.at(1).at(1) = std::cos(nx * zeta);
    }
    return M;
}
