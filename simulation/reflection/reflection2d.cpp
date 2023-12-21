#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>

// FILES
#define OUTPUT_FILE "./RT2d.dat"

// PHYS & MATH CONSTANTS
std::complex<double> I(0, 1.);
#define h 6.62607015e-34
#define pi 3.14159265
#define J_per_eV 1.6022e-19
constexpr double hbar = h / (2 * pi);
constexpr double m = 1.675e-27;

constexpr double V_sio2 = 90.9e-09 * J_per_eV;
constexpr double V_ni = 224.e-09 * J_per_eV;
constexpr double V_ti = -40.e-9 * J_per_eV;

constexpr double delta_D = 1.;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;

constexpr int N_bilayer = 8;

constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 12e-10;
constexpr double d_lambda = 0.01e-10;

constexpr double theta_min = 0 * pi / 180;
constexpr double theta_max = 1.5 * pi / 180;
constexpr double d_theta = .01 * pi / 180;

constexpr int N_loop_lambda = (int)((lambda_max - lambda_min) / d_lambda);
constexpr int N_loop_theta = (int)((theta_max - theta_min) / d_theta);

Eigen::Matrix2d mat_layer(double k0, double theta, double V, double D)
{
    Eigen::Matrix2d M;
    double kx0 = k0 * std::sin(theta);
    std::complex<double> kx = std::sqrt(kx0 * kx0 - 2. * m * V / (hbar * hbar));
    std::complex<double> nx = kx / kx0;
    double delta = kx0 * D;

    M(0, 0) = std::cos(nx * delta).real();
    M(0, 1) = (std::sin(nx * delta) / nx).real();
    M(1, 0) = (-std::sin(nx * delta) * nx).real();
    M(1, 1) = std::cos(nx * delta).real();

    return M;
}

std::complex<double> reflect(double k0, double theta, double Vg, Eigen::Matrix2d M)
{
    double kx0 = k0 * std::sin(theta);
    std::complex<double> kxg = std::sqrt(kx0 * kx0 - 2. * m * Vg / (hbar * hbar));
    double n0 = 1.;
    std::complex<double> ng = kxg / kx0;
    std::complex<double> numer, denom;

    numer = -n0 * ng * M(0, 1) - M(1, 0) + I * (ng * M(0, 0) - n0 * M(1, 1));
    denom = -n0 * ng * M(0, 1) + M(1, 0) + I * (-ng * M(0, 0) - n0 * M(1, 1));

    return numer / denom;
}

std::complex<double> transparent(double k0, double theta, double Vg, double zeta, Eigen::Matrix2d M)
{
    double kx0 = k0 * std::sin(theta);
    std::complex<double> kxg = std::sqrt(kx0 * kx0 - 2. * m * Vg / (hbar * hbar));
    double n0 = 1.;
    std::complex<double> ng = kxg / kx0;

    std::complex<double> eikgz = std::cos(kxg * zeta) + I * std::sin(kxg * zeta);
    std::complex<double> numer, denom;
    numer = 2. * I * n0 / eikgz;
    denom = n0 * ng * M(0, 1) - M(1, 0) + I * (n0 * M(1, 1) + ng * M(0, 0));

    return numer / denom;
}

int reflection2d()
{
    double lambda;
    double theta;
    double k0;
    double R2norm;
    double T2norm;
    double zeta;
    std::complex<double> R, T;
    std::complex<double> R_not_normed, T_not_normed;
    Eigen::Matrix2d bilayer_mat_ni_ti;
    Eigen::Matrix2d M;
    std::ofstream output(OUTPUT_FILE);

    for (int j = 0; j < N_loop_lambda; j++)
    {
        lambda = lambda_min + j * d_lambda;
        k0 = 2. * pi / lambda;

        for (int l = 0; l < N_loop_theta; l++)
        {
            theta = theta_min + l * d_theta;
            zeta = 0.;

            bilayer_mat_ni_ti = mat_layer(k0, theta, V_ti, D_ti) * mat_layer(k0, theta, V_ni, D_ni);
            M = Eigen::Matrix2d::Identity();
            for (int l = 0; l < N_bilayer; l++)
            {
                M = bilayer_mat_ni_ti * M;
                zeta += D_ni + D_ti;
            }

            R_not_normed = reflect(k0, theta, V_sio2, M);
            T_not_normed = transparent(k0, theta, V_sio2, zeta, M);
            R2norm = std::norm(R_not_normed);
            T2norm = std::norm(T_not_normed);
            R = R_not_normed / (sqrt(R2norm) + sqrt(T2norm));
            T = T_not_normed / (sqrt(R2norm) + sqrt(T2norm));

            output << lambda << " " << theta << " " << R.real() << " " << R.imag() << " " << T.real() << " " << T.imag() << " " << R2norm << " " << T2norm << std::endl;
        }

        output << std::endl;
    }

    output.close();
    return 0;
}

int main()
{
    return reflection2d();
}
