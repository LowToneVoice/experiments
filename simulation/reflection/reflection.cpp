#include<iostream>
#include<fstream>
#include<Eigen/Dense>
#include<complex>

// FILES
#define OUTPUT_R "./reflection.dat"
#define OUTPUT_T "./transparency.dat"
#define OUTPUT_N "./norm.dat"

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

constexpr double delta_D = .93;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;

constexpr int N_bilayer = 8;

constexpr double lambda_min = 3.0e-10;
constexpr double lambda_max = 9.e-10;
constexpr double d_lambda = 0.01e-10;

constexpr double theta = 1.05 * pi / 180;

constexpr int N_loop = (int)((lambda_max - lambda_min) / d_lambda);

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

int reflection()
{
    double lambda[N_loop];
    double k0;
    double R2norm[N_loop];
    double T2norm[N_loop];
    double zeta;
    std::complex<double> R[N_loop], T[N_loop];
    std::complex<double> R_normed, T_normed;
    Eigen::Matrix2d layermatrix_ni_ti;
    Eigen::Matrix2d refmat;
    std::ofstream fileR(OUTPUT_R);
    std::ofstream fileT(OUTPUT_T);
    std::ofstream fileN(OUTPUT_N);

    for (int j = 0; j < N_loop; j++)
    {
        lambda[j] = lambda_min + d_lambda * j;
        k0 = 2. * pi / lambda[j];
        zeta = 0;

        layermatrix_ni_ti = mat_layer(k0, theta, V_ti, D_ti) * mat_layer(k0, theta, V_ni, D_ni);
        refmat = Eigen::Matrix2d::Identity();
        for (int l = 0; l < N_bilayer; l++)
        {
            refmat = layermatrix_ni_ti * refmat;
            zeta += D_ni + D_ti;
        }

        R[j] = reflect(k0, theta, V_sio2, refmat);
        T[j] = transparent(k0, theta, V_sio2, zeta, refmat);
        R2norm[j] = std::norm(R[j]);
        T2norm[j] = std::norm(T[j]);
        R_normed = R[j] / sqrt(std::norm(R[j]) + std::norm(T[j]));
        T_normed = T[j] / sqrt(std::norm(R[j]) + std::norm(T[j]));

        fileR << lambda[j] << " " << R[j].real() << " " << R[j].imag() << " " << R2norm[j] << std::endl;
        fileT << lambda[j] << " " << T[j].real() << " " << T[j].imag() << " " << T2norm[j] << std::endl;
        fileN << lambda[j] << " " << R2norm[j] + T2norm[j] << " " << std::norm(R_normed) + std::norm(T_normed) << std::endl;
    }

    fileR.close();
    fileT.close();
    fileN.close();
    return 0;
}

int main()
{
    return reflection();
}
