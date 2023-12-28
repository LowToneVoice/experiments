#include <iostream>
#include <fstream>
#include <random>
#include <complex>
#include <Eigen/Dense>

#define OUTPUT "prob.dat"

// PHYS & MATH CONSTANTS
std::complex<double> I(0, 1.);
#define h 6.62607015e-34
#define pi M_PI
#define J_per_eV 1.6022e-19
#define m 1.675e-27
#define NA 6.02214e23
constexpr double g = 9.8;
constexpr double hbar = h / (2 * pi);

constexpr double V_sio2 = 90.9e-09 * J_per_eV;
constexpr double V_ni = 224.e-09 * J_per_eV;
constexpr double V_ti = -40.e-9 * J_per_eV;
constexpr double m_Ni = 58.6934e-3 / NA;
constexpr double m_Ti = 47.867e-3 / NA;
constexpr double rho_Si = 5.00e22; // used cm3 unit in Seki's paper
constexpr double rho_Ni = 8.908e-3 / m_Ni;
constexpr double rho_Ti = 4.506e-3 / m_Ti;
constexpr double bc_Ni = 10.3e-15;
constexpr double bc_Ti = -3.36e-15;

// ALIGNMENT CONSTANTS
constexpr double delta_D = 1.;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;
constexpr double D_ni_error = .05e-10 / delta_D;
constexpr double D_ti_error = .03e-10 / delta_D;
constexpr double gap = 189e-6;
constexpr double gap_error = .1e-6;

constexpr int N_bilayer = 8;

constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 10e-10;
constexpr double theta_min = .2 * pi / 180;
constexpr double theta_max = 1.5 * pi / 180;

constexpr int daq_freq = 62.5e6;

// ALIGNMENT VARIABLES
constexpr double lambda_min_used = lambda_min;
constexpr double lambda_max_used = lambda_max;
constexpr double theta = 1.05 * pi / 180;
constexpr double theta_error = .001 * pi / 180;
constexpr double angle_from_parallel = 1e-3 * pi / 180;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1.;
constexpr double total_length_error = .005;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
// constexpr double d_lambda = 0.00001e-10;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;

constexpr int beam_count = 1000;
constexpr int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);
constexpr int N_loop_theta = (int)((theta_max - theta_min) / d_theta);

// matrix M of one layer
Eigen::Matrix2cd mat_layer(double k0, double theta, double V, double D)
{
    Eigen::Matrix2cd M;
    double kx0 = k0 * std::sin(theta);
    std::complex<double> kx = std::sqrt(kx0 * kx0 - 2. * m * V / (hbar * hbar));
    std::complex<double> nx = kx / kx0;
    double delta = kx0 * D;

    M(0, 0) = std::cos(nx * delta);
    M(0, 1) = std::sin(nx * delta) / nx;
    M(1, 0) = -std::sin(nx * delta) * nx;
    M(1, 1) = std::cos(nx * delta);

    return M;
}
Eigen::Matrix2cd mat_layer_error(double k0, double theta, double V, double D, double D_error)
{
    Eigen::Matrix2cd M;
    double kx0 = k0 * std::sin(theta);
    std::complex<double> kx = std::sqrt(kx0 * kx0 - 2. * m * V / (hbar * hbar));
    std::complex<double> nx = kx / kx0;
    double delta = kx0 * D;
    double delta_error = kx0 * D_error;

    M(0, 0) = -nx * delta_error * std::sin(nx * delta);
    M(0, 1) = delta_error * std::cos(nx * delta);
    M(1, 0) = -nx * nx * delta_error * std::cos(nx * delta);
    M(1, 1) = -nx * delta_error * std::sin(nx * delta);

    return M;
}

// complex reflection factor
std::complex<double> reflect(double k0, double theta, double Vg, Eigen::Matrix2cd M)
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

// complex transparent factor
std::complex<double> transparent(double k0, double theta, double Vg, double zeta, Eigen::Matrix2cd M)
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

int main()
{
    double Phi_g_main_perLambda = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta);
    double Phi_a_main_toLambda = 4 * pi * gap * angle_from_parallel;
    double Phi_g_sub_perLambda = 2 * pi * g * pow(m / h, 2) * angle_from_parallel / 2 * pow(gap / sin(theta), 2);
    double Phi_a_sub_perLambda = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * angle_from_parallel;
    double delta_g_main_perLambda = abs(Phi_g_main_perLambda * (gap_error / gap + total_length_error / total_length + 2 * theta_error / pow(theta, 2)));
    double delta_g_sub_perLambda = abs(Phi_g_sub_perLambda * (2 * gap_error / gap + 2 * theta_error / theta));
    double delta_a_main_toLambda = abs(Phi_a_main_toLambda * gap_error / gap);
    double delta_a_sub_perLambda = abs(Phi_a_sub_perLambda * (gap_error / gap + 2 * theta_error / theta));

    std::cout << "G main: (" << Phi_g_main_perLambda << " +- " << delta_g_main_perLambda << ") * lambda " << std::endl;
    std::cout << "G sub: (" << Phi_g_sub_perLambda << " +- " << delta_g_sub_perLambda << ") * lambda" << std::endl;
    std::cout << "G main: (" << Phi_a_main_toLambda << " +- " << delta_a_main_toLambda << ") / lambda " << std::endl;
    std::cout << "G sub: (" << Phi_a_sub_perLambda << " +- " << delta_a_sub_perLambda << ") * lambda" << std::endl;

    // matrix error
    double lambda, k0, zeta;
    Eigen::Matrix2cd mat_bilayer_ni_ti, M, mat_bilayer_ni_ti_error, M_mid;
    std::complex<double> R_not_normed, T_not_normed, R, T;
    double R2norm, T2norm, probO, probH, Phase;
    std::ofstream output(OUTPUT);

    for (int j = 0; j < N_loop_lambda; j++)
    {
        // initialization
        lambda = lambda_min_used + d_lambda * j;
        k0 = 2. * pi / lambda;
        zeta = 0;

        // reflection calculation without error
        mat_bilayer_ni_ti = mat_layer(k0, theta, V_ti, D_ti) * mat_layer(k0, theta, V_ni, D_ni);
        M = Eigen::Matrix2cd::Identity();
        for (int l = 0; l < N_bilayer; l++)
        {
            M = mat_bilayer_ni_ti * M;
            zeta += D_ni + D_ti;
        }
        R_not_normed = reflect(k0, theta, V_sio2, M);
        T_not_normed = transparent(k0, theta, V_sio2, zeta, M);
        R = R_not_normed / sqrt(std::norm(R_not_normed) + std::norm(T_not_normed));
        T = T_not_normed / sqrt(std::norm(R_not_normed) + std::norm(T_not_normed));
        R2norm = std::norm(R);
        T2norm = std::norm(T);

        Phase = Phi_g_main_perLambda * lambda * (1 + abs(Phi_g_sub_perLambda / Phi_g_main_perLambda));

        // probability calculation
        probO = std::norm(T * R * T * R + R * T * R * T * std::exp(I * Phase));
        probH = std::norm(T * R * T * T * T + R * T * R * R * T * std::exp(I * Phase) + R * T * T * std::exp(I * Phase));
        output << lambda << " " << probO << " " << probH;


        // reflection calculation with error
        mat_bilayer_ni_ti = (mat_layer(k0, theta, V_ti, D_ti) + mat_layer_error(k0, theta, V_ti, D_ti, D_ti_error)) * (mat_layer(k0, theta, V_ni, D_ni) + mat_layer_error(k0, theta, V_ni, D_ni, D_ni_error));
        M = Eigen::Matrix2cd::Identity();
        for (int l = 0; l < N_bilayer; l++)
        {
            M = mat_bilayer_ni_ti * M;
            zeta += D_ni + D_ti + D_ni_error + D_ti_error;
        }
        R_not_normed = reflect(k0, theta, V_sio2, M);
        T_not_normed = transparent(k0, theta, V_sio2, zeta, M);
        R = R_not_normed / sqrt(std::norm(R_not_normed) + std::norm(T_not_normed));
        T = T_not_normed / sqrt(std::norm(R_not_normed) + std::norm(T_not_normed));
        R2norm = std::norm(R);
        T2norm = std::norm(T);

        Phase = Phi_g_main_perLambda * lambda * (1 + abs(Phi_g_sub_perLambda / Phi_g_main_perLambda));

        // probability calculation
        probO = std::norm(T * R * T * R + R * T * R * T * std::exp(I * Phase));
        probH = std::norm(T * R * T * T * T + R * T * R * R * T * std::exp(I * Phase) + R * T * T * std::exp(I * Phase));
        output << " " << probO << " " << probH << std::endl;

    }

    output.close();

    return 0;
}
