#include <iostream>
#include <fstream>
#include <random>
#include <complex>

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
constexpr double gap = 189e-6;

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
constexpr double theta_error = 1e-3 * pi / 180;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1.;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
// constexpr double d_lambda = 0.00001e-10;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;

constexpr int beam_count = 1000;
constexpr int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);
constexpr int N_loop_theta = (int)((theta_max - theta_min) / d_theta);

int main()
{
    double Phi_g_main_perLambda = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta);
    double Phi_a_main_toLambda = 4 * pi * gap * theta_error;
    double Phi_g_sub_perLambda = 2 * pi * g * pow(m / h, 2) * theta_error / 2 * pow(gap / sin(theta), 2);
    double Phi_a_sub_perLambda = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * theta_error;

    std::cout << "G main: " << Phi_g_main_perLambda << " * lambda" << std::endl;
    std::cout << "G sub: " << Phi_g_sub_perLambda << " * lambda" << std::endl;
    std::cout << "A main: " << Phi_a_main_toLambda << " / lambda" << std::endl;
    std::cout << "A sub: " << Phi_a_sub_perLambda << " * lambda" << std::endl;
}
