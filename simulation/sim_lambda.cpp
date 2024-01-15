#include <iostream>
#include <fstream>
#include <random>
#include <Eigen/Dense>
#include <complex>
#include <TFile.h>
#include <TTree.h>

// FILES
#define OUTPUT_TREE "./dat/montecarlo/ref/lambda/g_mix_5deg.root"
#define OUTPUT_O "./dat/montecarlo/ref/lambda/g_mix_5deg_O.dat"
#define OUTPUT_H "./dat/montecarlo/ref/lambda/g_mix_5deg_H.dat"
#define OUTPUT_PROB "./dat/theoretical/ref/lambda/g_mix_5deg.dat"

#define OUT_INTERVAL 1000

// CHANNELS
#define H_TDC 0
#define O_TDC 1
#define H_ADC 2
#define O_ADC 3

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
constexpr double rho_Si = 5.00e22;  // used cm3 unit in Seki's paper
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

constexpr double lambda_fade_width = .01 * (lambda_max - lambda_min);
// constexpr double lambda_fade_width = 1e-12;

// ALIGNMENT VARIABLES
constexpr double lambda_min_used = lambda_min;
constexpr double lambda_max_used = lambda_max;
constexpr double theta = 1.05 * pi / 180;
constexpr double delta = 5 * pi / 180;
constexpr double angle_from_parallel = .1 * pi / 180;
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
constexpr int N_loop_lambda_fade = (int)(lambda_fade_width / d_lambda);

// matrix M of one layer
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

// complex reflection factor
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

// complex transparent factor
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

int sim_lambda()
{
    double lambda;
    double k0;
    double zeta;
    std::complex<double> R, T;
    std::complex<double> R_not_normed, T_not_nnormed;
    std::complex<double> wave_fncO, wave_fncH;
    double R2norm, T2norm, R_phase, T_phase;
    Eigen::Matrix2d mat_bilayer_ni_ti;
    Eigen::Matrix2d M;
    std::ofstream fileO(OUTPUT_O);
    std::ofstream fileH(OUTPUT_H);
    std::ofstream fileP(OUTPUT_PROB);
    if (!fileO.is_open() || !fileH.is_open() || !fileP.is_open())
    {
        std::cerr << "FILE DID NOT OPENED" << std::endl;
        return 1;
    }

    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double probability, probO, probH;
    std::mt19937 mt{(std::random_device{}())};
    std::uniform_real_distribution<double> rand(0., 1.);
    Double_t lmd;
    Int_t channel;
    TFile file(OUTPUT_TREE, "recreate");
    TTree tree("tree", "tree");
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lmd);
    int count[2] = {};
    int l;

    for (int j = 0; j < N_loop_lambda; j++)
    {
        lmd = lambda_min_used + d_lambda * j;
        wave_fncH = 0;
        wave_fncO = 0;
        l = 0;

        do
        {
            // initialization
            l++;
            lambda = lmd - lambda_fade_width / 2 + l * d_lambda;
            k0 = 2. * pi / lambda;
            zeta = 0;

            // reflection calculation
            mat_bilayer_ni_ti = mat_layer(k0, theta, V_ti, D_ti) * mat_layer(k0, theta, V_ni, D_ni);
            M = Eigen::Matrix2d::Identity();
            for (int l = 0; l < N_bilayer; l++)
            {
                M = mat_bilayer_ni_ti * M;
                zeta += D_ni + D_ti;
            }
            R_not_normed = reflect(k0, theta, V_sio2, M);
            T_not_nnormed = transparent(k0, theta, V_sio2, zeta, M);
            R = R_not_normed / sqrt(std::norm(R_not_normed) + std::norm(T_not_nnormed));
            T = T_not_nnormed / sqrt(std::norm(R_not_normed) + std::norm(T_not_nnormed));
            R2norm = std::norm(R);
            R_phase = std::arg(R);
            T2norm = std::norm(T);
            T_phase = std::arg(T);

            // Phase calculation
            Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * lambda * sin(delta);
            Phi_a_main = 4 * pi * gap / lambda * angle_from_parallel;
            Phi_g_sub = 2 * pi * g * pow(m / h, 2) * angle_from_parallel / 2 * pow(gap / sin(theta), 2) * lambda * sin(delta);
            Phi_a_sub = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * lambda * angle_from_parallel; // maximize the value
            Phase = Phi_g_main + Phi_g_sub;

            // probability calculation
            // R = 1. / sqrt(2);
            // T = I / sqrt(2);
            wave_fncH += (T * R * T * T * T + R * T * R * R * T * std::exp(I * Phase) + R * T * T * std::exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncO += (T * R * T * R + R * T * R * T * std::exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
        } while (l < N_loop_lambda_fade);

        probO = std::norm(wave_fncO);
        probH = std::norm(wave_fncH);

        fileP << lambda << " " << probO << " " << probH << std::endl;

        for (int beam = 0; beam < beam_count; beam++)
        {
            probability = rand(mt);
            if (probability < probO)
            {
                fileO << lambda << std::endl;
                channel = O_TDC;
                tree.Fill();
                count[0]++;
            }
            else if (1 - probability <= probH)
            {
                fileH << lambda << std::endl;
                channel = H_TDC;
                tree.Fill();
                count[1]++;
            }
        }
    }

    tree.Write();
    fileO.close();
    fileH.close();
    fileP.close();
    file.Close();

    std::cout << "counts: " << count[0] << " and " << count[1] << std::endl;
    std::cout << "Phi_g_main = " << -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * sin(delta) << " * lambda" << std::endl;

    return 0;
}

int main()
{
    return sim_lambda();
}
