#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <TFile.h>
#include <TTree.h>
#include <Eigen/Dense>

#define MAX_DATA_SIZE 1e6
#define OUTPUT_TREE "./dat/RT.root"
#define OUTPUT_O "./dat/transparent.dat"
#define OUTPUT_H "./dat/reflection.dat"

// PHYS & MATH CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define pi M_PI
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
#define J_per_eV 1.6022e-19
constexpr std::complex<double> i = std::complex<double>(0, 1);
constexpr double V_sio2 = 90.9e-09 * J_per_eV;
constexpr double V_ni = 224.e-09 * J_per_eV;
constexpr double V_ti = -40.e-9 * J_per_eV;
constexpr double rho_Ni = 9.1e28;
constexpr double bc_Ni = 3.4e-15;

// CHANNELS
#define H_TDC 0
#define O_TDC 1
#define H_ADC 2
#define O_ADC 3

// ALIGNMENTS CONSTANTS: !!!!! DO NOT CHANGE !!!!!
constexpr double daq_freq = 62.5e6;
constexpr double gap = 189e-6;
constexpr double LMD_MIN = 2e-10;
constexpr double LMD_MAX = 9e-10;
constexpr double beam_intensity = 9.4e7; // /s/cm2
constexpr double beam_size = 5.5 * 4.5;  // cm2
constexpr double delta_D = 1;
// constexpr double delta_D = 0.93;
constexpr double D_SiO2 = 14e-3;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;
constexpr int N_bilayer = 8;

// ALIGNMENT VARIABLES
constexpr int daq_downsizing = 1;
constexpr double theta = 1.05 * pi / 180;
constexpr double theta0 = 1.e-3 * pi / 180;
constexpr double angle_delta = pi / 2;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1;
constexpr double beam_time = 30 * 60; // sec
// constexpr int beam_count = beam_time * BEAM_SIZE;
constexpr int beam_count = 10000;
constexpr int particle_count = 100000;
constexpr double dt = 1 / daq_freq * daq_downsizing;
constexpr double dLMD = h / total_length / m * dt;
// constexpr int L_max = (int)((LMD_MAX - LMD_MIN) / dLMD);
constexpr int L_max = 10000;

// select lambda following the given probability
double lambda_selection()
{
    // random var.
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> rand(0.0, 1.0);

    // Approximate the wavelength distribution based on the two points (700, 0.388) and (870, 0.16625)
    double p = rand(mt);
    double alpha = (0.388 - 0.16625) / (870e-12 - 700e-12);
    double L_min = 700e-12, L_max = 878.5e-12;
    double y0 = .388;
    double S = .5 * (y0 - alpha * (L_max - L_min) + y0) * (L_max - L_min);
    return L_min + (y0 - sqrt(y0 * y0 - 2 * alpha * p * S)) / alpha;
}

// matrix M_j of a layer
Eigen::Matrix<std::complex<double>, 2, 2> mat_layer(double k0_vertical, double thickness, double V)
{
    double kj_square = pow(k0_vertical, 2) - 2 * m * V / pow(h / 2 / pi, 2);
    double n, zeta;
    Eigen::Matrix<std::complex<double>, 2, 2> matrix;

    // in case kj in real
    if (kj_square >= 0)
    {
        n = sqrt(kj_square) / k0_vertical;
        matrix << std::complex<double>(cos(n * thickness), 0), std::complex<double>(sin(n * thickness) / n, 0), std::complex<double>(-sin(n * thickness) * n, 0), std::complex<double>(cos(n * thickness), 0);
    }
    else
    {
        n = sqrt(-kj_square) / k0_vertical;
        matrix << std::complex<double>(cosh(n * thickness), 0), std::complex<double>(sinh(n * thickness) / n, 0), std::complex<double>(-sinh(n * thickness) * n, 0), std::complex<double>(cosh(n * thickness), 0);
    }

    return matrix;
}

// complex reflection factor
std::vector<std::complex<double>> refl_trans_factor_complex(double lambda, double theta)
{
    std::vector<std::complex<double>> vector;
    Eigen::Matrix<std::complex<double>, 2, 2> matrix, mat_SiO2, mat_Ni, mat_Ti;
    double k0_vertical = 2 * pi / lambda * sin(theta0);
    std::complex<double> R, T;
    double zeta = 0;

    matrix << mat_layer(k0_vertical, D_SiO2, V_sio2);
    mat_Ni << mat_layer(k0_vertical, D_ni, V_ni);
    mat_Ti << mat_layer(k0_vertical, D_ti, V_ti);
    zeta += D_SiO2;

    for (int j = 0; j < N_bilayer; j++)
    {
        matrix = mat_Ni * matrix;
        zeta += D_ni;
        matrix = mat_Ti * matrix;
        zeta += D_ti;
    }

    std::complex<double> M11 = matrix(0, 0);
    std::complex<double> M12 = matrix(0, 1);
    std::complex<double> M21 = matrix(1, 0);
    std::complex<double> M22 = matrix(1, 1);
    std::complex<double> ik = i * k0_vertical;
    R = ((M21 + ik * M22) - ik * (M22 - M11)) / ((pow(k0_vertical, 2) * M12 - M21) + ik * (M22 + M11));
    T = 2 * (cos(k0_vertical * zeta), sin(k0_vertical * zeta)) * ik * (M11 * M22 - M12 * M21) / ((pow(k0_vertical, 2) * M12 - M21) + ik * (M22 + M11));
    vector[0] = R;
    vector[1] = T;

    return vector;
}

// TODO: #1 building alpha_calc in sim copy.cpp
// double alpha_calc(double lambda)
// {
//     return 0.5;
// }

int sim_copy()
{
    // PREPARATIONS
    std::cout << "///// PREPARATIONS /////" << std::endl;

    // variables
    double alpha = .5, gamma = .5;
    int count[2] = {};
    int j;
    double lmd;
    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double r;
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> rand(0.0, 1.0);

    // open files
    std::ofstream file_O(OUTPUT_O);
    std::ofstream file_H(OUTPUT_H);
    if (!file_O.is_open() || !file_H.is_open())
    {
        std::cerr << "Failed to open the input files." << std::endl;
        return 1;
    }
    else
    {
        std::cout << "File has safely opened." << std::endl;
    }

    // SIMULATION
    std::cout << "///// SIMULATIONS /////" << std::endl;

    // root file
    Double_t lambda;
    Int_t channel;
    TFile *file = new TFile(OUTPUT_TREE, "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("channel", &channel);
    tree->Branch("lambda", &lambda);
    double probO, probH;
    std::complex<double> R, T;
    std::vector<std::complex<double>> vec;

    // loop
    for (j = 0; j < particle_count; j++)
    {
        lambda = lambda_selection();
        vec = refl_trans_factor_complex(lambda, theta);
        R = vec[0];
        T = vec[1];
        // alpha = alpha_calc(lambda);
        // gamma = 1 - alpha;

        // phase calculation
        Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance * sin(angle_delta) / tan(2 * theta) * lambda;
        Phi_a_main = 4 * pi * gap / lambda * theta0;
        Phi_g_sub = 2 * pi * g * pow(m / h, 2) * theta0 / 2 * pow(gap / sin(theta), 2) * lambda;
        Phi_a_sub = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * lambda * theta0;
        Phase = Phi_g_main + Phi_g_sub;

        // probability calculation
        probO = std::norm(T);
        probH = std::norm(R);
        r = rand(mt);
        if (r <= probO)
        {
            file_O << lambda << std::endl;
            channel = O_TDC;
            tree->Fill();
            count[0]++;
        }
        else if (r >= 1 - probH)
        {
            file_H << lambda << std::endl;
            channel = H_TDC;
            tree->Fill();
            count[1]++;
        }

        if (j % 1000 == 999)
        {
            std::cout << "loop " << j + 1 << " / " << particle_count << " has ended." << std::endl;
        }
        if (count[0] + count[1] > MAX_DATA_SIZE)
        {
            std::cout << "EXCEEDED MAX DATASIZE" << std::endl;
            break;
        }
    }

    std::cout << "O hit: " << count[0] << std::endl;
    std::cout << "H hit: " << count[1] << std::endl;

    // ENDING
    tree->Write();
    file_O.close();
    file_H.close();

    return 0;
}

int main()
{
    return sim_copy();
}
