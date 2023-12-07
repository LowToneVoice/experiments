#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <Eigen/Dense>

#define MAX_DATA_SIZE 1e6
#define OUTPUT_TREE "./test.root"
#define OUTPUT_T "./test_trans.dat"
#define OUTPUT_R "./test_refl.dat"
#define DISPLAY_INTERVAL 1000

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

// matrix M_j of a layer
Eigen::Matrix2cd mat_layer(double k0_vertical, double thickness, double V)
{
    double kj_square = pow(k0_vertical, 2) - 2 * m * V / pow(h / 2 / pi, 2);
    double n, zeta, nd;
    Eigen::Matrix2cd matrix;

    // in case kj in real
    if (kj_square >= 0)
    {
        n = sqrt(kj_square) / k0_vertical;
        nd = n * thickness;
        matrix(0, 0) = std::complex<double>(cos(nd), 0);
        matrix(0, 1) = std::complex<double>(sin(nd) / n, 0);
        matrix(1, 0) = std::complex<double>(-n * sin(nd), 0);
        matrix(1, 1) = std::complex<double>(cos(nd), 0);
    }
    else
    {
        n = sqrt(-kj_square) / k0_vertical;
        nd = n * thickness;
        matrix(0, 0) = std::complex<double>(cosh(nd), 0);
        matrix(0, 1) = std::complex<double>(sinh(nd) / n, 0);
        matrix(1, 0) = std::complex<double>(-n * sinh(nd), 0);
        matrix(1, 1) = std::complex<double>(cosh(nd), 0);
    }

    return matrix;
}

// complex reflection factor
Eigen::Vector2cd refl_trans_factor_complex(double lambda, double theta, int Si_layer_up)
{
    Eigen::Vector2cd vector;
    Eigen::Matrix2cd matrix, mat_SiO2, mat_Ni, mat_Ti;
    double k0_vertical = 2 * pi / lambda * sin(theta0);
    std::complex<double> R, T;
    double zeta = 0;

    // front SiO2 layer
    if (Si_layer_up)
    {
        matrix << mat_layer(k0_vertical, D_SiO2, V_sio2);
        zeta += D_SiO2;
    }
    else
    {
        matrix << std::complex<double>(1, 0), std::complex<double>(0, 0), std::complex<double>(1, 0), std::complex<double>(0, 0);
    }

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

    // back SiO2 layer
    if (!Si_layer_up)
    {
        matrix << mat_layer(k0_vertical, D_SiO2, V_sio2);
        zeta += D_SiO2;
    }

    std::complex<double> M11 = matrix(0, 0);
    std::complex<double> M12 = matrix(0, 1);
    std::complex<double> M21 = matrix(1, 0);
    std::complex<double> M22 = matrix(1, 1);
    std::complex<double> ik = i * k0_vertical;
    std::complex<double> k2 = k0_vertical * k0_vertical;
    double kz = k0_vertical * zeta;

    R = ((M21 + ik * M22) - ik * (M11 + ik * M12)) / (ik * (M11 - ik * M12) - (M21 - ik * M22));
    T = 2. * ik * std::complex<double>(cos(kz * zeta), -sin(kz * zeta)) / (ik * (M11 - ik * M12) - (M21 - ik * M22));
    vector[0] = R;
    vector[1] = T;

    return vector;
}

int main()
{
    double lambda;

    int count[2] = {};
    int j;
    double lmd;
    double r;
    double maxR = 0, minR = 1, maxT = 0, minT = 1;

    std::ofstream file_reflect(OUTPUT_R);
    std::ofstream file_transparent(OUTPUT_T);
    if(!file_reflect.is_open()||!file_transparent.is_open())
    {
        std::cerr << "FILE FAILED." << std::endl;
        return 1;
    }

    double probR, probT;
    std::complex<double> R, T;
    Eigen::Vector2cd vec;
    for (lambda = LMD_MIN; lambda < LMD_MAX; lambda += dLMD)
    {
        vec = refl_trans_factor_complex(lambda, theta0, 0);
        R = vec[0];
        T = vec[1];
        file_reflect << lambda << " " << R.real() << " " << R.imag() <<" "<< std::norm(R) << std::endl;
        file_transparent << lambda << " " << T.real() << " " << T.imag() << " " << std::norm(T) << std::endl;
    }

    file_reflect.close();
    file_transparent.close();

    return 0;
}
