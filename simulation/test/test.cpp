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
constexpr double hbar = h / 2 / pi;
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
Eigen::Matrix2d mat_layer(double k0_vertical, double thickness, double V)
{
    double kj_square = pow(k0_vertical, 2) - 2 * m * V / (hbar * hbar);
    double n, zeta, nd;
    Eigen::Matrix2d matrix;

    // in case kj in real
    if (kj_square >= 0)
    {
        n = sqrt(kj_square) / k0_vertical;
        nd = n * thickness;
        matrix(0, 0) = cos(nd);
        matrix(0, 1) = sin(nd) / n;
        matrix(1, 0) = -sin(nd) * n;
        matrix(1, 1) = cos(nd);
    }
    else
    {
        n = sqrt(-kj_square) / k0_vertical;
        nd = n * thickness;
        matrix(0, 0) = cosh(nd);
        matrix(0, 1) = sinh(nd) / n;
        matrix(1, 0) = n * sinh(nd);
        matrix(1, 1) = cosh(nd);
    }

    return matrix;
}

// complex reflection and transparent factor
std::complex<double> refl_trans_factor_substance_changed(double lambda, double theta, double Vg, double A, double B, double C, double D)
{
    double k0_vertical = 2 * pi / lambda * sin(theta);
    double kg_sq = k0_vertical * k0_vertical - 2. * m * Vg / (hbar * hbar);
    std::complex<double> R, denom;

    if (kg_sq >= 0.)
    {
        double kg = sqrt(kg_sq);
        double n0 = 1.;
        double ng = kg / k0_vertical;
        R = std::complex<double> (n0 * ng * B + C, n0 * D - ng * A);
        denom = std::complex<double> (n0 * ng * B - C, ng * A + n0 * D);
        R /= denom;
    }
    else
    {
        double kg = sqrt(-kg_sq);
        double n0 = 1.;
        double ng = kg / k0_vertical;
        R = std::complex<double> (ng * A + C, n0 * D + n0 * ng * B);
        denom = std::complex<double> (-ng * A - D, n0 * ng * B + n0 * C);
        R /= denom;
    }

    return R;
}

// complex reflection factor
Eigen::Vector2cd refl_trans_factor_complex(double lambda, double theta, bool Si_layer_up, bool Si_layer_down)
{
    Eigen::Vector2cd vector;
    Eigen::Matrix2d matrix, mat_SiO2, mat_Ni, mat_Ti;
    double k0_vertical = 2. * pi / lambda * sin(theta);
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
        matrix = Eigen::Matrix2d::Identity();
    }

    mat_Ni = mat_layer(k0_vertical, D_ni, V_ni);
    mat_Ti = mat_layer(k0_vertical, D_ti, V_ti);

    for (int j = 0; j < N_bilayer; j++)
    {
        matrix = mat_Ni * matrix;
        zeta += D_ni;
        matrix = mat_Ti * matrix;
        zeta += D_ti;
    }

    // back SiO2 layer
    if (Si_layer_down)
    {
        matrix = mat_layer(k0_vertical, D_SiO2, V_sio2) * matrix;
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
    double r;
    double maxR = 0, minR = 1, maxT = 0, minT = 1;

    std::ofstream file_reflect(OUTPUT_R);
    std::ofstream file_transparent(OUTPUT_T);
    if(!file_reflect.is_open()||!file_transparent.is_open())
    {
        std::cerr << "FILE FAILED." << std::endl;
        return 1;
    }

    std::complex<double> R, T;
    Eigen::Matrix2d M, M_ni, M_ti;
    Eigen::Vector2cd vec;
    double k0_v;

    for (lambda = LMD_MIN; lambda < LMD_MAX; lambda += dLMD)
    {
        // vec = refl_trans_factor_complex(lambda, theta0, false, true);
        // R = vec[0];
        // T = vec[1];
        M = Eigen::Matrix2d::Identity();
        k0_v = 2 * pi / lambda * sin(theta);
        M_ni = mat_layer(k0_v, D_ni, V_ni);
        M_ti = mat_layer(k0_v, D_ti, V_ti);
        for ( j = 0; j < N_bilayer; j++)
        {
            M = M_ti * M_ni * M;
        }

        R = refl_trans_factor_substance_changed(lambda, theta, V_sio2, M(0, 0), M(0, 1), M(1, 0), M(1, 1));

        file_reflect
            << lambda << " " << R.real() << " " << R.imag() << " " << std::norm(R) << std::endl;
        file_transparent << lambda << " " << T.real() << " " << T.imag() << " " << std::norm(T) << std::endl;
    }

    file_reflect.close();
    file_transparent.close();

    return 0;
}
