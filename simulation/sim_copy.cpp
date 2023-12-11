#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <random>
#include <TFile.h>
#include <TTree.h>
// #include <Eigen/Dense>

#define MAX_DATA_SIZE 1e6
#define OUTPUT_TREE "./dat/RT.root"
#define OUTPUT_O "./dat/transparent.dat"
#define OUTPUT_H "./dat/reflection.dat"
#define DISPLAY_INTERVAL 1000

// PHYS & MATH CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define pi M_PI
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
#define J_per_eV 1.6022e-19
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
// constexpr int beam_count = beam_time * BEAM_SIZE;
constexpr int beam_count = 10000;
constexpr int particle_count = 100000;
constexpr double dt = 1 / daq_freq * daq_downsizing;
constexpr double dLMD = h / total_length / m * dt;
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
double lambda_selection_constant(double min, double max)
{
    // random var.
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    double p = rand(mt);
    return p * min + (1 - p) * max;
}

std::vector<std::vector<double>> mat_layer(double k0, double theta, double V, double D)
{
    std::vector<std::vector<double>> M(2, std::vector<double>(2, 0.));

    double kx_0 = k0 * std::sin(theta);
    double kx_in_square = kx_0 * kx_0 - 2. * m * V / (hbar * hbar);

    if (kx_in_square < 0.0) // n in imaginary
    {
        double kx_in = sqrt(-kx_in_square);
        double nx = kx_in / kx_0;
        double delta = kx_0 * D;

        M.at(0).at(0) = std::cosh(nx * delta);
        M.at(0).at(1) = std::sinh(nx * delta) / nx;
        M.at(1).at(0) = std::sinh(nx * delta) * nx;
        M.at(1).at(1) = std::cosh(nx * delta);
    }
    else    // n in real
    {
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

/** 屈折率のRe,Imの判定　使い所がよくわかってないので残す **/
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
/* 真空から入社して Vg の物質へ抜けるときの反射率 */
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
        return ue / sita;
    }
    else
    {
        // ngが実数の時
        double kx_g = sqrt(kx_g_square);
        double n0 = 1.;
        double ng = kx_g / kx_0;

        std::complex<double> ue(-n0 * ng * B - C, ng * A - n0 * D);
        std::complex<double> sita(-n0 * ng * B + C, -ng * A - n0 * D);
        return ue / sita;
    }
}

/* transparent factor */
std::complex<double> transparent(double k0, double theta, double Vg, double A, double B, double C, double D, double z, std::complex<double> R)
{
    double kx_0 = k0 * std::sin(theta);
    std::complex<double> kx_g, eik0z, eikgz, T;
    std::complex<double> I(0, 1.);
    if (kx_0 * kx_0 - 2. * m * Vg / (hbar * hbar) < 0.)
    {
        std::cerr << "transparent kx_g is imaginary" << std::endl;
        return 1;
    }
    kx_g = std::sqrt(kx_0 * kx_0 - 2. * m * Vg / (hbar * hbar));

    eik0z = std::exp(I * kx_0 * z);
    eikgz = std::exp(I * kx_g * z);

    T = (A * (eik0z + R / eik0z) + I * kx_0 * B * (eik0z - R / eik0z)) / eikgz;

    return T;
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

int sim_copy()
{
    // PREPARATIONS
    std::cout << "///// PREPARATIONS /////" << std::endl;

    // variables
    int count[2] = {};
    int j, k;
    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double r;
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    double maxR = 0, minR = 1, maxT = 0, minT = 1;

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
    // Eigen::Vector2cd vec;

    // loop
    for (j = 0; j < particle_count; j++)
    {
        lambda = lambda_selection_constant(LMD_MIN, LMD_MAX);
        double k0 = 2 * pi / lambda;

        // mirror matrix
        double zeta;
        std::vector<std::vector<double>> M_ni = mat_layer(k0, theta, V_ni, D_ni); // Niの行列
        std::vector<std::vector<double>> M_ti = mat_layer(k0, theta, V_ti, D_ti); // Tiの行列
        std::vector<std::vector<double>> layermatrix_ni_ti(2, std::vector<double>(2, 0.)); // 二個分の行列
        std::vector<std::vector<double>> refmatrix_ni_ti(2, std::vector<double>(2, 0.));   // 多層膜の行列
        layermatrix_ni_ti = DOT(M_ti, M_ni);
        // まずは単位行列
        refmatrix_ni_ti.at(0).at(0) = 1.;
        refmatrix_ni_ti.at(0).at(1) = 0.;
        refmatrix_ni_ti.at(1).at(0) = 0.;
        refmatrix_ni_ti.at(1).at(1) = 1.;
        for (int l = 0; l < N_bilayer; l++)
        {
            refmatrix_ni_ti = DOT(layermatrix_ni_ti, refmatrix_ni_ti);
            zeta += D_ni + D_ti;
        }
        double E = refmatrix_ni_ti.at(0).at(0);
        double F = refmatrix_ni_ti.at(0).at(1);
        double G = refmatrix_ni_ti.at(1).at(0);
        double H = refmatrix_ni_ti.at(1).at(1);

        R = reflect(k0, theta, V_sio2, E, F, G, H);
        T = transparent(k0, theta, V_sio2, E, F, G, H, zeta + D_SiO2, R);

        // phase calculation
        // Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance * sin(angle_delta) / tan(2 * theta) * lambda;
        // Phi_a_main = 4 * pi * gap / lambda * theta0;
        // Phi_g_sub = 2 * pi * g * pow(m / h, 2) * theta0 / 2 * pow(gap / sin(theta), 2) * lambda;
        // Phi_a_sub = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * lambda * theta0;
        // Phase = Phi_g_main + Phi_g_sub;

        // probability calculation
        probO = std::norm(T);
        probH = std::norm(R);
        if (probO > maxT)
            maxT = probO;
        if (probO < minT)
            minT = probO;
        if (probH > maxR)
            maxR = probH;
        if (probH < minR)
            minR = probH;
        r = rand(mt);
        if (r <= probO)
        {
            file_O << lambda << std::endl;
            channel = O_TDC;
            tree->Fill();
            count[0]++;
        }
        if (r <= probH)
        {
            file_H << lambda << std::endl;
            channel = H_TDC;
            tree->Fill();
            count[1]++;
        }

        if (j % DISPLAY_INTERVAL == DISPLAY_INTERVAL - 1)
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
    std::cout << "R is between " << minR << " and " << maxR << std::endl;
    std::cout << "T is between " << minT << " and " << maxT << std::endl;

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
