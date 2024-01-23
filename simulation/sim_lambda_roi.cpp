#include <iostream>
#include <fstream>
#include <random>
#include <Eigen/Dense>
#include <complex>
#include <TFile.h>
#include <TTree.h>

using namespace std;
using namespace Eigen;

#define OUT_INTERVAL 1000

// data labels
string PHASE_CONTRIB = "mix";
string MAIN_SUB = "mix";
string ANGLE_DELTA_DEG = "30";
string TIME_MIN = "60";
string ANGLE_FROM_PARALLEL_DEG = "1e-1";
string LMD_USED_MIN = "7e-10";
string LMD_USED_MAX = "10e-10";

const bool REFLECTION = true;

// CHANNELS
#define H_TDC 0
#define O_TDC 1
#define H_ADC 2
#define O_ADC 3

// PHYS & MATH CONSTANTS
complex<double> I(0, 1.);
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

constexpr double ethalone_height = 12e-3;
constexpr double slit_width = 260e-6;

constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 10e-10;
constexpr double theta_min = .2 * pi / 180;
constexpr double theta_max = 1.5 * pi / 180;

constexpr int daq_freq = 62.5e6;

constexpr double lambda_fade_width = .01 * (lambda_max - lambda_min);
// constexpr double lambda_fade_width = 1e-12;

constexpr double whole_beam_intensity_perSec_perArea = 5.4e8; // /s/m2
constexpr double TDC_area = 5065049 / whole_beam_intensity_perSec_perArea;

// ALIGNMENT VARIABLES
constexpr double theta = 1.05 * pi / 180;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1.;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
// constexpr double d_lambda = 0.00001e-10;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;

constexpr int N_loop_lambda_fade = (int)(lambda_fade_width / d_lambda);

// beam count per sec
double beam_count_per_sec(double lambda)
{
    // in case lambda is over 0.9 nm
    if (lambda > 9e-10)
        return 0;

    const double L_tof_sample = 18022e-03;
    const double tof_max_sample = 45000e-6;
    const double sample_division = 1000;
    const double d_lambda_sample = (tof_max_sample / sample_division) * h / L_tof_sample / m;

    double position[][2] = {{.216e-9, 0}, {.265e-9, 2.57}, {.4e-9, 1.5}, {.7e-9, .38}, {.8e-9, .22}, {.94e-9, 0}};

    // count_per_sec := k1 * lambda + k2 at d_lambda = d_lambda = 9.88e-13
    double k1[5], k2[5];
    double count;
    for (int i = 0; i < 5; i++)
    {
        k1[i] = (position[i + 1][1] - position[i][1]) / (position[i + 1][0] - position[i][0]);
        k2[i] = position[i][1] - k1[i] * position[i][0];

        if (position[i][0] <= lambda && lambda < position[i + 1][0])
        {
            count = k1[i] * lambda * (d_lambda / d_lambda_sample) + k2[i] * (d_lambda / d_lambda_sample);
            if (count > 1)
            {
                cerr << "count is over 1 / s. lambda is " << lambda << endl;
            }
            return count * (ethalone_height * slit_width / TDC_area);
        }
    }

    cerr << "lambda is out of range" << endl;
    return 0;
}

// matrix M of one layer
Matrix2d mat_layer(double k0, double theta, double V, double D)
{
    Matrix2d M;
    double kx0 = k0 * sin(theta);
    complex<double> kx = sqrt(kx0 * kx0 - 2. * m * V / (hbar * hbar));
    complex<double> nx = kx / kx0;
    double delta = kx0 * D;

    M(0, 0) = cos(nx * delta).real();
    M(0, 1) = (sin(nx * delta) / nx).real();
    M(1, 0) = (-sin(nx * delta) * nx).real();
    M(1, 1) = cos(nx * delta).real();

    return M;
}

// complex reflection factor
complex<double> reflect(double k0, double theta, double Vg, Matrix2d M)
{
    double kx0 = k0 * sin(theta);
    complex<double> kxg = sqrt(kx0 * kx0 - 2. * m * Vg / (hbar * hbar));
    double n0 = 1.;
    complex<double> ng = kxg / kx0;
    complex<double> numer, denom;

    numer = -n0 * ng * M(0, 1) - M(1, 0) + I * (ng * M(0, 0) - n0 * M(1, 1));
    denom = -n0 * ng * M(0, 1) + M(1, 0) + I * (-ng * M(0, 0) - n0 * M(1, 1));

    return numer / denom;
}

// complex transparent factor
complex<double> transparent(double k0, double theta, double Vg, double zeta, Matrix2d M)
{
    double kx0 = k0 * sin(theta);
    complex<double> kxg = sqrt(kx0 * kx0 - 2. * m * Vg / (hbar * hbar));
    double n0 = 1.;
    complex<double> ng = kxg / kx0;

    complex<double> eikgz = cos(kxg * zeta) + I * sin(kxg * zeta);
    complex<double> numer, denom;
    numer = 2. * I * n0 / eikgz;
    denom = n0 * ng * M(0, 1) - M(1, 0) + I * (n0 * M(1, 1) + ng * M(0, 0));

    return numer / denom;
}

int sim_lambda_roi(
    string phase_contrib_input = PHASE_CONTRIB,
    string main_sub_input = MAIN_SUB,
    string time_min_input = TIME_MIN,
    string angle_delta_deg_input = ANGLE_DELTA_DEG,
    string angle_from_parallel_deg_input = ANGLE_FROM_PARALLEL_DEG,
    string lmd_used_min_input = LMD_USED_MIN,
    string lmd_used_max_input = LMD_USED_MAX
)
{
    // file names
    string FORMAT = phase_contrib_input + "_" + main_sub_input + "_" + angle_delta_deg_input + "deg_" + time_min_input + "min_ALPHA" + angle_from_parallel_deg_input + "_lmd" + lmd_used_min_input + "to" + lmd_used_max_input;
    string OUTPUT_TREE = "./dat/montecarlo/root/" + FORMAT + ".root";
    string OUTPUT_H = "./dat/montecarlo/dat/" + FORMAT + "_H.dat";
    string OUTPUT_O = "./dat/montecarlo/dat/" + FORMAT + "_O.dat";
    string OUTPUT_PROB = "./dat/theoretical/" + FORMAT + ".dat";

    // main / sub selection
    int G_MAIN = 0, G_SUB = 0, A_MAIN = 0, A_SUB = 0;
    bool flag = false;
    if (PHASE_CONTRIB == "mix" || PHASE_CONTRIB == "g")
    {
        if (MAIN_SUB == "mix" || MAIN_SUB == "main")
        {
            G_MAIN = 1;
            flag = true;
        }
        if (MAIN_SUB == "mix" || MAIN_SUB == "sub")
        {
            G_SUB = 1;
            flag = true;
        }
    }
    if (PHASE_CONTRIB == "mix" || PHASE_CONTRIB == "a")
    {
        if (MAIN_SUB == "mix" || MAIN_SUB == "main")
        {
            A_MAIN = 1;
            flag = true;
        }
        if (MAIN_SUB == "mix" || MAIN_SUB == "sub")
        {
            A_SUB = 1;
            flag = true;
        }
    }
    if (flag == false)
    {
        cerr << "Simulation type is not correct." << endl;
        return 1;
    }

    // input value
    const double lambda_min_used = stod(lmd_used_min_input);
    const double lambda_max_used = stod(lmd_used_max_input);
    const double delta = stod(angle_delta_deg_input) * pi / 180;
    const double beam_time_sec = stod(time_min_input) * 60;
    const double angle_from_parallel = stod(angle_from_parallel_deg_input) * pi / 180;

    const int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);

    double lambda;
    double k0;
    double zeta;
    complex<double> R, T;
    complex<double> R_not_normed, T_not_nnormed;
    complex<double> wave_fncO, wave_fncH;
    double R2norm, T2norm, R_phase, T_phase;
    Matrix2d mat_bilayer_ni_ti;
    Matrix2d M;
    ofstream fileO(OUTPUT_O);
    ofstream fileH(OUTPUT_H);
    ofstream fileP(OUTPUT_PROB);
    if (!fileO.is_open() || !fileH.is_open() || !fileP.is_open())
    {
        cerr << "FILE DID NOT OPENED" << endl;
        return 1;
    }

    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double probability, probO, probH;
    mt19937 mt{(random_device{}())};
    uniform_real_distribution<double> rand(0., 1.);
    Double_t lmd;
    Int_t channel;
    TFile file(OUTPUT_TREE.c_str(), "recreate");
    TTree tree("tree", "tree");
    tree.Branch("channel", &channel);
    tree.Branch("lambda", &lmd);
    int count[2] = {};
    int beam_count = 0;
    int l;
    double max_prob = 0;

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
            M = Matrix2d::Identity();
            for (int l = 0; l < N_bilayer; l++)
            {
                M = mat_bilayer_ni_ti * M;
                zeta += D_ni + D_ti;
            }
            R_not_normed = reflect(k0, theta, V_sio2, M);
            T_not_nnormed = transparent(k0, theta, V_sio2, zeta, M);
            R = R_not_normed / sqrt(norm(R_not_normed) + norm(T_not_nnormed));
            T = T_not_nnormed / sqrt(norm(R_not_normed) + norm(T_not_nnormed));
            R2norm = norm(R);
            R_phase = arg(R);
            T2norm = norm(T);
            T_phase = arg(T);

            // Phase calculation
            Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * lambda * sin(delta);
            Phi_a_main = 4 * pi * gap / lambda * angle_from_parallel;
            Phi_g_sub = 2 * pi * g * pow(m / h, 2) * angle_from_parallel / 2 * pow(gap / sin(theta), 2) * lambda * sin(delta);
            Phi_a_sub = -4 * pi * gap * rho_Ni * bc_Ni / 2 / pi / pow(theta, 2) * lambda * angle_from_parallel; // maximize the value
            Phase = Phi_g_main * G_MAIN + Phi_g_sub * G_SUB + Phi_a_main * A_MAIN + Phi_a_sub * A_SUB;

            // probability calculation
            wave_fncH += (T * R * T * T * T + R * T * R * R * T * exp(I * Phase) + R * T * T * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncO += (T * R * T * R + R * T * R * T * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
        } while (l < N_loop_lambda_fade);

        probO = norm(wave_fncO);
        probH = norm(wave_fncH);

        fileP << lmd << " " << probH << " " << probO << endl;

        for (int time_sec = 0; time_sec < beam_time_sec; time_sec++)
        {
            double beam_entry_probability = beam_count_per_sec(lmd);
            if (max_prob < probO * beam_entry_probability)
                max_prob = probO * beam_entry_probability;
            if (max_prob < probH * beam_entry_probability)
                max_prob = probH * beam_entry_probability;

            // when beam does not enter
            if (rand(mt) > beam_entry_probability)
                continue;
            else
            {
                beam_count++;

                probability = rand(mt);
                if (probability < probO)
                {
                    fileO << lambda << endl;
                    channel = O_TDC;
                    tree.Fill();
                    count[0]++;
                }
                else if (1 - probability <= probH)
                {
                    fileH << lambda << endl;
                    channel = H_TDC;
                    tree.Fill();
                    count[1]++;
                }
            }
        }

        // progress display
        if (j % OUT_INTERVAL == OUT_INTERVAL - 1)
        {
            std::cout << j + 1 << " / " << N_loop_lambda << " has finished." << endl;
        }
    }

    tree.Write();
    fileO.close();
    fileH.close();
    fileP.close();
    file.Close();

    std::cout << "counts: " << count[0] << " and " << count[1] << endl;
    std::cout << "Phi_g_main = " << -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * sin(delta) << " * lambda" << endl;
    std::cout << "used N = " << setprecision(3) << float(beam_count * N_loop_lambda) << endl;
    std::cout << "max probability: " << max_prob << endl;

    return 0;
}

int main()
{
    return sim_lambda_roi();
}
