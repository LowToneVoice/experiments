#include <iostream>
#include <fstream>
#include <random>
#include <Eigen/Dense>
#include <complex>
#include <TFile.h>
#include <TTree.h>

#define OUT_INTERVAL 1000

// data labels
std::string PHASE_CONTRIB = "mix";
std::string MAIN_SUB = "mix";
std::string ANGLE_DELTA_DEG = "30";
std::string TIME_MIN = "10";
std::string ANGLE_FROM_PARALLEL_DEG = "6e-3";

// datafiles
std::string FORMAT = PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG;
std::string OUTPUT_TREE = "./dat/montecarlo/root/" + FORMAT + ".root";
std::string OUTPUT_H = "./dat/montecarlo/dat/" + FORMAT + "_H.dat";
std::string OUTPUT_O = "./dat/montecarlo/dat/" + FORMAT + "_O.dat";
std::string OUTPUT_PROB = "./dat/theoretical/" + FORMAT + ".dat";

const bool REFLECTION = true;

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
constexpr double lambda_min_used = 6.9e-10;
constexpr double lambda_max_used = 8.4e-10;
constexpr double theta = 1.05 * pi / 180;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1.;
constexpr int daq_downsizing = 16;
constexpr double dt = 1. / daq_freq * daq_downsizing;
// constexpr double d_lambda = 0.00001e-10;
constexpr double d_lambda = h / total_length / m * dt;
constexpr double d_theta = .05 * pi / 180;

constexpr int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);
constexpr int N_loop_theta = (int)((theta_max - theta_min) / d_theta);
constexpr int N_loop_lambda_fade = (int)(lambda_fade_width / d_lambda);

int position();
int lambda_roi_histogram();
int read_lambda();
int oscillation_lambda();

int main()
{
}
