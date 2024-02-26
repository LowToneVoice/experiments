#include <iostream>
#include <fstream>
#include <random>
#include <complex>
#include <Eigen/Dense>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>

using namespace std;
using namespace Eigen;

// FITTING RANGE AND INITIAL CONDITIONS
constexpr double fit_width_peak = 1.5e11;
constexpr double fit_peak_height = 30;

// CHANNELS
#define H_TDC 0
#define O_TDC 1
#define H_TDC_Ccut 2
#define H_TDC_Dcut 3
#define O_TDC_Ccut 4
#define O_TDC_Dcut 5

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
// beam
constexpr double lambda_min = 2.0e-10;
constexpr double lambda_max = 10e-10;
constexpr int daq_freq = 62.5e6;
constexpr double lambda_fade_width = .01 * (lambda_max - lambda_min);
// constexpr double lambda_fade_width = 1e-12;
constexpr double whole_beam_intensity_perSec_perArea = 5.4e8; // /s/m2
// ethalones
constexpr double delta_D = 1.;
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;
constexpr double gap = 189e-6;
constexpr int N_bilayer = 8;
constexpr double ethalone_height = 12e-3;
constexpr double slit_width = 0.2e-3;
// TDC
// constexpr double TDC_area = 5065049 / whole_beam_intensity_perSec_perArea;

// ALIGNMENT VARIABLES
// ethalone
constexpr double theta = 1.05 * pi / 180;
constexpr double mirror_distance = 150e-3;
// TDC
constexpr double total_length = 1.;
constexpr int daq_downsizing = 256;
constexpr double dt = 1. / daq_freq * daq_downsizing;

