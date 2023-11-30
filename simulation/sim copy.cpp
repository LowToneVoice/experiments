#include<iostream>
#include<fstream>
#include<math.h>
#include<random>
#include<TFile.h>
#include<TTree.h>

#define MAX_DATA_SIZE 1e6
#define OUTPUT_TREE "tree.root"
#define OUTPUT_O "outO.dat"
#define OUTPUT_H "outH.dat"

// PHYS & MATH CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define pi M_PI
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
constexpr double V_sio2 = 90.9e-09 * 1.6022e-19;
constexpr double V_ni = 224.e-09 * 1.6022e-19;
constexpr double V_ti = -40.e-9 * 1.6022e-19;

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
constexpr double D_ni = 133.5e-10 / delta_D;
constexpr double D_ti = 98.3e-10 / delta_D;
constexpr int N_bilayer = 8;

// ALIGNMENT VARIABLES
constexpr int daq_downsizing = 1;
constexpr double theta = 1.05 * pi / 180;
constexpr double theta0 = 1.e-3 * pi / 180;
constexpr double mirror_distance = 150e-3;
constexpr double total_length = 1;
constexpr double beam_time = 30 * 60; // sec
// constexpr int beam_count = beam_time * BEAM_SIZE;
constexpr int beam_count = 10000;
constexpr double dt = 1 / daq_freq * daq_downsizing;
constexpr double dLMD = h / total_length / m * dt;
// constexpr int L_max = (int)((LMD_MAX - LMD_MIN) / dLMD);
constexpr int L_max = 10000;


// TODO: #1 building alpha_calc in sim copy.cpp
double alpha_calc(double lambda)
{
}

int sim()
{
    // PREPARATIONS
    std::cout << "///// PREPARATIONS /////" << std::endl;

    // variables
    double alpha = .5, gamma = .5;
    int count[2] = {};
    int i;
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
        return 0;
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

    // loop
    for (i = 0; i < beam_count; i++)
    {
        lmd = LMD_MIN;
        while (lmd <= LMD_MAX)
        {
            // phase calculation
            Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * lmd;
            Phi_a_main = 4 * pi * gap / lmd * theta0;
            Phi_g_sub = 2 * pi * g * pow(m / h, 2) * theta0 / 2 * pow(gap / sin(theta), 2) * lmd;
            Phi_a_sub = -4 * pi * gap * rho * bc / 2 / pi / pow(theta, 2) * lmd * theta0;
            Phase = Phi_g_main + Phi_g_sub + Phi_a_main + Phi_a_sub;

            // probability calculation
            r = rand(mt);
            if (r < alpha * (1 + cos(Phase)))
            {
                file_O << lmd << std::endl;
                channel = O_TDC;
                lambda = lmd;
                tree->Fill();
                count[0]++;
            }
            else if (r > 1 - (gamma - alpha * cos(Phase)))
            {
                file_H << lmd << std::endl;
                channel = H_TDC;
                lambda = lmd;
                tree->Fill();
                count[1]++;
            }

            lmd += dLMD;
        }
        // loop count
        if (i % 100 == 99)
        {
            std::cout << "loop " << i + 1 << " / " << (int)beam_count << " has ended." << std::endl;
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

    return 1;
}
