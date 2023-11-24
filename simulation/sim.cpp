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

// ALIGNMENTS CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define DAQ_FREQ 62.5e+6
#define THICKNESS 189e-6
#define LMD_MIN 2e-10
#define LMD_MAX 9e-10
#define BEAM_INTENSITY 9.4e+7   // /s/cm2
#define BEAM_SIZE (5.5 * 4.5)   // cm2

// PHYS & MATH CONSTANTS: !!!!! DO NOT CHANGE !!!!!
#define pi M_PI
#define m 1.6749e-27
#define g 9.8
#define h 6.62607015e-34
#define rho 9.1e+28
#define bc 3.4e-15

int sim()
{
    // PREPARATIONS
    // Alignment variables
    double daq_downsizing = 1;  // >= 1
    double theta = 1.05 * pi / 180, theta0 = 1e-3 * pi / 180;
    double mirror_distance = 150e-3, total_length = 1;
    double beam_time = 30 * 60; // sec
    double alpha = .5, gamma = .5;

    // Alignment calculations
    double beam_count = BEAM_INTENSITY * BEAM_SIZE * beam_time;
    double dt = 1 / DAQ_FREQ * daq_downsizing;
    double dLMD = h / total_length / m * dt;
    int L_max = (int)((LMD_MAX - LMD_MIN) / dLMD);

    // Variable resets
    beam_count = 1e+4;
    // L_max = 10000;
    std::cout << "///// ALIGNMENT CALCULATIONS /////" << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "dLMD = " << dLMD << std::endl;
    std::cout << "L_max = " << L_max << std::endl;
    std::cout << "beam_count = " << beam_count << std::endl;
    std::cout << "max datasize = " << MAX_DATA_SIZE << std::endl;

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
    // variables
    int count[2] = {};
    int i;
    double lmd;
    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double r;
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    int H_TDC = 0, O_TDC = 1, H_ADC = 2, O_ADC = 3;
    // root file
    Double_t lambda;
    Int_t channel;
    TFile *file = new TFile(OUTPUT_TREE, "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("channel", &channel);
    tree->Branch("lambda", &lambda);

    // loop
    for (i = 0; i < (int)beam_count; i++)
    {
        lmd = LMD_MIN;
        while (lmd <= LMD_MAX)
        {
            // phase calculation
            Phi_g_main = -2 * pi * g * pow(m / h, 2) * 2 * THICKNESS * mirror_distance / tan(2 * theta) * lmd;
            Phi_a_main = 4 * pi * THICKNESS / lmd * theta0;
            Phi_g_sub = 2 * pi * g * pow(m / h, 2) * theta0 / 2 * pow(THICKNESS / sin(theta), 2) * lmd;
            Phi_a_sub = -4 * pi * THICKNESS * rho * bc / 2 / pi / pow(theta, 2) * lmd * theta0;
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
