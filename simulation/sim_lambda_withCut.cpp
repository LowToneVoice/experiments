#include "simulate.h"

int sim_lambda_withCut(
    string phase_contrib_input = PHASE_CONTRIB,
    string main_sub_input = MAIN_SUB,
    string time_min_input = TIME_MIN,
    string angle_delta_deg_input = ANGLE_DELTA_DEG,
    string angle_from_parallel_deg_input = ANGLE_FROM_PARALLEL_DEG,
    string lmd_used_min_input = LMD_USED_MIN,
    string lmd_used_max_input = LMD_USED_MAX)
{
    // file names
    string file_format = FILE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);
    string OUTPUT_TREE = "./dat/montecarlo/root_cut/" + file_format + ".root";
    string OUTPUT_PROB = "./dat/theoretical_cut/" + file_format + ".dat";

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

    // const int N_loop_lambda = (int)((lambda_max_used - lambda_min_used) / d_lambda);
    const int N_loop_lambda = (int)((lambda_max - lambda_min) / d_lambda);

    double lambda;
    double k0;
    double zeta;
    complex<double> R, T;
    complex<double> R_not_normed, T_not_nnormed;
    complex<double> wave_fncO[3] = {}, wave_fncH[3] = {};
    double R2norm, T2norm, R_phase, T_phase;
    Matrix2d mat_bilayer_ni_ti;
    Matrix2d M;
    ofstream fileP(OUTPUT_PROB);
    if (!fileP.is_open())
    {
        cerr << "FILE DID NOT OPENED" << endl;
        return 1;
    }

    double Phi_g_main, Phi_a_main, Phi_g_sub, Phi_a_sub;
    double Phase;
    double probability, probO[3] = {}, probH[3] = {};
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
    double max_prob = 0, min_prob = 1;

    for (int j = 0; j < N_loop_lambda; j++)
    {
        lmd = lambda_min + d_lambda * j;
        // lmd = lambda_min_used + d_lambda * j;
        wave_fncH[0] = 0;
        wave_fncO[0] = 0;
        wave_fncH[1] = 0;
        wave_fncO[1] = 0;
        wave_fncH[2] = 0;
        wave_fncO[2] = 0;
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
            wave_fncH[0] += (T * R * T * T * T + R * T * R * R * T * exp(I * Phase) + R * T * T * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncH[1] += (T * R * T * T * T ) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncH[2] += (R * T * R * R * T * exp(I * Phase) + R * T * T * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncO[0] += (T * R * T * R + R * T * R * T * exp(I * Phase) + R * R * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncO[1] += (T * R * T * R) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
            wave_fncO[2] += (R * T * R * T * exp(I * Phase) + R * R * exp(I * Phase)) / (double)(N_loop_lambda_fade == 0 ? 1 : N_loop_lambda_fade);
        } while (l < N_loop_lambda_fade);

        probO[0] = norm(wave_fncO[0]);
        probH[0] = norm(wave_fncH[0]);
        probO[1] = norm(wave_fncO[1]);
        probH[1] = norm(wave_fncH[1]);
        probO[2] = norm(wave_fncO[2]);
        probH[2] = norm(wave_fncH[2]);

        fileP << lmd << " " << probH[0] << " " << probO[0] << " " << probH[1] << " " << probO[1] << " " << probH[2] << " " << probO[2] << endl;

        for (int particle = 0; particle < beam_time_sec * beam_count_per_sec(lmd); particle++)
        {
            beam_count++;
            probability = rand(mt);
            if (probability < min_prob)
                min_prob = probability;
            if (probability > min_prob)
                max_prob = probability;

            if (probability < probO[0])
            {
                channel = O_TDC;
                tree.Fill();
                count[0]++;
            }
            if (probability < probO[1])
            {
                channel = O_TDC_Ccut;
                tree.Fill();
            }
            if (probability < probO[2])
            {
                channel = O_TDC_Dcut;
                tree.Fill();
            }
            if (1 - probability <= probH[0])
            {
                channel = H_TDC;
                tree.Fill();
                count[1]++;
            }
            if (1 - probability <= probH[1])
            {
                channel = H_TDC_Ccut;
                tree.Fill();
            }
            if (1 - probability <= probH[2])
            {
                channel = H_TDC_Dcut;
                tree.Fill();
            }
        }

        // progress display
        if (j % OUT_INTERVAL == OUT_INTERVAL - 1)
        {
            std::cout << j + 1 << " / " << N_loop_lambda << " has finished." << endl;
        }
    }

    tree.Write();
    fileP.close();
    file.Close();

    std::cout << "counts: " << count[0] << " and " << count[1] << endl;
    std::cout << "Phi_g_main = " << -2 * pi * g * pow(m / h, 2) * 2 * gap * mirror_distance / tan(2 * theta) * sin(delta) << " * lambda" << endl;
    std::cout << "used N = " << setprecision(3) << beam_count << endl;
    std::cout << "min probability: " << min_prob << endl;
    std::cout << "max probability: " << max_prob << endl;

    return 0;
}

int main()
{
    sim_lambda_withCut();
}
