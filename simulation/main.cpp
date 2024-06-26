#include "montecarlo.h"

int main()
{
    // DATA LABELS
    string Phase_contrib = "mix";
    string Main_sub = "mix";
    string Angle_delta_deg = "0";
    string Time_min = "10";
    string Angle_from_parallel_deg = "3e-3";
    string Lmd_used_min = "7e-10";
    string Lmd_used_max = "8e-10";
    string File_extension = "pdf";

    sim_lambda_roi(Phase_contrib, Main_sub, Time_min, Angle_delta_deg, Angle_from_parallel_deg, Lmd_used_min, Lmd_used_max);

    // read_lambda_perSec(Phase_contrib, Main_sub, Time_min, Angle_delta_deg, Angle_from_parallel_deg, Lmd_used_min, Lmd_used_max, File_extension);

    oscillation_lambda(Phase_contrib, Main_sub, Time_min, Angle_delta_deg, Angle_from_parallel_deg, Lmd_used_min, Lmd_used_max, File_extension);

    // theoretical_lambda(Phase_contrib, Main_sub, Time_min, Angle_delta_deg, Angle_from_parallel_deg, Lmd_used_min, Lmd_used_max, File_extension);

    return 0;
}
