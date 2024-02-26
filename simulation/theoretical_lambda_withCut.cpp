#include "simulate.h"

int theoretical_lambda_withCut(
    string phase_contrib_input = PHASE_CONTRIB,
    string main_sub_input = MAIN_SUB,
    string time_min_input = TIME_MIN,
    string angle_delta_deg_input = ANGLE_DELTA_DEG,
    string angle_from_parallel_deg_input = ANGLE_FROM_PARALLEL_DEG,
    string lmd_used_min_input = LMD_USED_MIN,
    string lmd_used_max_input = LMD_USED_MAX,
    string file_extension = FILE_EXTENSION)
{
    // datafiles
    string file_format = FILE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);
    FILE *gpl_file;
    string input_file = "./dat/theoretical_cut/" + file_format + ".dat";
    string output_beam_normal = "./beam_count/theoretical/cut/" + file_format + "." + file_extension;
    string output_beam_zoom = "./beam_count/theoretical/cut_zoom/" + file_format + "." + file_extension;
    string output_oscil_normal = "./oscil_graph/theoretical/cut/" + file_format + "." + file_extension;
    string output_oscil_zoom = "./oscil_graph/theoretical/cut_zoom/" + file_format + "." + file_extension;
    gpl_file = popen("gnuplot", "w");

    // input data
    double lambda_min_used = stod(lmd_used_min_input);
    double lambda_max_used = stod(lmd_used_max_input);
    string title_format = TITLE_FORMAT(phase_contrib_input, main_sub_input, time_min_input, angle_delta_deg_input, angle_from_parallel_deg_input, lmd_used_min_input, lmd_used_max_input);

    fprintf(gpl_file, "input='%s'\n", input_file.c_str());
    fprintf(gpl_file, "set term %s\n", file_extension.c_str());

    // normal beam
    fprintf(gpl_file, "set title '%s'\n", title_format.c_str());
    fprintf(gpl_file, "set output '%s'\n", output_beam_normal.c_str());
    fprintf(gpl_file, "set xrange [%e:%e]\n", lambda_min, lambda_max * 1.2);
    fprintf(gpl_file, "set xlabel 'lambda'\n");
    fprintf(gpl_file, "set ylabel 'probability'\n");
    fprintf(gpl_file, "plot input u 1:2 w l title 'H beam', ");
    fprintf(gpl_file, "input u 1:3 w l title 'O beam', ");
    fprintf(gpl_file, "input u 1:($4+$6) w l title 'O beam with cut', ");
    fprintf(gpl_file, "input u 1:($5+$7) w l title 'H beam with cut'\n");

    // zoom beam
    fprintf(gpl_file, "set output '%s'\n", output_beam_zoom.c_str());
    fprintf(gpl_file, "set xrange [%e:%e]\n", lambda_min_used, lambda_max_used);
    fprintf(gpl_file, "plot input u 1:2 w l title 'H beam', ");
    fprintf(gpl_file, "input u 1:3 w l title 'O beam', ");
    fprintf(gpl_file, "input u 1:($4+$6) w l title 'O beam with cut', ");
    fprintf(gpl_file, "input u 1:($5+$7) w l title 'H beam with cut'\n");

    // oscillation normal
    fprintf(gpl_file, "set output '%s'\n", output_oscil_normal.c_str());
    fprintf(gpl_file, "set xrange [%e:%e]\n", lambda_min, lambda_max);
    fprintf(gpl_file, "set ylabel '(I_H-I_O) / (I_H+I_O)'\n");
    fprintf(gpl_file, "set nokey\n");
    fprintf(gpl_file, "set yrange [-2:2]\n");
    fprintf(gpl_file, "plot input u 1:(($2/($4+$6)-$3/($5+$7))/($2/($4+$6)+$3/($5+$7))) w l\n");

    // oscillation zoom
    fprintf(gpl_file, "set output '%s'\n", output_oscil_zoom.c_str());
    fprintf(gpl_file, "set xrange [%e:%e]\n", lambda_min_used, lambda_max_used);
    fprintf(gpl_file, "plot input u 1:(($2/($4+$6)-$3/($5+$7))/($2/($4+$6)+$3/($5+$7))) w l\n");

    pclose(gpl_file);

    return 0;
}

int main()
{
    theoretical_lambda_withCut();
}
