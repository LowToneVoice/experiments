#include"analysis.h"

int analysis()
{
    string input_tree = "../simulation/BL05/20231224225335_list_copy.root";
    string output_tree = "./data/20231224225335.root";
    string img_tof = "./img/count/tof/20231224225335.pdf";
    string img_lmd = "./img/count/lmd/20231224225335.pdf";
    lambda_roi_hist(1, .41, .44, .34, .48, .52, .54, .32, .48, input_tree, output_tree, img_tof, img_lmd);

    return 0;
}
