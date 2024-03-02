#include <iostream>
#include <fstream>

int main()
{
    double start_theta = 179;
    double end_theta = 181;
    double step_theta = .05;
    int time_perStep_sec_theta = 60;


    std::ofstream outFile("theta_etalon_survey.m");

    std::string target_obj = "theta_etalon";
    std::string command = "setPosition";
    int numx = (int)((end_X - start_X) / step_X);
    int num_theta = (int)((end_theta - start_theta) / step_theta);

    int KPperTime = 25;
    int count = 0;

    if (outFile.is_open())
    {
        outFile << "// target_obj:yaw_etalon" << std::endl;
        outFile << "// command:" + command << std::endl;

        outFile << "system startLogging;" << std::endl;
        outFile << std::endl;

        for (int i = 0; i < numx; i++)
        {
            for (int j = 0; j < num_theta; j++)
            {
                /* code */
            }

            outFile << target_obj << " " << command << ":" << std::to_string(start + step * i) << ";" << std::endl;
            outFile << "niki startDetection;" << std::endl;
            outFile << "gate waitWhileKP:" << std::to_string(KPperTime * time_perStep_sec) << ";" << std::endl;
            count += KPperTime * time_perStep_sec;
            outFile << "niki stopDetection;" << std::endl;

            outFile << std::endl;
        }

        outFile << "system stopLogging;" << std::endl;
        outFile << "gate GatenetLogRename:;" << std::endl;

        std::cout
            << "Completed writing files.\n";
        std::cout << "whole time: " << count / 60 / 25 << std::endl;
    }
    else
    {
        std::cerr << "Cannot open the file.\n";
        return 1;
    }
}
