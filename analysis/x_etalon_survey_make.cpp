#include <iostream>
#include <fstream>

int main()
{
    double start = 51.0;
    double end = 55.0;
    double step = .1;
    int time_perStep = 20;

    std::ofstream outFile("x_etalon_survey.m");

    std::string target_obj = "x_etalon";
    std::string command = "setPosition";
    int numx = (int)((end - start) / step);

    int KPperTime = 25;
    int count = 0;

    if (outFile.is_open())
    {
        outFile << "// target_obj:yaw_etalon" << std::endl;
        outFile << "// command:" + command << std::endl;

        outFile << "system startLogging;" << std::endl;
        outFile << "niki startDetection;" << std::endl;
        outFile << std::endl;

        for (int i = 0; i < numx; i++)
        {
            outFile << target_obj << " " << command << ":" << std::to_string(start + step * i) << ";" << std::endl;
            outFile << "gate waitWhileKP:" << std::to_string(KPperTime * time_perStep) << ";" << std::endl;
            count += KPperTime * time_perStep;

            outFile << std::endl;
        }

        outFile << "niki stopDetection;" << std::endl;
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
