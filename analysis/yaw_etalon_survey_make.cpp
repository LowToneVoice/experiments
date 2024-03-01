#include <iostream>
#include <fstream>

int main()
{
    double start = 6.465;
    double end = 6.495;
    double step = .001;
    int time_perStep = 5 * 60;
    double reference = 6.48;

    std::ofstream outFile("yaw_etalon_survey.m");

    std::string target_obj = "yaw_etalon";
    std::string command = "setPosition";
    int numx = (int)((end - start) / step);

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
            outFile << target_obj << " " << command << ":" << std::to_string(reference) << ";" << std::endl;
            outFile << "niki startDetection;" << std::endl;
            outFile << "gate waitWhileKP:" << std::to_string(KPperTime * time_perStep) << ";" << std::endl;
            count += KPperTime * time_perStep;
            outFile << "niki stopDetection;" << std::endl;

            outFile << std::endl;

            outFile << target_obj << " " << command << ":" << std::to_string(start + step * i) << ";" << std::endl;
            outFile << "niki startDetection;" << std::endl;
            outFile << "gate waitWhileKP:" << std::to_string(KPperTime * time_perStep) << ";" << std::endl;
            count += KPperTime * time_perStep;
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
