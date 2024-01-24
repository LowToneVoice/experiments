#include"montecarlo.h"

int main()
{
    if (sim_lambda_roi() != 0)
    {
        cerr << "Error" << endl;
        return 1;
    }

    read_lambda();

    oscillation_lambda();

    return 0;
}
