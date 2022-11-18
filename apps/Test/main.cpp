// Philipp Neufeld, 2021-2022

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include <QSim/Execution/ThreadPool.h>

using namespace QSim;
using namespace Eigen;

using namespace std::chrono_literals;

int main()
{
    ThreadPool pool(5);

    pool.Submit([](){ std::this_thread::sleep_for(2s); });
    for (int i=0; i<4; i++) pool.Submit([](){ std::this_thread::sleep_for(5s); });

    pool.WaitUntilReadyForTask();
    std::cout << "Ready for task" << std::endl;

    pool.WaitUntilFinished();
    std::cout << "Finished!" << std::endl;

    return 0;
}
