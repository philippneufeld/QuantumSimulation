// Philipp Neufeld, 2021-2022

#include <iostream>
#include <QSim/Util/DataFile2.h>

int main(int argc, const char* argv[])
{
    QSim::DataFileDriver driver;

    std::cout << (driver.OpenFile("test.dat") ? "success" : "failure") << std::endl;

#ifdef QSIM_PYTHON3
    std::cout << "QSIM_PYTHON3" << std::endl;
#endif

    /*auto file = driver.GetOSHandle().lock();

    std::cout << file->GetCursorPos() << std::endl;

    std::cout << (file->SetCursorPos(2) ? "success" : "failure") << std::endl;
    std::cout << file->GetCursorPos() << std::endl;

    file->Write("aaa", 3);

    std::cout << file->GetCursorPos() << std::endl;
    
    std::cout << (file->SetCursorPos(20) ? "success" : "failure") << std::endl;
    std::cout << file->GetCursorPos() << std::endl;
    
    file->Write("bbb", 3);


    std::cout << file->GetCursorPos() << std::endl;*/

    return 0;
}
