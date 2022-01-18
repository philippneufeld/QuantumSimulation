// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_CalcApp2_H_
#define QSim_Util_CalcApp2_H_

#include <string>
#include <map>
#include <vector>
#include <cstdint>
#include <type_traits>
#include <memory>

#include "../Math/Matrix.h"
#include "DataFile3.h"

namespace QSim
{
    
    class SimulationApp
    {
    public:
        SimulationApp();
        ~SimulationApp();

        virtual void Init(DataFileGroup& simdata) = 0;
        virtual void Continue(DataFileGroup& simdata)  = 0;
        virtual bool IsFinnished(DataFileGroup& simdata) = 0;
        virtual void Plot(DataFileGroup& simdata) = 0;

        int Run(int argc, const char** argv);

    private:
        std::string ExtractProgramName(int argc, const char** argv);
        std::string GenerateDefaultFilename(int argc, const char** argv);
    };

}

#endif
