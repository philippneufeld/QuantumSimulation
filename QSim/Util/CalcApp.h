// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_CalcApp_H_
#define QSim_Util_CalcApp_H_

#include <string>
#include <map>
#include <vector>
#include <cstdint>
#include <type_traits>

#include "../Math/Matrix.h"
#include "DataFile.h"

namespace QSim
{
    
    class CalcApp
    {
        constexpr static std::uint64_t MajorVersion_v = 1;
        constexpr static std::uint64_t MinorVersion_v = 0;
    public:
        CalcApp();
        ~CalcApp();

        virtual void DoCalculation() { };
        virtual void Plot() { }
        int Run(int argc, const char** argv);

        // Data control
        void ClearData() { m_data.Clear(); }
        void RemoveData(const std::string& name) { m_data.RemoveData(name); }
        bool ContainsData(const std::string& name) const { return m_data.Contains(name); }
        std::vector<std::string> GetDataIndex() const { return m_data.GetIndex(); }

        template<typename MT, typename=void/*std::enable_if_t<std::is_same<Internal::TMatrixElementType_t<MT>, double>::value>*/>
        void StoreMatrix(const std::string& name, const TMatrix<MT>& data);
        TDynamicMatrix<double> LoadMatrix(const std::string& name);

    private:
        void StoreMatrixHelper(const std::string& name, std::uint64_t rows, 
            std::uint64_t cols, const double* data);

    private:
        DataFile m_data;
    };

    template<typename MT, typename>
    void CalcApp::StoreMatrix(const std::string& name, const TMatrix<MT>& data)
    {
        StoreMatrixHelper(name, (~data).Rows(), (~data).Cols(), (~data).Data());
    }

}

#endif
