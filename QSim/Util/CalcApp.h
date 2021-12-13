// Philipp Neufeld, 2021

#ifndef QSim_Util_CalcApp_H_
#define QSim_Util_CalcApp_H_

#include <string>
#include <map>
#include <vector>
#include <cstdint>

#include "../Math/Matrix.h"

namespace QSim
{
    
    struct CalcAppFigureDesc
    {
        std::string caption = "";
        std::string xLabel = "";
        std::string yLabel = "";
    };

    struct CalcAppPlotDesc
    {
        std::string figName = "";
        std::string xDataName = "";
        std::string yDataName = "";
        std::size_t xDataRow = 0;
        std::size_t yDataRow = 0;
    };


    class CalcApp
    {
        constexpr static std::uint64_t MajorVersion_v = 1;
        constexpr static std::uint64_t MinorVersion_v = 0;
    public:
        CalcApp();
        ~CalcApp();

        virtual void DoCalculation() { };
        int Run(int argc, const char** argv);

        // Data control
        template<typename MT>
        void SetData(const std::string& name, const TMatrix<MT>& data) { m_data[name] = data; }
        void RemoveData(const std::string& name) { m_data.erase(name); }

        // Plotting control
        void SetFigureDesc(const std::string& figName, 
            const CalcAppFigureDesc& desc) { m_figureDescs[figName] = desc; }
        void AddPlot(const CalcAppPlotDesc& desc) { m_plots.push_back(desc); }

    private:
        bool SaveDataToFile(const std::string& path) const;
        bool LoadDataFromFile(const std::string& path);
        
        std::vector<char> SerializeData() const;
        std::vector<char> SerializePlots() const;

        std::map<std::string, TDynamicMatrix<double>> DeserializeData(
            const std::vector<char>& ser) const;
        
    private:
        std::vector<CalcAppPlotDesc> m_plots;
        std::map<std::string, CalcAppFigureDesc> m_figureDescs;
        std::map<std::string, TDynamicMatrix<double>> m_data;
    };

}

#endif
