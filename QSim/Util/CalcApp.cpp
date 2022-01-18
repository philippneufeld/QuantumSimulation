// Philipp Neufeld, 2021-2022

#include <iostream>
#include <cstdio>
#include <chrono>

#include "CalcApp.h"
#include "Argparse.h"

namespace QSim
{

    CalcApp::CalcApp()
    {

    }

    CalcApp::~CalcApp()
    {

    }

    int CalcApp::Run(int argc, const char** argv)
    {
        // parse command line arguents
        QSim::ArgumentParser parser;
        parser.AddOptionDefault("f,file", "I/O data filename.", "./data.dat");
        parser.AddOption("h,help", "Print this help string.");
        parser.AddOption("noplot", "Don't plot the results.");
        parser.AddOption("nocalc", "Just load the data from the specified file and skip the calculation");
        auto cmdArgs = parser.Parse(argc, argv);

        if (cmdArgs.IsError())
        {
            std::cout << cmdArgs.GetError() << std::endl;
            return -1;
        }
        else if (cmdArgs.IsOptionPresent("help"))
        {
            std::cout << parser.GetHelpString() << std::endl;
            return 0;
        }

        // load data file
        std::string filePath = cmdArgs.GetOptionStringValue("file");
        m_data.LoadFromFile(filePath);

        // run calculation
        if (!cmdArgs.IsOptionPresent("nocalc"))
        {
            std::cout << "Starting calculation..." << std::endl;
            auto start_ts = std::chrono::high_resolution_clock::now();

            this->DoCalculation();

            auto end_ts = std::chrono::high_resolution_clock::now();
            std::cout << "Finished calculation. (Elapsed time: " 
                << (std::chrono::high_resolution_clock::now() - start_ts).count() / 1.0e9 
                << "s)" << std::endl;
            
            m_data.SaveToFile(filePath);
        }

        // plotting
        if (!cmdArgs.IsOptionPresent("noplot"))
            this->Plot();

        return 0;
    }

    TDynamicMatrix<double> CalcApp::LoadMatrix(const std::string& name)
    {
        TDynamicMatrix<double> mat;
        if (m_data.GetDataTypeId(name) == 1)
        {
            const Memory& memory = m_data.GetData(name);
            auto ui64Ptr = static_cast<const std::uint64_t*>(memory.GetData());
            mat.Resize(ui64Ptr[0], ui64Ptr[1]);
            std::copy_n(reinterpret_cast<const double*>(ui64Ptr + 2), mat.Size(), &mat[0]);
        }
        return mat;
    }

    void CalcApp::StoreMatrixHelper(const std::string& name, std::uint64_t rows, 
        std::uint64_t cols, const double* data)
    {
        Memory memory;
        memory.Allocate(2*sizeof(std::uint64_t) + rows*cols*sizeof(double));
        auto ui64Ptr = static_cast<std::uint64_t*>(memory.GetData());
        ui64Ptr[0] = rows;
        ui64Ptr[1] = cols;
        std::copy_n(data, rows*cols, reinterpret_cast<double*>(ui64Ptr + 2));
        
        m_data.SetDataEx(name, memory.GetData(), memory.GetSize(), 1);
    }

    /*std::vector<char> CalcApp::SerializeData() const
    {
        // check size to allocate
        std::size_t size = 0;
        for (const auto& data: m_data)
        {
            const auto& name = data.first;
            const auto& dat = data.second;
            size += 3 * sizeof(std::uint64_t) + name.length() + sizeof(double) * dat.Rows() * dat.Cols();
        }

        std::vector<char> ser;
        ser.reserve(size);
        auto it = std::back_inserter(ser);

        for (const auto& data: m_data)
        {
            const auto& name = data.first;
            const auto& dat = data.second;

            // write data name
            std::uint64_t nameLen = name.length();
            std::copy_n(reinterpret_cast<const char*>(&nameLen), sizeof(nameLen), it);
            std::copy_n(name.data(), nameLen, it);

            // write rows, cols
            std::uint64_t rows = dat.Rows();
            std::uint64_t cols = dat.Cols();
            std::copy_n(reinterpret_cast<const char*>(&rows), sizeof(rows), it);
            std::copy_n(reinterpret_cast<const char*>(&cols), sizeof(cols), it);
            
            // write raw data
            std::copy_n(reinterpret_cast<const char*>(dat.Data()), sizeof(double)*rows*cols, it);
        }
        
        return ser;
    }

    std::map<std::string, TDynamicMatrix<double>> CalcApp::DeserializeData(
            const std::vector<char>& ser) const
    {
        std::map<std::string, TDynamicMatrix<double>> dataMap;

        for (auto it = ser.begin(); it < ser.end();)
        {
            // read data name
            std::uint64_t size = 0;
            std::copy_n(it, sizeof(size), reinterpret_cast<char*>(&size));
            std::string name(size, 0);
            it += sizeof(size);
            std::copy_n(it, size, &name[0]);
            it += size;

            // read rows, cols
            std::uint64_t rows = 0;  
            std::copy_n(it, sizeof(rows), reinterpret_cast<char*>(&rows));
            it += sizeof(rows);
            std::uint64_t cols = 0;
            std::copy_n(it, sizeof(cols), reinterpret_cast<char*>(&cols));
            it += sizeof(cols);
            
            // read raw data
            TDynamicMatrix<double> data(rows, cols, reinterpret_cast<const double*>(&(*it)));
            it += sizeof(double) * rows * cols;

            // insert to map
            dataMap[name] = std::move(data);
        }
        
        return dataMap;
    }

    bool CalcApp::SaveDataToFile(const std::string& path) const
    {
        std::FILE* pFile = std::fopen(path.c_str(), "wb");
        if (!pFile)
            return false;

        // write version
        std::uint64_t version[] = {MajorVersion_v, MinorVersion_v};
        std::fwrite(version, sizeof(std::uint64_t), 2, pFile);
        
        // write data block
        std::vector<char> data = SerializeData();
        std::uint64_t dataSize = data.size();
        std::fwrite(&dataSize, sizeof(dataSize), 1, pFile);
        std::fwrite(data.data(), 1, data.size(), pFile);

        std::fclose(pFile);
        return true;
    }
    
    bool CalcApp::LoadDataFromFile(const std::string& path)
    {
        std::FILE* pFile = std::fopen(path.c_str(), "rb");
        if (!pFile)
            return false;

        std::uint64_t size = 0;
        std::vector<char> ser;

        // read version
        std::uint64_t version[2];
        std::fread(version, sizeof(std::uint64_t), 2, pFile);

        // read data block
        std::fread(&size, sizeof(size), 1, pFile);
        ser.resize(size);
        std::fread(&ser[0], 1, size, pFile);
        m_data = DeserializeData(ser);

        return true;
    }*/

}
