// Philipp Neufeld, 2021

#include <iostream>
#include <cstdio>
#include <iterator>
#include <algorithm>

#include "DataFile.h"

namespace QSim
{

    DataFile::DataFile() { }
    DataFile::~DataFile() { }

    void DataFile::Clear()
    {
        m_data.clear();
    }

    void DataFile::SetData(const std::string& name, const void* data, std::size_t n)
    {
        SetDataEx(name, data, n, 0);
    }

    void DataFile::SetDataEx(const std::string& name, const void* data, std::size_t n, std::uint64_t typeId)
    {
        m_data[name] = std::make_pair(typeId, std::move(Memory(data, n)));
    }
    
    const Memory& DataFile::GetData(const std::string& name) const
    {
        static Memory nullMem; // no memory allocated (static because it is returned by reference)
        auto it = m_data.find(name);
        return it != m_data.end() ? it->second.second : nullMem;
    }

    std::uint64_t DataFile::GetDataTypeId(const std::string& name) const
    {
        auto it = m_data.find(name);
        return it != m_data.end() ? it->second.first : -1;
    }

    void DataFile::RemoveData(const std::string& name)
    {
        auto it = m_data.find(name);
        if (it != m_data.end())
            m_data.erase(it);
    }

    bool DataFile::Contains(const std::string& name) const
    {
        return (m_data.find(name) != m_data.end());
    }

    std::vector<std::string> DataFile::GetIndex() const
    {
        std::vector<std::string> keys;
        keys.reserve(m_data.size());
        auto it = std::back_inserter(keys);
        for (auto& el: m_data)
            *it++ = el.first;
        return keys;
    }

    bool DataFile::SaveToFile(const std::string& path) const
    {
        std::FILE* pFile = std::fopen(path.c_str(), "wb");
        if (!pFile)
            return false;

        // write version
        std::uint64_t version[] = {MajorVersion_v, MinorVersion_v};
        std::fwrite(version, sizeof(std::uint64_t), 2, pFile);
        
        // write number of data blocks
        std::uint64_t blockCnt = m_data.size();
        std::fwrite(reinterpret_cast<const void*>(&blockCnt), sizeof(blockCnt), 1, pFile);

        // write data blocks
        for (const auto& el: m_data)
        {
            const auto& name = el.first;
            std::uint64_t typeId = el.second.first;
            const auto& data = el.second.second;

            // write name
            std::uint64_t nameSize = name.size();
            std::fwrite(reinterpret_cast<const void*>(&nameSize), sizeof(nameSize), 1, pFile);
            std::fwrite(name.data(), sizeof(char), nameSize, pFile);

            // write type id
            std::fwrite(reinterpret_cast<const void*>(&typeId), sizeof(typeId), 1, pFile);

            // write data
            std::uint64_t dataSize = data.GetSize();
            std::fwrite(reinterpret_cast<const void*>(&dataSize), sizeof(dataSize), 1, pFile);
            std::fwrite(data.GetData(), 1, dataSize, pFile);
        }

        std::fclose(pFile);
        return true;
    }
    
    bool DataFile::LoadFromFile(const std::string& path)
    {
        std::FILE* pFile = std::fopen(path.c_str(), "rb");
        if (!pFile)
            return false;

        // read version
        std::uint64_t version[2];
        std::fread(version, sizeof(std::uint64_t), 2, pFile);

        // read block count
        std::uint64_t blockCnt = 0;
        std::fread(reinterpret_cast<void*>(&blockCnt), sizeof(blockCnt), 1, pFile);

        // read data blocks
        for (std::uint64_t i = 0; i < blockCnt; i++)
        {
            // read name
            std::uint64_t nameSize = 0;
            std::fread(reinterpret_cast<void*>(&nameSize), sizeof(nameSize), 1, pFile);
            std::string name(nameSize, 0);
            std::fread(&(name[0]), sizeof(char), nameSize, pFile);

            // read type id
            std::uint64_t typeId = 0;
            std::fread(reinterpret_cast<void*>(&typeId), sizeof(typeId), 1, pFile);

            // read data
            Memory data;
            std::uint64_t dataSize = 0;
            std::fread(reinterpret_cast<void*>(&dataSize), sizeof(dataSize), 1, pFile);
            data.Allocate(dataSize);
            std::fread(data.GetData(), 1, dataSize, pFile);

            m_data[name] = std::make_pair(typeId, std::move(data));
        }

        return true;
    }

    /*std::vector<char> CalcApp::SerializeMatrices() const
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

    std::map<std::string, TDynamicMatrix<double>> CalcApp::DeserializeMatrices(
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
    }*/

}
