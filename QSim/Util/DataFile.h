// Philipp Neufeld, 2021

#ifndef QSim_Util_DataFile_H_
#define QSim_Util_DataFile_H_

#include <string>
#include <map>
#include <vector>
#include <cstdint>
#include <memory>

#include "Memory.h"

namespace QSim
{
    
    class DataFile
    {
        constexpr static std::uint64_t MajorVersion_v = 1;
        constexpr static std::uint64_t MinorVersion_v = 0;
    public:
        DataFile();
        ~DataFile();

        // Data control
        void Clear();
        void SetData(const std::string& name, const void* data, std::size_t n);
        void SetDataEx(const std::string& name, const void* data, std::size_t n, std::uint64_t typeId);
        const Memory& GetData(const std::string& name) const;
        std::uint64_t GetDataTypeId(const std::string& name) const;
        void RemoveData(const std::string& name);
        bool Contains(const std::string& name) const;
        std::vector<std::string> GetIndex() const;
        
        // File I/O operations
        bool SaveToFile(const std::string& path) const;
        bool LoadFromFile(const std::string& path);

    private:
        std::map<std::string, std::pair<std::uint64_t, Memory>> m_data;
    };

}

#endif
