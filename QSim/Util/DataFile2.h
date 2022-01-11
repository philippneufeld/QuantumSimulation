// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_DataFile2_H_
#define QSim_Util_DataFile2_H_

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>


namespace QSim
{
    
    class DataFileOSHandle
    {
    public:
        DataFileOSHandle(std::FILE* pFile);
        ~DataFileOSHandle();

        DataFileOSHandle(const DataFileOSHandle&) = delete;
        DataFileOSHandle(DataFileOSHandle&& rhs);

        DataFileOSHandle& operator=(const DataFileOSHandle&) = delete;
        DataFileOSHandle& operator=(DataFileOSHandle&& rhs);

        std::FILE* GetNative();

    private:
        std::FILE* m_pFile;
    };

    class DataFileDriver
    {
    public:
        DataFileDriver();

    private:
        std::shared_ptr<DataFileOSHandle> m_pFile;
    };

}

#endif
