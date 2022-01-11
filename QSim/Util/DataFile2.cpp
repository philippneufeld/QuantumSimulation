// Philipp Neufeld, 2021-2022

#include "DataFile2.h"

namespace QSim
{
    DataFileOSHandle::DataFileOSHandle(std::FILE* pFile)
        : m_pFile(pFile) { }

    DataFileOSHandle::~DataFileOSHandle()
    {
        if(m_pFile)
            std::fclose(m_pFile);
    }

    DataFileOSHandle::DataFileOSHandle(DataFileOSHandle&& rhs)
        : m_pFile(rhs.m_pFile)
    {
        rhs.m_pFile = nullptr;
    }

    DataFileOSHandle& DataFileOSHandle::operator=(DataFileOSHandle&& rhs)
    {
        std::swap(m_pFile, rhs.m_pFile);
        return *this;
    }

    std::FILE* DataFileOSHandle::GetNative() 
    { 
        return m_pFile; 
    }



    DataFileDriver::DataFileDriver()
        : m_pFile(nullptr)
    {

    }
}
