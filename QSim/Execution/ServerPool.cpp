// Philipp Neufeld, 2021-2022

#include <algorithm>
#include "ServerPool.h"

namespace QSim
{
    //
    // SocketDataPackageBin
    //

    SocketDataPackageBin::SocketDataPackageBin() 
        : m_size(0), m_pData(nullptr) {}
    
    SocketDataPackageBin::~SocketDataPackageBin() 
    { 
        Allocate(0);
    }

    SocketDataPackageBin::SocketDataPackageBin(const SocketDataPackageBin& rhs)
        : SocketDataPackageBin()
    {
        if (rhs.m_size > 0)
        {
            Allocate(rhs.m_size);
            std::copy(rhs.m_pData, rhs.m_pData+rhs.m_size, m_pData);
        }
    }
    
    SocketDataPackageBin::SocketDataPackageBin(SocketDataPackageBin&& rhs)
        : m_pData(rhs.m_pData), m_size(rhs.m_size)
    {
        rhs.m_pData = nullptr;
        rhs.m_size = 0;
    }

    SocketDataPackageBin& SocketDataPackageBin::operator=(const SocketDataPackageBin& rhs)
    {
        SocketDataPackageBin tmp(rhs);
        std::swap(*this, tmp);
    }

    SocketDataPackageBin& SocketDataPackageBin::operator=(SocketDataPackageBin&& rhs)
    {
        std::swap(m_pData, rhs.m_pData);
        std::swap(m_size, rhs.m_size);
    }
    
    bool SocketDataPackageBin::Allocate(std::uint64_t size)
    {
        if (size == m_size)
            return true;

        if (m_pData)
            delete[] m_pData;
        m_pData = nullptr;
        m_size = 0;      
        
        if (size > 0)
        {
            m_pData = new std::uint8_t[size];
            if (!m_pData)
            {
                m_pData = nullptr;
                m_size = 0;
            }
        }

        return (size == m_size);
    }


}
