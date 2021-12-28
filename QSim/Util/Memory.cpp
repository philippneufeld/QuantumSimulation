// Philipp Neufeld, 2021-2022

#include <cstdlib>
#include <algorithm>

#include "Memory.h"

namespace QSim
{

    Memory::Memory(const void* data, std::size_t n) : Memory()
    {
        Allocate(n);
        std::copy_n(static_cast<const std::uint8_t*>(data), n, m_data);
    }

    Memory::Memory(const Memory& rhs) : Memory()
    {
        Allocate(rhs.m_size);
        std::copy_n(rhs.m_data, rhs.m_size, m_data);
    }

    Memory::Memory(Memory&& rhs) : Memory()
    {
        std::swap(m_data, rhs.m_data);
        std::swap(m_size, rhs.m_size);
    }

    Memory& Memory::operator=(const Memory& rhs)
    {
        Allocate(rhs.m_size);
        std::copy_n(rhs.m_data, rhs.m_size, m_data);
        return *this;
    }

    Memory& Memory::operator=(Memory&& rhs)
    {
        std::swap(m_data, rhs.m_data);
        std::swap(m_size, rhs.m_size);
        return *this;
    }

    bool Memory::Allocate(std::size_t size)
    {
        if (size != m_size)
            Release();

        if (size > 0)
        {
            m_data = static_cast<std::uint8_t*>(std::malloc(size));
            m_size = size;
            return !!m_data;
        }
        else
            return true;
    }

    void Memory::Release()
    {
        if (m_data)
            std::free(m_data);
        m_data = nullptr;
        m_size = 0;
    }
}
