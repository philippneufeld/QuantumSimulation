// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_Memory
#define QSim_Util_Memory

#include <cstdint>

#include "../Platform.h"

namespace QSim
{

    class Memory
    {
    public:
        Memory() : m_size(0), m_data(nullptr) { }
        Memory(const void* data, std::size_t n);
        ~Memory() { Release(); }

        // copy operations
        Memory(const Memory& rhs);
        Memory(Memory&& rhs);
        Memory& operator=(const Memory& rhs);
        Memory& operator=(Memory&& rhs);

        bool Allocate(std::size_t size);
        void Release();

        void* GetData() { return m_data; }
        const void* GetData() const { return m_data; }
        std::size_t GetSize() const { return m_size; }

    private:
        std::size_t m_size;
        std::uint8_t* m_data;
    };

}

#endif
