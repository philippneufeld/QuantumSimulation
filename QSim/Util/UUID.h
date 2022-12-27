// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_UUID_H_
#define QSim_Util_UUID_H_

#include <cstdint>
#include <string>

// SSE and SSE2 are always available
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2

namespace QSim
{

    class UUIDv4
    {
    private:
        UUIDv4(__m128i uuid);

    public:
        UUIDv4();

        UUIDv4(const UUIDv4& rhs) = default;
        UUIDv4& operator=(const UUIDv4& rhs) = default;

        bool operator==(const UUIDv4& rhs) const;

        void ToString(char* buffer) const;
        void ToString(std::string& buffer) const;
        std::string ToString() const;

        static UUIDv4 FromString(const char* uuidStr);
        static UUIDv4 FromString(const std::string& uuidStr);

    private:
        static __m128i Generate();

    private:
        std::uint8_t m_data[16];
    };

}

#endif
