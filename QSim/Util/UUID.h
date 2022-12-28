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
        bool operator!=(const UUIDv4& rhs) const;
        bool operator<(const UUIDv4& rhs) const;
        bool operator<=(const UUIDv4& rhs) const;
        bool operator>(const UUIDv4& rhs) const;
        bool operator>=(const UUIDv4& rhs) const;

        std::size_t Hash() const;

        void ToString(char* buffer) const;
        void ToString(std::string& buffer) const;
        std::string ToString() const;

        static UUIDv4 FromString(const char* uuidStr);
        static UUIDv4 FromString(const std::string& uuidStr);

        void StoreToBufferLE(void* buffer) const;
        void StoreToBufferBE(void* buffer) const;

    private:
        __m128i Load() const;
        void Store(__m128i uuid);

        static __m128i Generate();
        static __m128i PackUUID(__m128i hi, __m128i lo);
        static __m128i SwapEndianess(__m128i v);

    private:
        std::uint8_t m_data[16];
    };

}

namespace std
{
    template<>
    struct hash<QSim::UUIDv4>
    {
        std::size_t operator()(const QSim::UUIDv4& uuid) const
        {
            return uuid.Hash();
        }
    };
}

#endif
