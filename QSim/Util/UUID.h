// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_UUID_H_
#define QSim_Util_UUID_H_

#include <cstdint>
#include <string>

namespace QSim
{

    class UUIDv4
    {
    public:
        UUIDv4();

        UUIDv4(const UUIDv4& rhs) = default;
        UUIDv4& operator=(const UUIDv4& rhs) = default;

        void ToString(char* buffer) const;
        void ToString(std::string& buffer) const;
        std::string ToString() const;

    private:
        static void Generate(std::uint8_t* buffer);

    private:
        std::uint8_t m_data[16];
    };

}

#endif
