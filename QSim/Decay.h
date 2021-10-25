// Philipp Neufeld, 2021

#ifndef QSim_Decay_H_
#define QSim_Decay_H_

#include <cstdint>
#include <cstring>

namespace QSim
{

    template<typename Ty>
    class TDecay
    {
    public:
        TDecay(std::size_t initial, std::size_t final, Ty rate=0)
            : m_initial(initial), m_final(final), m_rate(rate) { }

        // copy operations
        TDecay(const TDecay& rhs) = default;
        TDecay& operator=(const TDecay& rhs) = default;

        // Getter
        std::size_t GetInitialIndex() const { return m_initial; }
        std::size_t GetFinalIndex() const { return m_final; }
        Ty GetRate() const { return m_rate; }

        // Setter
        void SetInitialIndex(std::size_t initial) { m_initial = initial; }
        void SetFinalIndex(std::size_t final) { m_final = final; }
        void SetRate(Ty rate) { m_rate = rate; }

        // Comparison operators that enable this class to be used as keys in a map
        bool operator==(const TDecay& rhs) const { return std::memcmp(this, &rhs, sizeof(TDecay)) == 0; }
        bool operator<(const TDecay& rhs) const { return std::memcmp(this, &rhs, sizeof(TDecay)) < 0; }
        bool operator>(const TDecay& rhs) const { return std::memcmp(this, &rhs, sizeof(TDecay)) > 0; }
        bool operator>=(const TDecay& rhs) const { return !((*this) < rhs); }
        bool operator<=(const TDecay& rhs) const { return !((*this) > rhs); }
        bool operator!=(const TDecay& rhs) const { return !((*this) == rhs); }

    private:
        std::size_t m_initial;
        std::size_t m_final;
        Ty m_rate;
    }; 

}

#endif
