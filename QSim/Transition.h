// Philipp Neufeld, 2021

#ifndef QSim_Transition_H_
#define QSim_Transition_H_

#include <cstdint>
#include <cstring>

namespace QSim
{

    template<typename Ty>
    class TTransition
    {
    public:
        TTransition(std::size_t lvl1, std::size_t lvl2, Ty rabi=0)
            : m_l1(lvl1), m_l2(lvl2), m_rabi(rabi) { }

        // copy operations
        TTransition(const TTransition& rhs) = default;
        TTransition& operator=(const TTransition& rhs) = default;

        // Getter
        std::size_t GetLevel1Index() const { return m_l1; }
        std::size_t GetLevel2Index() const { return m_l2; }
        Ty GetRabi() const { return m_rabi; }

        // Setter
        void SetLevel1Index(std::size_t l1) { m_l1 = l1; }
        void SetLevel2Index(std::size_t l2) { m_l2 = l2; }
        void SetRabi(Ty rabi) { m_rabi = rabi; }

        // Comparison operators that enable this class to be used as keys in a map
        bool operator==(const TTransition& rhs) const { return std::memcmp(this, &rhs, sizeof(TTransition)) == 0; }
        bool operator<(const TTransition& rhs) const { return std::memcmp(this, &rhs, sizeof(TTransition)) < 0; }
        bool operator>(const TTransition& rhs) const { return std::memcmp(this, &rhs, sizeof(TTransition)) > 0; }
        bool operator>=(const TTransition& rhs) const { return !((*this) < rhs); }
        bool operator<=(const TTransition& rhs) const { return !((*this) > rhs); }
        bool operator!=(const TTransition& rhs) const { return !((*this) == rhs); }

    private:
        std::size_t m_l1;
        std::size_t m_l2;
        Ty m_rabi;
    }; 

}

#endif
