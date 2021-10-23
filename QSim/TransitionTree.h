// Philipp Neufeld

#ifndef QSim_TransitionTree_H_
#define QSim_TransitionTree_H_

#include <cstdint>
#include <cstring>
#include <map>
#include <vector>

namespace QSim
{

    template<typename Ty>
    class TTransition
    {
    public:
        TTransition(std::size_t lvl1, std::size_t lvl2, Ty rabi);
            : m_l1(lvl1), m_l2(lvl2), m_rabi(rabi) { }

        // copy operations
        TTransition(const TTransition& rhs) = default;
        TTransition& operator=(const TTransition& rhs) = default;

        // Getter
        std::size_t GetLevel1Index() const { return m_l1; }
        std::size_t GetLevel2Index() const { return m_l2; }
        Ty GetRabi() const { return m_rabi; }

        // Setter
        void SetLevel1Index(std::size_t l1) const { m_l1 = l1; }
        void SetLevel2Index(std::size_t l2) const { m_l2 = l2; }
        void SetRabi(Ty rabi) { m_rabi = rabi; }

        // Comparison operators that enable this class to be used as keys in a map
        bool operator==(const TTransition& rhs) { return std::memcmp(this, &rhs, sizeof(TTransition)) == 0; }
        bool operator<(const TTransition& rhs) { return std::memcmp(this, &rhs, sizeof(TTransition)) < 0; }
        bool operator>(const TTransition& rhs) { return std::memcmp(this, &rhs, sizeof(TTransition)) > 0; }
        bool operator>=(const TTransition& rhs) { return !((*this) < rhs); }
        bool operator<=(const TTransition& rhs) { return !((*this) > rhs); }
        bool operator!=(const TTransition& rhs) { return !((*this) == rhs); }

    private:
        std::size_t m_l1;
        std::size_t m_l2;
        Ty m_rabi;
    }; 


    template<typename Ty>
    class TTransitionTreeNode
    {
    public:
        TTransitionTreeNode() : m_subNodes() { }

        // copy operators
        TTransitionTreeNode(const TTransitionTreeNode& rhs) = default;
        TTransitionTreeNode& operator=(const TTransitionTreeNode& rhs) = default;

        // node operations
        void AddNode(TTransition<Ty> trans) { m_subNodes.emplace(trans); }
        TTransitionTreeNode& GetNode(TTransition<Ty> trans) { return m_subNodes[trans]; }
        const TTransitionTreeNode& GetNode(TTransition<Ty> trans) const { return m_subNodes[trans]; }

    private:
        std::map<TTransition<Ty>, TTransitionTreeNode> m_subNodes;
    };


    template<typename Ty>
    class TTransitionTree
    {
    public:
        TTransitionTree(std::vector<TTransition<Ty>> transitions);

    private:

    };

}

#endif
