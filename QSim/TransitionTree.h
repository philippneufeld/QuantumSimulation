// Philipp Neufeld, 2021

#ifndef QSim_TransitionTree_H_
#define QSim_TransitionTree_H_

#include <cstdint>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>

#include "Transition.h"

namespace QSim
{

    template<typename Ty>
    class TTransitionTree
    {
    public:
        TTransitionTree(std::size_t headLvlIdx) 
            : m_headLevelIdx(headLvlIdx), m_subNodes() { }

        // copy operators
        TTransitionTree(const TTransitionTree& rhs) = default;
        TTransitionTree(TTransitionTree&& rhs) = default;
        TTransitionTree& operator=(const TTransitionTree& rhs) = default;
        TTransitionTree& operator=(TTransitionTree&& rhs) = default;

        bool BuildTree(const std::vector<TTransition<Ty>>& transitions);

        // Getter
        const std::map<TTransition<Ty>, TTransitionTree>& GetNodes() const { return m_subNodes; }
        std::set<std::size_t> GetTreeLevelIndices() const { return m_levelIndices; }

        //
        template<typename VT, typename MT>
        bool AddPhotonBasis(const std::vector<TTransition<Ty>>& allTrans, const TMatrix<VT>& levels, TMatrix<MT>& inout) const;

    private:
        std::size_t m_headLevelIdx;
        std::map<TTransition<Ty>, TTransitionTree> m_subNodes;
        std::set<std::size_t> m_levelIndices;
    };


    template<typename Ty>
    bool TTransitionTree<Ty>::BuildTree(const std::vector<TTransition<Ty>>& transitions)
    {
        m_subNodes.clear();
        m_levelIndices.clear();

        m_levelIndices.insert(m_headLevelIdx);

        for (const auto& trans: transitions)
        {
            auto lvl1 = trans.GetLevel1Index();
            auto lvl2 = trans.GetLevel2Index();
            if (m_headLevelIdx == lvl1 || m_headLevelIdx == lvl2)
            {
                // find all remaining transitions (excluding the current one)
                std::vector<TTransition<Ty>> nextTransitions;
                std::copy_if(
                    transitions.begin(), 
                    transitions.end(), 
                    std::back_inserter(nextTransitions),
                    [&](auto t) { return t != trans; });

                // Build next tree level
                auto nextLevelIdx = (m_headLevelIdx == lvl1) ? lvl2 : lvl1;
                TTransitionTree<Ty> node(nextLevelIdx);
                bool success = node.BuildTree(nextTransitions);

                // check for circular transition paths
                if (success)
                    success = (node.m_levelIndices.find(m_headLevelIdx) == node.m_levelIndices.end());

                // If any of the above operations failed, abort
                if (!success)
                {
                    m_subNodes.clear();
                    m_levelIndices.clear();
                    return false;
                }

                // update members
                m_subNodes.insert(std::make_pair(trans, node));
                m_levelIndices.insert(node.m_levelIndices.begin(), node.m_levelIndices.end());
            }
        }
        return true;
    }

    template<typename Ty>
    template<typename VT, typename MT>
    bool TTransitionTree<Ty>::AddPhotonBasis(const std::vector<TTransition<Ty>>& allTrans, const TMatrix<VT>& levels, TMatrix<MT>& inout) const
    {
        // Validate matrices
        if ((~inout).Rows() != (~levels).Rows())
            return false;
        if ((~inout).Cols() != allTrans.size())
            return false;
        if ((~levels).Cols() != 1)
            return false;

        for (const auto& el: m_subNodes)
        {
            const auto& trans = el.first;
            auto it = std::find(allTrans.begin(), allTrans.end(), trans);
            if (it == allTrans.end())
                return false;
            std::size_t transIdx = it - allTrans.begin();

            Ty rel = 1;
            if ((~levels)(m_headLevelIdx, 0) < std::max((~levels)(trans.GetLevel1Index(), 0), (~levels)(trans.GetLevel2Index(), 0)))
                rel = -rel;

            const auto& subtree = el.second;
            for (std::size_t lvlIdx: subtree.m_levelIndices)
                (~inout)(lvlIdx, transIdx) += rel;
            
            if (!subtree.AddPhotonBasis(allTrans, levels, inout))
                return false;
        }

        return true;
    }

}

#endif
