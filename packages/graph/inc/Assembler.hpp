#ifndef CIE_FEM_ASSEMBLER_HPP
#define CIE_FEM_ASSEMBLER_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/maths/inc/Comparison.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <packages/macros/inc/checks.hpp>
#include <unordered_map>
#include <optional> // optional
#include <ranges> // transform_view
#include <set> // set


namespace cie::fem {


class Assembler
{
private:
    using DefaultMesh = Graph<void,void>;

    using DoFMap = std::unordered_map<
        DefaultMesh::VertexID,                      // <== ID of the cell
        DynamicArray<std::optional<std::size_t>>    // <== DoF indices of its basis functions
    >;

    static auto makeIndexView(Ref<const DoFMap::mapped_type> rIndices)
    {
        return std::ranges::transform_view(
            rIndices,
            [](DoFMap::mapped_type::value_type maybeIndex){return maybeIndex.value();}
        );
    }

public:
    using DoFPairVector = DynamicArray<std::pair<std::size_t,std::size_t>>;

    using DoFPairIterator = std::back_insert_iterator<DoFPairVector>;

    Assembler() noexcept;

    explicit Assembler(std::size_t dofBegin) noexcept;

    template <class TVertexData,
              class TEdgeData,
              concepts::FunctionWithSignature<std::size_t,Ref<const typename Graph<TVertexData,TEdgeData>::Vertex>> TDoFCounter,
              concepts::FunctionWithSignature<void,Ref<const typename Graph<TVertexData,TEdgeData>::Edge>,DoFPairIterator> TDoFPairFunctor>
    void addGraph(Ref<const Graph<TVertexData,TEdgeData>> rGraph,
                  TDoFCounter&& rDoFCounter,
                  TDoFPairFunctor&& rDoFMatcher);

    std::size_t dofCount() const noexcept;

    template <class TIndex, class TValue>
    void makeCSRMatrix(Ref<TIndex> rRowCount,
                       Ref<TIndex> rColumnCount,
                       Ref<DynamicArray<TIndex>> rRowExtents,
                       Ref<DynamicArray<TIndex>> rColumnIndices,
                       Ref<DynamicArray<TValue>> rNonzeros,
                       OptionalRef<mp::ThreadPoolBase> rThreadPool = {}) const;

    auto keys() const
    {
        return std::ranges::views::keys(_dofMap);
    }

    auto values() const
    {
        return std::ranges::transform_view(
            _dofMap,
            [](const auto& rPair) {return Assembler::makeIndexView(rPair.second);}
        );
    }

    auto items() const
    {
        return std::ranges::transform_view(
            _dofMap,
            [](const auto& rPair) {
                return std::make_pair(rPair.first,
                                      Assembler::makeIndexView(rPair.second));
            }
        );
    }

    auto operator[](DefaultMesh::VertexID vertexID) const
    {
        const auto it = _dofMap.find(vertexID);
        CIE_OUT_OF_RANGE_CHECK(it != _dofMap.end())
        return this->makeIndexView(it->second);
    }

private:
    std::size_t _dofCounter;

    DoFMap _dofMap;

    using IndexPair = StaticArray<std::size_t,2>;
    std::set<IndexPair,utils::Ordering<IndexPair>> _sparsityPattern;
}; // class Assembler


} // namespace cie::fem

#include "packages/graph/impl/Assembler_impl.hpp"

#endif
