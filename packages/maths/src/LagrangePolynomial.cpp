// --- FEM Includes ---
#include "packages/maths/inc/LagrangePolynomial.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/NChooseK.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <iterator>


namespace cie::fem::maths {


template <class TValue>
LagrangePolynomial<TValue>::LagrangePolynomial(Ptr<const TValue> pNodeBegin,
                                               Ptr<const TValue> pNodeEnd,
                                               Size iBase)
    : Polynomial<TValue>()
{
    CIE_BEGIN_EXCEPTION_TRACING

    const Size polynomialOrder = std::distance(pNodeBegin, pNodeEnd);

    if (polynomialOrder) [[likely]] {
        DynamicArray<TValue> localNodes;
        localNodes.reserve(polynomialOrder - 1);
        this->_coefficients.resize(polynomialOrder);

        // Get the base node
        CIE_OUT_OF_RANGE_CHECK(iBase < polynomialOrder)
        const auto itBaseNode = pNodeBegin + iBase;
        const TValue baseNode = *itBaseNode;

        // Collect and negate the rest.
        // At the same time, compute the denominator
        TValue denominator = static_cast<TValue>(1);
        const auto transform = [&denominator, baseNode](TValue node) {
            denominator *= baseNode - node;
            return -node;
        };
        CIE_DIVISION_BY_ZERO_CHECK(denominator)

        auto itLocal = std::back_inserter(localNodes);
        std::transform(pNodeBegin, itBaseNode, itLocal, transform);

        if (iBase < polynomialOrder - 1)
            std::transform(itBaseNode + 1, pNodeEnd, itLocal, transform);

        // Loop backward through the exponents
        const auto itLocalNodeBegin = localNodes.begin();
        for (Size numberOfSelectedNodes=0; numberOfSelectedNodes<polynomialOrder; ++numberOfSelectedNodes) {
            const Size exponent = polynomialOrder - 1 - numberOfSelectedNodes;
            utils::NChooseK permutation(polynomialOrder - 1, numberOfSelectedNodes);
            TValue coefficient = 0;
            do {
                TValue term = static_cast<TValue>(1);
                for (auto iNode : permutation)
                    term *= *(itLocalNodeBegin + iNode);
                coefficient += term;
            } while (permutation.next());

            this->_coefficients[exponent] = coefficient / denominator;
        } // for numberOfSelectedNodes in range(polynomialOrder)
    } // if polynomialOrder

    CIE_END_EXCEPTION_TRACING
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(LagrangePolynomial);


} // namespace cie::fem::maths
