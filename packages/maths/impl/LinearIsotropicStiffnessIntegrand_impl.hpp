#ifndef CIE_FEM_LINEAR_ISOTROPIC_STIFFNESS_INTEGRAND_IMPL_HPP
#define CIE_FEM_LINEAR_ISOTROPIC_STIFFNESS_INTEGRAND_IMPL_HPP

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/maths/inc/LinearIsotropicStiffnessIntegrand.hpp"


namespace cie::fem::maths {


template <Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand()
    : LinearIsotropicStiffnessIntegrand(0, nullptr)
{
}


template <Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const std::shared_ptr<const TAnsatzDerivatives>> pAnsatzDerivatives)
    : _modulus(modulus),
      _pAnsatzDerivatives(pAnsatzDerivatives),
      _buffer()
{
}


template <Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const std::shared_ptr<const TAnsatzDerivatives>> pAnsatzDerivatives,
                                                                                         std::span<Value> buffer)
    : LinearIsotropicStiffnessIntegrand(modulus, pAnsatzDerivatives)
{
    this->setBuffer(buffer);
}


template <Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::evaluate(ConstIterator itArgumentBegin,
                                                                     ConstIterator itArgumentEnd,
                                                                     Iterator itOut) const
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())
    CIE_CHECK_POINTER(_pAnsatzDerivatives)

    Ref<const TAnsatzDerivatives> rAnsatzDerivatives = *_pAnsatzDerivatives;
    const unsigned derivativeComponentCount = rAnsatzDerivatives.size();
    const unsigned ansatzCount = derivativeComponentCount / Dimension;
    rAnsatzDerivatives.evaluate(itArgumentBegin, itArgumentEnd, _buffer.data());

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    EigenAdaptor derivativeAdaptor(_buffer.data(), Dimension, ansatzCount);
    EigenAdaptor outputAdaptor(itOut, ansatzCount, ansatzCount);

    outputAdaptor = derivativeAdaptor.transpose() * _modulus * derivativeAdaptor;
}


template <Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::size() const
{
    const auto derivativeComponentCount = _pAnsatzDerivatives->size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::setBuffer(std::span<Value> buffer)
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
}


template <Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::getMinBufferSize() const noexcept
{
    return _pAnsatzDerivatives->size();
}


} // namespace cie::fem::maths


#endif
