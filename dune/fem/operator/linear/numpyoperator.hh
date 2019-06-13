#ifndef DUNE_FEMPY_NUMPYOPERATOR_HH
#define DUNE_FEMPY_NUMPYOPERATOR_HH

// dune-fem includes
#include <dune/fem/operator/linear/spoperator.hh>

// local includes
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace Fem
  {

    //! NumpyLinearOperator
    template< class DomainFunction, class RangeFunction >
    struct NumpyLinearOperator
    : public SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType,
                                    typename RangeFunction::DiscreteFunctionSpaceType,
                                    SparseRowMatrix< double, size_t,
                                                     pybind11::array_t< double >,
                                                     pybind11::array_t<size_t> > >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef NumpyLinearOperator< DomainFunction, RangeFunction > ThisType;
      // for numpy backend we need to use different storage classes
      typedef SparseRowMatrix< double, size_t, pybind11::array_t< double >, pybind11::array_t<size_t> > Matrix;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, Matrix > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;

      NumpyLinearOperator( const std::string & ,
                           const DomainSpaceType &domainSpace,
                           const RangeSpaceType &rangeSpace,
                           const SolverParameter& param = SolverParameter() ) :
        BaseType( domainSpace, rangeSpace, param )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      const BaseType &systemMatrix() const
      {
        return *this;
      }

      BaseType &systemMatrix()
      {
        return *this;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPOPERATOR_HH
