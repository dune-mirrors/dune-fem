#ifndef DUNE_FEM_SPOPERATOR_HH
#define DUNE_FEM_SPOPERATOR_HH

// system includes
#include <string>

// local includes
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! SparseRowLinearOperator
    template< class DomainFunction, class RangeFunction,
              class Matrix = SparseRowMatrix< typename DomainFunction::DiscreteFunctionSpaceType::RangeFieldType > >
    struct SparseRowLinearOperator
    : public SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType,
                                    typename RangeFunction::DiscreteFunctionSpaceType,
                                    Matrix >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef SparseRowLinearOperator< DomainFunction, RangeFunction, Matrix > ThisType;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, Matrix > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;
      using BaseType::exportMatrix;

      SparseRowLinearOperator( const std::string & ,
                               const DomainSpaceType &domainSpace,
                               const RangeSpaceType &rangeSpace,
                               const SolverParameter& param = SolverParameter() ) :
        BaseType( domainSpace, rangeSpace, param )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      virtual void finalize () { BaseType::compress(); }

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPOPERATOR_HH
