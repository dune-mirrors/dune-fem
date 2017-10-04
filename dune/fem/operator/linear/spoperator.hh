#ifndef DUNE_FEM_SPOPERATOR_HH
#define DUNE_FEM_SPOPERATOR_HH

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! SparseRowLinearOperator
    template< class DomainFunction, class RangeFunction >
    struct SparseRowLinearOperator
    : public SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef SparseRowLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      typedef typename BaseType::MatrixType::field_type FieldType;

      static constexpr bool assembled = true;

      using BaseType::apply;
      using BaseType::matrix;

      SparseRowLinearOperator( const std::string & ,
                               const DomainSpaceType &domainSpace,
                               const RangeSpaceType &rangeSpace,
                               const MatrixParameter& param = SparseRowMatrixParameter() ) :
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

      void communicate()
      {}

      void maskRows ( const RangeFunction &maskFunction, FieldType diagonal = FieldType( 0 ) )
      {
        const auto &slaveDofs = maskFunction.space().slaveDofs();
        for( typename BaseType::size_type i = 0; i < matrix().rows(); ++i )
        {
          if( maskFunction.dofVector()[ i ] != FieldType( 0 ) )
            continue;

          matrix().clearRow( i );
          if( !slaveDofs.isSlave( i ) )
            matrix().set( i, i, diagonal );
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPOPERATOR_HH
