#ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_DENSE_HH
#define DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_DENSE_HH

#include <cassert>
#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include "localrieszprojection.hh"

namespace Dune
{

  namespace Fem
  {

    // DenseLocalRieszProjection
    // -------------------------

    template< class BasisFunctionSet, class Quadrature >
    class DenseLocalRieszProjection
    : public LocalRieszProjection< BasisFunctionSet, DenseLocalRieszProjection< BasisFunctionSet, Quadrature > >
    {
      typedef DenseLocalRieszProjection< BasisFunctionSet, Quadrature > ThisType;
      typedef LocalRieszProjection< BasisFunctionSet, DenseLocalRieszProjection< BasisFunctionSet, Quadrature > > BaseType;

    public:
      /** \copydoc Dune::Fem::LocalRieszProjection::BasisFunctionSet */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      /** \name Construction
       *  \{
       */

      explicit DenseLocalRieszProjection ( const BasisFunctionSetType &basisFunctionSet )
        : basisFunctionSet_( basisFunctionSet ),
          inverse_( basisFunctionSet.size(), basisFunctionSet.size() )
      {
        assemble( basisFunctionSet, inverse_ );
        inverse_.invert();
      }

      explicit DenseLocalRieszProjection ( BasisFunctionSetType &&basisFunctionSet )
        : basisFunctionSet_( std::forward< BasisFunctionSetType >( basisFunctionSet ) ),
          inverse_( basisFunctionSet.size(), basisFunctionSet.size() )
      {
        assemble( basisFunctionSet, inverse_ );
        inverse_.invert();
      }

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      DenseLocalRieszProjection ( const ThisType &other ) = default;

      DenseLocalRieszProjection ( ThisType &&other )
        : basisFunctionSet_( std::move( other.basisFunctionSet_ ) ),
          inverse_( std::move( other.inverse_ ) )
      {}

      ThisType &operator= ( const ThisType &other ) = default;

      ThisType &operator= ( ThisType &&other )
      {
        basisFunctionSet_ = std::move( other.basisFunctionSet_ );
        inverse_ = std::move( other.inverse_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::LocalRieszProjection::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet () const { return basisFunctionSet_; }

      /** \copydoc Dune::Fem::LocalRieszProjection::apply */
      template< class F, class LocalDofVector >
      void apply ( const F&f, LocalDofVector &dofs ) const
      {
        inverse_.mv( f, dofs );
      }

      /** \} */

    private:
      template< class DenseMatrix >
      void assemble ( const BasisFunctionSetType &basisFunctionSet, DenseMatrix &matrix )
      {
        const std::size_t size = basisFunctionSet.size();
        std::vector< typename BasisFunctionSetType::RangeType > phi( size );

        matrix = static_cast< typename DenseMatrix::value_type >( 0 );
        assert( matrix.N() == basisFunctionSet.size() );
        assert( matrix.M() == basisFunctionSet.size() );

        const auto &geometry = basisFunctionSet.entity().geometry();
        const Quadrature quadrature( geometry.type(), 2*basisFunctionSet.order() );
        const std::size_t nop = quadrature.nop();
        for( std::size_t qp = 0; qp < nop; ++qp )
        {
          basisFunctionSet.evaluateAll( quadrature[ qp ], phi );

          const auto weight = quadrature.weight( qp )*geometry.integrationElement( quadrature.point( qp ) );
          for( std::size_t i = 0; i < size; ++i )
          {
            matrix[ i ][ i ] += weight*( phi[ i ] * phi[ i ] );
            for( std::size_t j = i+1; j < size; ++j )
            {
              const auto value = weight*( phi[ i ]* phi[ j ] );
              matrix[ i ][ j ] += value;
              matrix[ j ][ i ] += value;
            }
          }
        }
      }

      BasisFunctionSetType basisFunctionSet_;
      Dune::DynamicMatrix< typename BasisFunctionSetType::FunctionSpaceType::RangeFieldType > inverse_;
    };

  } // namespace Fem

} // namepsace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_DENSE_HH

