#ifndef DUNE_FEM_SPACE_LAGRANGE_INTERPOLATION_HH
#define DUNE_FEM_SPACE_LAGRANGE_INTERPOLATION_HH

#include <cstddef>

#include <utility>

#include "lagrangepoints.hh"

namespace Dune
{

  namespace Fem
  {

    // LagrangeLocalInterpolation
    // --------------------------

    template< class GridPart, int maxOrder, class BasisFunctionSet >
    class LagrangeLocalInterpolation
    {
      typedef LagrangeLocalInterpolation< GridPart, maxOrder, BasisFunctionSet > ThisType;

    public:
      /** \brief basis function set type */
      typedef BasisFunctionSet BasisFunctionSetType;
      /** \brief point set type */
      typedef LagrangePointSet< GridPart, maxOrder > LagrangePointSetType;

    private:
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

    public:
      /** \name Construction
       *  \{
       */

      LagrangeLocalInterpolation ()
        : pointSet_( nullptr )
        , basisFunctionSet_()
      {
      }

      LagrangeLocalInterpolation ( const LagrangePointSetType &pointSet,
                                   const BasisFunctionSetType &basisFunctionSet )
        : pointSet_( &pointSet ),
          basisFunctionSet_( basisFunctionSet )
      {}

      LagrangeLocalInterpolation ( const LagrangePointSetType &pointSet,
                                   BasisFunctionSetType &&basisFunctionSet )
        : pointSet_( &pointSet ),
          basisFunctionSet_( std::forward< BasisFunctionSetType >( basisFunctionSet ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      LagrangeLocalInterpolation ( const ThisType & ) = default;

      /** \brief move constructor */
      LagrangeLocalInterpolation ( ThisType &&other )
        : pointSet_( std::move( other.pointSet_ ) ),
          basisFunctionSet_( std::move( other.basisFunctionSet_ ) )
      {}

      /** \brief assignment operator */
      LagrangeLocalInterpolation &operator= ( const ThisType & ) = default;

      /** \brief move assignment operator */
      LagrangeLocalInterpolation &operator= ( ThisType &&other )
      {
        pointSet_ = other.pointSet_ ;
        basisFunctionSet_ = std::move( other.basisFunctionSet_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \brief return basis function set */
      BasisFunctionSetType basisFunctionSet () const
      {
        return basisFunctionSet_;
      }

      /** \brief apply interpolation */
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      /** \brief apply interpolation */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        const LagrangePointSetType &pointSet = this->pointSet();

        int k = 0;
        const std::size_t nop = pointSet.nop();
        for( std::size_t pt = 0; pt < nop; ++pt )
        {
          typename FunctionSpaceType::RangeType phi;
          localFunction.evaluate( pointSet[ pt ], phi );
          for( int i = 0; i < FunctionSpaceType::dimRange; ++i )
            localDofVector[ k++ ] = phi[ i ];
        }
      }

      /** \} */

      void unbind()
      {
        pointSet_ = nullptr;
        // basisFunctionSet_ = BasisFunctionSetType();
      }

    protected:
      const LagrangePointSetType &pointSet () const
      {
        assert( pointSet_ );
        return *pointSet_;
      }

      const LagrangePointSetType* pointSet_ = nullptr;
      BasisFunctionSetType basisFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_INTERPOLATION_HH
