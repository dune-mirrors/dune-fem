#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/projection/local/l2projection.hh>
#include <dune/fem/operator/projection/local/riesz.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridPart, class BasisFunctionSet,
              class Quadrature = CachingQuadrature< GridPart, BasisFunctionSet::EntityType::codimension > >
    class DiscontinuousGalerkinLocalL2Projection;



    // DiscontinuousGalerkinLocalL2Projection
    // --------------------------------------

    template< class GridPart, class BasisFunctionSet, class Quadrature >
    class DiscontinuousGalerkinLocalL2Projection
    : public LocalL2Projection< BasisFunctionSet, DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > >
    {
      typedef DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > ThisType;
      typedef LocalL2Projection< BasisFunctionSet, DiscontinuousGalerkinLocalL2Projection< GridPart, BasisFunctionSet, Quadrature > > BaseType;

    public:
      /** \copydoc Dune::Fem::LocalL2Projection::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    private:
      typedef GridPart GridPartType;
      typedef typename GridPartType::GridType GridType;
      static const bool cartesian = true;//Dune::Capabilities::isCartesian< GridType >::v;

      typedef typename std::conditional< cartesian,
          OrthonormalLocalRieszProjection< BasisFunctionSetType >,
          DenseLocalRieszProjection< BasisFunctionSetType, Quadrature >
        >::type LocalRieszProjectionType;

      typedef DefaultLocalL2Projection< LocalRieszProjectionType, Quadrature > Implementation;

    public:
      /** \name Construction
       *  \{
       */

      explicit DiscontinuousGalerkinLocalL2Projection ( const BasisFunctionSetType &basisFunctionSet )
        : impl_( LocalRieszProjectionType( basisFunctionSet ) )
      {}

      explicit DiscontinuousGalerkinLocalL2Projection ( BasisFunctionSetType &&basisFunctionSet )
        : impl_( LocalRieszProjectionType( std::forward< BasisFunctionSetType >( basisFunctionSet ) ) )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      DiscontinuousGalerkinLocalL2Projection ( const ThisType & ) = default;

      DiscontinuousGalerkinLocalL2Projection ( ThisType &&other )
        : impl_( std::move( other.impl_ ) )
      {}

      ThisType &operator= ( const ThisType & ) = default;

      ThisType &operator= ( ThisType &&other )
      {
        impl_ = std::move( other.impl_ );
        return *this;
      }

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::LocalL2Projection::basisFunctionSet */
      BasisFunctionSet basisFunctionSet () const
      {
        return impl_.basisFunctionSet();
      }

      /** \copydoc Dune::Fem::LocalL2Projection::apply */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        impl_( localFunction, localDofVector );
      }

      /** \} */

    private:
      Implementation impl_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_INTERPOLATION_HH
