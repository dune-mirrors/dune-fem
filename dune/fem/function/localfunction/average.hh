#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_AVERAGE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_AVERAGE_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/discontinuousgalerkin/declaration.hh>
#include <dune/fem/space/finitevolume/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class LocalFunction, class GridPart >
    class LocalAverage;



    // LocalAverageHelper
    // ------------------

    struct LocalAverageHelper
    {
      template< class LocalFunction, class Quadrature >
      static void applyQuadrature ( const LocalFunction &localFunction,
                                    const typename LocalFunction::EntityType::Geometry &geometry,
                                    const Quadrature &quadrature,
                                    typename LocalFunction::RangeType &value )
      {
        typedef typename LocalFunction::RangeType RangeType;
        typedef typename RangeType::value_type RangeFieldType;

        value = RangeType( 0 );
        RangeFieldType volume( 0 );
        const int nop = quadrature.nop();
        for( int qp = 0; qp < nop; ++qp )
        {
          RangeType tmp;
          localFunction.evaluate( quadrature[ qp ], tmp );
          RangeFieldType weight = quadrature.weight( qp )*geometry.integrationElement( quadrature.point( qp ) );
          volume += weight;
          value.axpy( weight, tmp );
        }
        value /= volume;
      }
    };



    // LocalAverageImpl
    // ----------------

    template< class LocalFunction, class GridPart, class DiscreteFunctionSpace >
    struct LocalAverageImpl
    {
      static void apply ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average )
      {
        // get entity
        typedef typename LocalFunction::EntityType EntityType;
        const EntityType &entity = localFunction.entity();

        // create quadrature
        typedef CachingQuadrature< GridPart, EntityType::codimension > QuadratureType;
        QuadratureType quadrature( entity, localFunction.order() );
        LocalAverageHelper::applyQuadrature( localFunction, entity.geometry(), quadrature, average );
      }
    };



    // LocalAverageImpl for DiscontinuousGalerkinSpace
    // -----------------------------------------------

    template< class LocalFunction, class GridPart, class FunctionSpace, int order, class Storage >
    class LocalAverageImpl< LocalFunction, GridPart, DiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage > >
    {
      static const bool cartesian = GridPartCapabilities::isCartesian< GridPart >::v;

      static void applyAffine ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average )
      {
        // otherwise use direct computation
        const int dimRange = LocalFunction::dimRange;
        for( int i = 0; i < dimRange; ++i )
          average[ i ] = localFunction[ i ];

        const typename LocalFunction::BasisFunctionSetType &basisFunctionSet = localFunction.basisFunctionSet();
        average /= std::sqrt( basisFunctionSet.referenceElement().volume() );
      }

    public:
      static void apply ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average )
      {
        if( cartesian )
          return applyAffine( localFunction, average );

        // get entity and geometry
        typedef typename LocalFunction::EntityType EntityType;
        const EntityType &entity = localFunction.entity();
        const typename EntityType::Geometry geometry = entity.geometry();

        if( geometry.affine() )
          return applyAffine( localFunction, average );

        // for not affine elements, use quadrature
        typedef CachingQuadrature< GridPart, EntityType::codimension > QuadratureType;
        QuadratureType quadrature( entity, localFunction.order() );
        LocalAverageHelper::applyQuadrature( localFunction, geometry, quadrature, average );
      }
    };



    // LocalAverageImpl for FiniteVolumeSpace
    // --------------------------------------

    template< class LocalFunction, class GridPart, class FunctionSpace, int codim, class Storage >
    class LocalAverageImpl< LocalFunction, GridPart, FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
    public:
      static void apply ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average )
      {
        localFunction.evaluate( typename LocalFunction::DomainType( 0 ), average );
      }
    };



    // LocalAverage
    // ------------

    template< class LocalFunction, class GridPart >
    class LocalAverage
    {
      typedef typename ExportsDiscreteFunctionSpaceType< LocalFunction >::Type DiscreteFunctionSpaceType;
      typedef LocalAverageImpl< LocalFunction, GridPart, DiscreteFunctionSpaceType > Impl;

    public:
      static void apply ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average )
      {
        Impl::apply( localFunction, average );
      }

      void operator() ( const LocalFunction &localFunction, typename LocalFunction::RangeType &average ) const
      {
        Impl::apply( localFunction, average );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_AVERAGE_HH
