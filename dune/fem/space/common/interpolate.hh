#ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
#define DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH

#include <algorithm>
#include <type_traits>
#include <utility>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/bindguard.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/localinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------

    /**
     * \function interpolate
     * \ingroup  DiscreteFunctionSpace
     * \brief    perform native interpolation of a discrete function space
     *
     * By definition of its degrees of freedom, each discrete function space
     * has a native interpolation, which can be computed very quickly.
     *
     * For example, the native interpolation of a Lagrange discrete function
     * space is the evaluation in its Lagrange points.
     * An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
     *
     * The actual implementation must locally be provided by the discrete
     * function space through the method
     * \code
     * template< class LocalFunction, class LocalDofVector >
     * void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
     * \endcode
     *
     * \param[in]   u  grid function to interpolate
     * \param[out]  v  discrete function to represent the interpolation
     */
    template< class GridFunction, class DiscreteFunction >
    static inline void interpolate ( const GridFunction &u, DiscreteFunction &v )
    {
      // just call interpolate for the all partition
      interpolate( u, v, Partitions::all );
    }

    template< class Function, class DiscreteFunction, unsigned int partitions >
    static inline std::enable_if_t< !std::is_convertible< Function, HasLocalFunction >::value >
    interpolate ( const Function &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: GridPartType  GridPartType;
      typedef GridFunctionAdapter< Function, GridPartType > GridFunctionType;

      GridFunctionType uGrid( "uGrid", u, v.space().gridPart() );

      interpolate( uGrid, v, ps );
    }

    template< class GridFunction, class DiscreteFunction, unsigned int partitions >
    static inline std::enable_if_t< std::is_convertible< GridFunction, HasLocalFunction >::value && Capabilities::hasInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >::v >
    interpolate ( const GridFunction &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      ConstLocalFunction< GridFunction > uLocal( u );
      LocalContribution< DiscreteFunction, Assembly::Set > vLocal( v );
      LocalInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >
        interpolation( v.space() );

      // iterate over selected partition
      for( const auto& entity : elements( v.gridPart(), ps ) )
      {
        // initialize u to entity
        auto uGuard = bindGuard( uLocal, entity );

        // bind v to entity
        auto vGuard = bindGuard( vLocal, entity );

        // bind interpolation to entity
        auto iGuard = bindGuard( interpolation, entity );

        // perform local interpolation
        interpolation( uLocal, vLocal );
      }
    }



    namespace Impl
    {

      template< class Entity, class FunctionSpace, class Weight >
      struct WeightLocalFunction
      {
        typedef Entity EntityType;
        typedef FunctionSpace FunctionSpaceType;

        typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
        typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

        typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

        static constexpr int dimDomain = FunctionSpaceType::dimDomain;
        static constexpr int dimRange = FunctionSpaceType::dimRange;

        explicit WeightLocalFunction ( Weight &weight, int order ) : weight_( weight ), order_(order) {}

        void bind ( const EntityType &entity ) { entity_ = entity; weight_.setEntity( entity ); }
        void unbind () {}

        const EntityType entity() const { return entity_; }

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          const RangeFieldType weight = weight_( coordinate( x ) );
          for( int i = 0; i < dimRange; ++i )
            value[ i ] = weight;
        }

        template< class Quadrature, class Values >
        auto evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
          -> std::enable_if_t< std::is_same< decltype( values[ 0 ] ), RangeType & >::value >
        {
          for( const auto &qp : quadrature )
            evaluate( qp, values[ qp.index() ] );
        }

        int order() const
        { return order_; }
      private:
        Weight &weight_;
        int order_;
        Entity entity_;
      };

    } // namespace Impl



    // interpolate with weights
    // ------------------------

    template< class GridFunction, class DiscreteFunction, class Weight >
    inline static auto interpolate ( const GridFunction &u, DiscreteFunction &v, Weight &&weight )
      -> std::enable_if_t< !std::is_const< std::remove_reference_t< Weight > >::value  >
    {
      DiscreteFunction w( v );
      interpolate( u, v, std::forward< Weight >( weight ), w );
    }

    template< class GridFunction, class DiscreteFunction, class Weight >
    inline static auto interpolate ( const GridFunction &u, DiscreteFunction &v, Weight &&weight )
      -> std::enable_if_t< !std::is_const< std::remove_reference_t< Weight > >::value && !std::is_base_of< HasLocalFunction, GridFunction >::value >
    {
      interpolate( gridFunctionAdapter( u, v.gridPart() ), v, std::forward< Weight >( weight ) );
    }

    template< class GridFunction, class DiscreteFunction, class Weight >
    inline static auto interpolate ( const GridFunction &u, DiscreteFunction &v, Weight &&weight, DiscreteFunction &w )
      -> std::enable_if_t< std::is_base_of< HasLocalFunction, GridFunction >::value && Capabilities::hasInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >::v,
                           void_t< decltype( std::declval< Weight >().setEntity( std::declval< const typename DiscreteFunction::DiscreteFunctionSpaceType::EntityType & >() ) ) > >
    {
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType::EntityType EntityType;

      const auto &space = w.space();
      Impl::WeightLocalFunction< EntityType, std::remove_reference_t< typename DiscreteFunction::FunctionSpaceType >, Weight > localWeight( weight, w.order() );
      LocalInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType > interpolation( space );
      interpolate( u, v, [ &interpolation, &localWeight ] ( const EntityType &entity, AddLocalContribution< DiscreteFunction > &w ) {
          auto weightGuard = bindGuard( localWeight, entity );
          auto iGuard = bindGuard( interpolation, entity );
          interpolation( localWeight, w );
        }, w );
    }

    template< class GridFunction, class DiscreteFunction, class Weight >
    inline static auto interpolate ( const GridFunction &u, DiscreteFunction &v, Weight &&weight, DiscreteFunction &w )
      -> std::enable_if_t< std::is_base_of< HasLocalFunction, GridFunction >::value && Capabilities::hasInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >::v,
                           void_t< decltype( std::declval< Weight >()( std::declval< typename DiscreteFunction::DiscreteFunctionSpaceType::EntityType & >(), std::declval< AddLocalContribution< DiscreteFunction > & >() ) ) > >
    {
      typedef typename DiscreteFunction::DofType DofType;

      v.clear();
      w.clear();

      {
        ConstLocalFunction< GridFunction > uLocal( u );
        AddLocalContribution< DiscreteFunction > vLocal( v ), wLocal( w );
        LocalInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType >
          interpolation( v.space() );

        for( const auto &entity : v.space() )
        {
          auto uGuard = bindGuard( uLocal, entity );
          auto vGuard = bindGuard( vLocal, entity );
          auto wGuard = bindGuard( wLocal, entity );
          auto iGuard = bindGuard( interpolation, entity );

          // interpolate u and store in v
          interpolation( uLocal, vLocal );

          // evaluate DoF-wise weight
          weight( entity, wLocal );

          // multiply interpolated values by weight
          std::transform( vLocal.begin(), vLocal.end(), wLocal.begin(), vLocal.begin(), std::multiplies< DofType >() );
        }
      } // ensure the local contributions go out of scope, here (communication)

      std::transform( v.dbegin(), v.dend(), w.dbegin(), v.dbegin(), [] ( DofType v, DofType w ) {
          // weights are non-negative, so cancellation cannot occur
          return (w > DofType( 0 ) ? v / w : v);
        } );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
