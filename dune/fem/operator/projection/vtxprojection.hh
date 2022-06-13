#ifndef DUNE_FEM_VTXPROJECTION_HH
#define DUNE_FEM_VTXPROJECTION_HH

#include <cstddef>
#include <cmath>

#include <algorithm>
#include <functional>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune
{
  namespace Fem
  {

    template <class GridPartType>
    struct WeightDefault {
      typedef typename GridPartType::template Codim<0>::EntityType EntityType;
      typedef FieldVector<double,GridPartType::dimensionworld> WorldDomainType;
      typedef FieldVector<double,GridPartType::dimension> DomainType;
      void setEntity(const EntityType& en) {
        en_ = en;
        volume_ = en.geometry().volume();
        baryCenter_ = en.geometry().center();
      }
      double operator()(const DomainType& point) {
        // return volume_;
        auto tau = en_.geometry().global(point);
        tau -= baryCenter_;
        return tau.two_norm() / volume_;
      }
      double volume_;
      WorldDomainType baryCenter_;
      EntityType en_;
    };

    struct VtxProjectionImpl
    {
      template< class DiscreteFunction, class Intersection >
      struct OutsideLocalFunction
      {
        typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

        typedef typename LocalFunctionType::EntityType EntityType;
        typedef typename LocalFunctionType::FunctionSpaceType FunctionSpaceType;

        typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;

        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
        typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

        typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

        static const int dimDomain = FunctionSpaceType::dimDomain;
        static const int dimRange = FunctionSpaceType::dimRange;

        typedef typename Intersection::LocalGeometry GeometryType;

        explicit OutsideLocalFunction ( const DiscreteFunction &df ) : order_(df.order()), localFunction_( df ) {}

        void init ( const EntityType &outside, const GeometryType &geoIn, const GeometryType &geoOut )
        {
          localFunction_.init( outside );
          geoIn_ = &geoIn;
          geoOut_ = &geoOut;
          entity_ = outside;
        }

        const EntityType &entity() const { return entity_; }

        int order() const { return order_; }

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          localFunction_.evaluate( geoOut_->global( geoIn_->local( coordinate( x ) ) ), value );
        };

        template< class Point >
        void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
        {
          DUNE_THROW( NotImplemented, "Vertex projection weights do not provide Jacobians." );
        }

        template< class Point >
        void hessian ( const Point &x, HessianRangeType &hessian ) const
        {
          DUNE_THROW( NotImplemented, "Vertex projection weights do not provide Hessians." );
        }

        template< class Quadrature, class Values >
        void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
        {
          for( const auto &qp : quadrature )
            doEvaluate( qp, values[ qp.index() ] );
        }

      private:
        template< class Point >
        void doEvaluate ( const Point &x, RangeType &value ) const
        {
          evaluate( x, value );
        };

        template< class Point >
        void doEvaluate ( const Point &x, JacobianRangeType &value ) const
        {
          jacobian( x, value );
        };

        template< class Point >
        void doEvaluate ( const Point &x, HessianRangeType &value ) const
        {
          hessian( x, value );
        };

        int order_;
        LocalFunctionType localFunction_;
        const GeometryType *geoIn_ = nullptr;
        const GeometryType *geoOut_ = nullptr;
        EntityType entity_;
      };


      template< class DiscreteFunction >
      static void makeConforming ( DiscreteFunction &u )
      {
        typedef typename DiscreteFunction::GridPartType GridPartType;

        typedef typename GridPartType::IntersectionType IntersectionType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType
          DiscreteFunctionSpaceType;

        if( u.space().gridPart().isConforming() || !Fem::GridPartCapabilities::hasGrid< GridPartType >::v )
          return;

        const auto &blockMapper = u.space().blockMapper();

        std::vector< typename DiscreteFunction::DofType > ldu, ldw;
        const std::size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;
        const std::size_t maxNumBlocks = blockMapper.maxNumDofs();
        ldu.reserve( maxNumBlocks * localBlockSize );
        ldw.reserve( maxNumBlocks * localBlockSize );

        std::vector< bool > onSubEntity;
        onSubEntity.reserve( maxNumBlocks );

        OutsideLocalFunction< DiscreteFunction, IntersectionType > uOutside( u );

        LocalInterpolation< DiscreteFunctionSpaceType > interpolation( u.space() );

        for( const EntityType &inside : u.space() )
        {
          for( const IntersectionType &intersection : intersections( u.gridPart(), inside ) )
          {
            if( intersection.conforming() || !intersection.neighbor() )
              continue;

            // skip this intersection, if inside is finer than outside
            const EntityType outside = intersection.outside();
            if( inside.level() <= outside.level() )
              continue;

            // initialize "outside local function"
            // note: intersection geometries must live until end of intersection loop!
            const auto geoIn = intersection.geometryInInside();
            const auto geoOut = intersection.geometryInOutside();
            uOutside.init( outside, geoIn, geoOut );

            // resize local DoF vectors
            const std::size_t numBlocks = blockMapper.numDofs( inside );
            ldu.resize( numBlocks * localBlockSize );
            ldw.resize( numBlocks * localBlockSize );

            // interpolate "outside" values
            auto guard = bindGuard( interpolation, inside );
            interpolation( uOutside, ldw );

            // fetch inside local DoFs
            u.getLocalDofs( inside, ldu );

            // patch DoFs assigned to the face
            blockMapper.onSubEntity( inside, intersection.indexInInside(), 1, onSubEntity );
            for( std::size_t i = 0; i < numBlocks; ++i )
            {
              if( !onSubEntity[ i ] )
                continue;
              for( std::size_t j = 0; j < localBlockSize; ++j )
                ldu[ i*localBlockSize + j ] = ldw[ i*localBlockSize + j ];
            }

            // write back inside local DoFs
            u.setLocalDofs( inside, ldu );
          }
        }
      }

      template< class Function, class DiscreteFunction, class Weight >
      static auto project ( const Function &f, DiscreteFunction &u, Weight &weight )
      -> std::enable_if_t<std::is_void< decltype( interpolate(f,u,weight,u )) >::value>
      // -> std::enable_if_t<std::is_void< decltype( interpolate(f,u,weight,std::declval<std::remove_cv<std::remove_reference<DiscreteFunction>>&>()) ) >::value>
      {
        DiscreteFunction w ( "weight", u.space() );
        interpolate( f, u, weight, w );
        makeConforming( u );
      }
      template< class Function, class DiscreteFunction >
      static auto project ( const Function &f, DiscreteFunction &u )
      -> std::enable_if_t< std::is_void< decltype( project( f,u,std::declval<WeightDefault<typename DiscreteFunction::GridPartType>&>() ) ) >::value>
      {
        WeightDefault<typename DiscreteFunction::GridPartType> weight;
        project(f, u, weight);
      }
    };

    /*======================================================================*/
    /*! \ingroup VtxProjectionOperator
     *  \class VtxProjection
     *  \brief The Projection class which average
     *         discontinuous function using the interpolation of the space
     *         (e.g. the Lagrangepoints)
     */
    /*======================================================================*/
    template < typename DType, typename RType >
    class VtxProjection : public Operator< DType, RType >
    {
     public:
      typedef DType DomainType;
      typedef RType RangeType;
      typedef typename DomainType :: RangeFieldType DomainFieldType;
      typedef typename RType :: RangeFieldType  RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef typename RType :: DiscreteFunctionSpaceType :: GridPartType GridPartType;
      //! Constructor
      VtxProjection() {}

      //! apply projection
      template <class WeightType>
      void operator() (const DomainType& f, RangeType& discFunc,
                       WeightType& weight) const
      {
        VtxProjectionImpl::project(f,discFunc,weight);
      }
      //! apply projection with default weight
      void operator() (const DomainType& f, RangeType& discFunc) const
      {
        VtxProjectionImpl::project(f,discFunc);
      }

    private:
    };

  } // namespace Fem

} // name space Dune

#endif // #ifndef DUNE_FEM_VTXPROJECTION_HH
