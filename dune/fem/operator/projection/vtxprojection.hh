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
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/temporary.hh>

namespace Dune
{
  namespace Fem
  {

    template <class GridPartType>
    struct WeightDefault {
      typedef FieldVector<double,GridPartType::dimension> DomainType;
      template <class EntityType>
      void setEntity(const EntityType& en) {
        volume_ = en.geometry().volume();
      }
      double operator()(const DomainType& point) {
        return volume_;
      }
      double volume_;
    };

    struct VtxProjectionImpl
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

        static const int dimDomain = FunctionSpaceType::dimDomain;
        static const int dimRange = FunctionSpaceType::dimRange;

        explicit WeightLocalFunction ( Weight &weight ) : weight_( weight ) {}

        void bind ( const EntityType &entity ) { weight_.setEntity( entity ); }
        void unbind () {}

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          RangeFieldType weight = weight_( coordinate( x ) );
          for( int i = 0; i < dimRange; ++i )
            value[ i ] = weight;
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

        Weight &weight_;
      };

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

        explicit OutsideLocalFunction ( const DiscreteFunction &df ) : localFunction_( df ) {}

        void bind ( const EntityType &outside, const GeometryType &geoIn, const GeometryType &geoOut )
        {
          localFunction_.bind( outside );
          geoIn_ = &geoIn;
          geoOut_ = &geoOut;
        }
        void unbind() { localFunction_.unbind(); }

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

        ConstLocalFunction<DiscreteFunction> localFunction_;
        const GeometryType *geoIn_ = nullptr;
        const GeometryType *geoOut_ = nullptr;
      };

      template< class Function, class DiscreteFunction, class Weight >
      static void project ( const Function &f, DiscreteFunction &u, Weight &weight )
      {
        typedef typename Function::FunctionSpaceType FunctionSpaceType;

        typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename DiscreteFunction::GridPartType GridPartType;
        typedef typename DiscreteFunction::DofType DofType;

        typedef typename GridPartType::IntersectionType IntersectionType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        const DiscreteFunctionSpaceType &space = u.space();
        const auto &blockMapper = space.blockMapper();
        DiscreteFunction w ( "weight", space );

        u.clear();
        w.clear();

        {
          WeightLocalFunction< EntityType, FunctionSpaceType, Weight > localWeight( weight );
          Dune::Fem::ConstLocalFunction< Function > lf( f );
          Dune::Fem::AddLocalContribution< DiscreteFunction> lw( w );
          Dune::Fem::AddLocalContribution< DiscreteFunction> lu( u );

          for ( const auto& entity : space )
          {
            // initialize local function and local weight
            auto lfGuard = bindGuard( lf, entity );
            auto luGuard = bindGuard( lu, entity );
            auto lwGuard = bindGuard( lw, entity );
            auto localWeightGuard = bindGuard( localWeight, entity );

            // interpolate function values and weights
            const auto interpolation = space.interpolation( entity );
            interpolation( lf, lu );
            interpolation( localWeight, lw );

            // multiply function values by weights
            // Q: add begin/end directly to local contributions?
            std::transform( lu.localDofVector().begin(), lu.localDofVector().end(),
                            lw.localDofVector().begin(), lu.localDofVector().begin(),
                            std::multiplies< DofType >() );
          }

        } // make sure that everything is communicated

        // divide DoFs by according weights
        std::transform( u.dbegin(), u.dend(), w.dbegin(), u.dbegin(), [] ( DofType u, DofType w ) {
            using std::abs;
            typename Dune::FieldTraits< DofType >::field_type weight = abs( w );
            return (weight > 1e-12 ? u / weight : u);
          } );

        // make function continuous over hanging nodes
        if( !GridPartType::Traits::conforming && Fem::GridPartCapabilities::hasGrid< GridPartType >::v)
        {
          std::vector< bool > onSubEntity;
          onSubEntity.reserve( blockMapper.maxNumDofs() );
          std::vector< bool > mask;

          OutsideLocalFunction< DiscreteFunction, IntersectionType > uOutside( u );
          Dune::Fem::SetLocalContribution< DiscreteFunction> lu( u );
          std::vector<double> lutmp( space.maxNumDofs() );

          for( const auto &inside : space )
          {
            lu.bind(inside);
            mask.resize( lu.size() );
            std::fill( mask.begin(), mask.end(), false );
            const std::size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;
            const std::size_t numBlocks = blockMapper.numDofs( inside );

            for( const auto &intersection : intersections( space.gridPart(), inside ) )
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

              auto uOutsideGuard = bindGuard( uOutside, outside, geoIn, geoOut );

              // interpolate "outside" values
              const auto interpolation = space.interpolation( inside );
              interpolation( uOutside, lutmp );

              // patch DoFs assigned to the face
              blockMapper.onSubEntity( inside, intersection.indexInInside(), 1, onSubEntity );
              for( std::size_t i = 0; i < numBlocks; ++i )
              {
                if( !onSubEntity[ i ] )
                  continue;
                for( std::size_t j = 0; j < localBlockSize; ++j )
                {
                  lu[ i*localBlockSize + j ] = lutmp[ i*localBlockSize + j ];
                  mask[ i*localBlockSize + j ] = true;
                }
              }
            }
            lu.unbind(mask);
          }
        }
      }
    };

    /*======================================================================*/
    /*! \ingroup VtxProjectionOperator
     *  \class VtxProjection
     *  \brief The Projection class which average
     *         discontinuous function in the Lagrangepoints
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
        WeightDefault<GridPartType> weight;
        VtxProjectionImpl::project(f,discFunc, weight);
      }

    private:
    };

  } // namespace Fem

} // name space Dune

#endif // #ifndef DUNE_FEM_VTXPROJECTION_HH
