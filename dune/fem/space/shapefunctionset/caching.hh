#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH

// C++ includes
#include <cstddef>
#include <vector>
#include <type_traits>

#include <dune/geometry/quadraturerules.hh>

// dune-fem includes
#include <dune/fem/common/utility.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/quadrature/caching/registry.hh>
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{

  namespace Fem
  {

    // CachingShapeFunctionSet
    // -----------------------

    template< class ShapeFunctionSet >
    class CachingShapeFunctionSet
    : private QuadratureStorageRegistry::StorageInterface
    {
      typedef CachingShapeFunctionSet< ShapeFunctionSet > ThisType;

    public:
      typedef ShapeFunctionSet  ShapeFunctionSetType;

      typedef typename ShapeFunctionSet::FunctionSpaceType FunctionSpaceType;

      typedef typename ShapeFunctionSet::DomainType DomainType;
      typedef typename ShapeFunctionSet::RangeType RangeType;
      typedef typename ShapeFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename ShapeFunctionSet::HessianRangeType HessianRangeType;

      typedef std::vector< RangeType >          RangeVectorType ;
      typedef std::vector< JacobianRangeType >  JacobianRangeVectorType ;

      typedef std::vector< RangeVectorType >         RangeCacheVectorType;
      typedef std::vector< JacobianRangeVectorType > JacobianCacheVectorType;

    public:
      // point set id if available (otherwise -1)
      static const int pointSetId = detail::SelectPointSetId< ShapeFunctionSetType >::value;

      explicit CachingShapeFunctionSet ( const GeometryType &type,
                                         const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : type_( type )
      , shapeFunctionSet_( shapeFunctionSet )
      {
        QuadratureStorageRegistry::registerStorage( *this );
      }

      ~CachingShapeFunctionSet ();

      int order () const
      {
        return shapeFunctionSet_.order();
      }

      std::size_t size () const
      {
        return shapeFunctionSet_.size();
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.evaluateEach( x, functor );
      }

      template< class Quadrature, class Functor >
      void evaluateEach ( const QuadraturePointWrapper< Quadrature > &x, Functor functor ) const
      {
        const bool cacheable = std::is_convertible< Quadrature, CachingInterface >::value;
        evaluateEach( x.quadrature(), x.index(), functor, std::integral_constant< bool, cacheable >() );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.jacobianEach( x, functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const QuadraturePointWrapper< Quadrature > &x, Functor functor ) const
      {
        const bool cacheable = std::is_convertible< Quadrature, CachingInterface >::value;
        jacobianEach( x.quadrature(), x.index(), functor, std::integral_constant< bool, cacheable >() );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.hessianEach( x, functor );
      }

      GeometryType type () const {  return type_; }

      template < class QuadratureType >
      const RangeVectorType& rangeCache( const QuadratureType& quadrature ) const
      {
        return ReturnCache< QuadratureType, std::is_base_of< CachingInterface, QuadratureType >::value > ::
          ranges( *this, quadrature, rangeCaches_, localRangeCache_ );
      }

      template < class QuadratureType >
      const JacobianRangeVectorType& jacobianCache( const QuadratureType& quadrature ) const
      {
        return ReturnCache< QuadratureType, std::is_base_of< CachingInterface, QuadratureType >::value > ::
          jacobians( *this, quadrature, jacobianCaches_, localJacobianCache_ );
      }

      const ThisType& scalarShapeFunctionSet() const { return *this; }
      const ThisType& impl() const { return *this; }

    private:
      template< class Quad, bool cacheable /* false */ >
      struct ReturnCache
      {
        static const RangeVectorType&
        ranges( const ThisType& shapeFunctionSet,
                const Quad& quad,
                const RangeCacheVectorType&,
                RangeVectorType& storage )
        {
          // this feature was disabled in default_codegen.hh
          // this method should therefore not be called.
          assert( false );
          std::abort();

          // evaluate all basis functions and multiply with dof value
          const unsigned int nop  = quad.nop();
          const unsigned int size = shapeFunctionSet.size();

          // make sure cache has the appropriate size
          storage.resize( size * nop );
          RangeType* data = storage.data();

          for( unsigned int qp = 0 ; qp < nop; ++ qp )
          {
            const int cacheQp = quad.cachingPoint( qp );
            AssignFunctor< RangeType* > funztor( data + (cacheQp * size) );
            shapeFunctionSet.evaluateEach( quad[ qp ], funztor );
          }
          return storage;
        }

        static const JacobianRangeVectorType&
        jacobians( const ThisType& shapeFunctionSet,
                   const Quad& quad,
                   const JacobianCacheVectorType&,
                   JacobianRangeVectorType& storage )
        {
          // this feature was disabled in default_codegen.hh
          // this method should therefore not be called.
          assert( false );
          std::abort();

          // evaluate all basis functions and multiply with dof value
          const unsigned int nop  = quad.nop();
          const unsigned int size = shapeFunctionSet.size();

          // make sure cache has the appropriate size
          storage.resize( size * nop );
          JacobianRangeType* data = storage.data();

          for( unsigned int qp = 0 ; qp < nop; ++ qp )
          {
            const int cacheQp = quad.cachingPoint( qp );
            AssignFunctor< JacobianRangeType* > funztor( data + ( cacheQp * size ) );
            shapeFunctionSet.jacobianEach( quad[ qp ], funztor );
          }
          return storage;
        }
      };

      template< class Quad >
      struct ReturnCache< Quad, true >
      {
        static const RangeVectorType&
        ranges( const ThisType& shapeFunctionSet,
                const Quad& quad,
                const RangeCacheVectorType& cache,
                const RangeVectorType& )
        {
          return cache[ quad.id() ];
        }

        static const JacobianRangeVectorType&
        jacobians( const ThisType& shapeFunctionSet,
                   const Quad& quad,
                   const JacobianCacheVectorType& cache,
                   const JacobianRangeVectorType& )
        {
          return cache[ quad.id() ];
        }
      };


      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          std::integral_constant< bool, false > ) const
      {
        evaluateEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          std::integral_constant< bool, true > ) const;

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          std::integral_constant< bool, false > ) const
      {
        jacobianEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          std::integral_constant< bool, true > ) const;


      void cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size );

      template< class PointVector >
      void cachePoints ( std::size_t id, const PointVector &points );

      GeometryType type_;
      ShapeFunctionSet shapeFunctionSet_;
      RangeCacheVectorType rangeCaches_;
      JacobianCacheVectorType jacobianCaches_;

      mutable RangeVectorType          localRangeCache_ ;
      mutable JacobianRangeVectorType  localJacobianCache_;
    };



    // Implementation of CachingShapeFunctionSet
    // -----------------------------------------

    template< class ShapeFunctionSet >
    inline CachingShapeFunctionSet< ShapeFunctionSet >::~CachingShapeFunctionSet ()
    {
      QuadratureStorageRegistry::unregisterStorage( *this );
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       std::integral_constant< bool, true > ) const
    {
      assert( (quadrature.id() < rangeCaches_.size()) && !rangeCaches_[ quadrature.id() ].empty() );
      const RangeType *cache = rangeCaches_[ quadrature.id() ].data();

      const unsigned int numShapeFunctions = size();
      const unsigned int cpt = quadrature.cachingPoint( pt );

      // for Lagrange-type basis evaluated on interpolation points
      // this is the Kronecker delta, there we only need
      // to evaluate the shapefunction with number 'pt'
      static const int quadPointSetId = SelectQuadraturePointSetId< Quadrature >::value;
      //std::cout << "QP:" << quadPointSetId << " " << pointSetId << std::endl;

      if constexpr ( quadPointSetId == pointSetId )
      {
        // std::cout << "QP matches: " << quadPointSetId << " " << pointSetId << " " << quadrature.nop() << " " << numShapeFunctions << std::endl;
        if( quadrature.isInterpolationQuadrature(numShapeFunctions) )
        {
          // negative values mean invalid point sets
          // we should not get here in this case
          assert( quadPointSetId >= 0 );
          assert( pointSetId >= 0 );

          //if( Quadrature::codimension == 1 )
          //  std::cout << "using interpolation point " << std::endl;
          // obtain interpolation point (different for face quadratures)
          const unsigned int point = quadrature.interpolationPoint( pt );
#ifndef NDEBUG
          for( unsigned int i = 0; i < numShapeFunctions; ++i )
            assert( (cache[ cpt*numShapeFunctions + i ] - RangeType(i==point)).two_norm() < 1e-8 ) ;
#endif
          // point should be 1
          functor( point, RangeType(1) );
          return;
        }
      }

      for( unsigned int i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ cpt*numShapeFunctions + i ] );
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       std::integral_constant< bool, true > ) const
    {
      assert( (quadrature.id() < jacobianCaches_.size()) && !jacobianCaches_[ quadrature.id() ].empty() );
      const JacobianRangeType *cache = jacobianCaches_[ quadrature.id() ].data();

      const unsigned int numShapeFunctions = size();
      const unsigned int cpt = quadrature.cachingPoint( pt );
      for( unsigned int i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ cpt*numShapeFunctions + i ] );
    }


    template< class ShapeFunctionSet >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size )
    {
      if( ! MPIManager::singleThreadMode () )
      {
        DUNE_THROW(SingleThreadModeError,"CachingShapeFunctionSet::cacheQuadrature: only call in single thread mode!");
      }

      if( id >= rangeCaches_.size() )
      {
        rangeCaches_.resize( id+1, RangeVectorType() );
        jacobianCaches_.resize( id+1, JacobianRangeVectorType() );
      }

      assert( rangeCaches_[ id ].empty() == jacobianCaches_[ id ].empty() );

      if( rangeCaches_[ id ].empty() )
      {
        typedef typename FunctionSpaceType::DomainFieldType ctype;
        const int dim = FunctionSpaceType::dimDomain;
        switch( codim )
        {
        case 0:
          cachePoints( id, PointProvider< ctype, dim, 0 >::getPoints( id, type_ ) );
          break;

        case 1:
          cachePoints( id, PointProvider< ctype, dim, 1 >::getPoints( id, type_ ) );
          break;

        default:
          DUNE_THROW( NotImplemented, "Caching for codim > 1 not implemented." );
        }
      }
    }


    template< class ShapeFunctionSet >
    template< class PointVector >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::cachePoints ( std::size_t id, const PointVector &points )
    {
      const unsigned int numShapeFunctions = size();
      const unsigned int numPoints = points.size();

      RangeVectorType& ranges = rangeCaches_[ id ];
      ranges.resize( numShapeFunctions * numPoints );

      JacobianRangeVectorType& jacobians = jacobianCaches_[ id ];
      jacobians.resize( numShapeFunctions * numPoints );

      if( ranges.empty() || jacobians.empty() )
        DUNE_THROW( OutOfMemoryError, "Unable to allocate shape function set caches." );

      for( unsigned int pt = 0; pt < numPoints; ++pt )
      {
        evaluateEach( points[ pt ], AssignFunctor< RangeType * >( ranges.data() + pt*numShapeFunctions ) );
        jacobianEach( points[ pt ], AssignFunctor< JacobianRangeType * >( jacobians.data() + pt*numShapeFunctions ) );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
