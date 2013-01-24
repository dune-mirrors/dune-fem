#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH

// C++ includes
#include <cstddef>
#include <vector>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

// dune-fem includes
#include <dune/fem/misc/functor.hh>
#include <dune/fem/quadrature/caching/registry.hh>
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/quadrature/quadrature.hh>

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
      typedef typename ShapeFunctionSet::FunctionSpaceType FunctionSpaceType;
      
      typedef typename ShapeFunctionSet::DomainType DomainType;
      typedef typename ShapeFunctionSet::RangeType RangeType;
      typedef typename ShapeFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename ShapeFunctionSet::HessianRangeType HessianRangeType;

      explicit CachingShapeFunctionSet ( const GeometryType &type,
                                         const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : type_( type ),
        shapeFunctionSet_( shapeFunctionSet )
      {
        QuadratureStorageRegistry::registerStorage( *this );
      }

      ~CachingShapeFunctionSet ();

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
        const bool cacheable = Conversion< Quadrature, CachingInterface >::exists;
        evaluateEach( x.quadrature(), x.point(), functor, integral_constant< bool, cacheable >() );
      }
     
      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.jacobianEach( x, functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const QuadraturePointWrapper< Quadrature > &x, Functor functor ) const
      {
        const bool cacheable = Conversion< Quadrature, CachingInterface >::exists;
        jacobianEach( x.quadrature(), x.point(), functor, integral_constant< bool, cacheable >() );
      }

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.hessianEach( x, functor );
      }

      GeometryType type () const DUNE_DEPRECATED {  return type_; }
      
      GeometryType geometryType () const DUNE_DEPRECATED {  return type(); }

    private:
      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, false > ) const
      {
        evaluateEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, true > ) const;

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, false > ) const
      {
        jacobianEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, true > ) const;


      void cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size );

      template< class PointVector >
      void cachePoints ( std::size_t id, const PointVector &points );

      // prohibit copying and assignment
      CachingShapeFunctionSet ( const ThisType & );
      const ThisType &operator= ( const ThisType & );

      typedef std::vector< RangeType * > ValueCacheVectorType;
      typedef std::vector< JacobianRangeType * > JacobianCacheVectorType;

      GeometryType type_;
      ShapeFunctionSet shapeFunctionSet_;
      ValueCacheVectorType valueCaches_;
      JacobianCacheVectorType jacobianCaches_;
    };



    // Implementation of CachingShapeFunctionSet
    // -----------------------------------------

    template< class ShapeFunctionSet >
    inline CachingShapeFunctionSet< ShapeFunctionSet >::~CachingShapeFunctionSet ()
    {
      QuadratureStorageRegistry::unregisterStorage( *this );
      for( typename ValueCacheVectorType::iterator it = valueCaches_.begin(); it != valueCaches_.end(); ++it )
        delete *it;
      for( typename JacobianCacheVectorType::iterator it = jacobianCaches_.begin(); it != jacobianCaches_.end(); ++it )
        delete *it;
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       integral_constant< bool, true > ) const
    {
      assert( (quadrature.id() < valueCaches_.size()) && valueCaches_[ quadrature.id() ] );
      const RangeType *cache = valueCaches_[ quadrature.id() ];

      const std::size_t numShapeFunctions = size();
      const std::size_t cpt = quadrature.cachingPoint( pt );
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ cpt*numShapeFunctions + i ] );
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       integral_constant< bool, true > ) const
    {
      assert( (quadrature.id() < jacobianCaches_.size()) && jacobianCaches_[ quadrature.id() ] );
      const JacobianRangeType *cache = jacobianCaches_[ quadrature.id() ];

      const std::size_t numShapeFunctions = size();
      const std::size_t cpt = quadrature.cachingPoint( pt );
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ cpt*numShapeFunctions + i ] );
    }


    template< class ShapeFunctionSet >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size )
    {
      if( id >= valueCaches_.size() )
      {
        valueCaches_.resize( id+1, nullptr );
        jacobianCaches_.resize( id+1, nullptr );
      }
      assert( bool( valueCaches_[ id ] ) == bool( jacobianCaches_[ id ] ) );

      if( !valueCaches_[ id ] )
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
      const std::size_t numShapeFunctions = size();
      const std::size_t numPoints = points.size();
      RangeType *values = new RangeType[ numShapeFunctions * numPoints ];
      JacobianRangeType *jacobians = new JacobianRangeType[ numShapeFunctions * numPoints ];
      if( !values || !jacobians )
        DUNE_THROW( OutOfMemoryError, "Unable to allocate shape function set caches." );

      for( std::size_t pt = 0; pt < numPoints; ++pt )
      {
        evaluateEach( points[ pt ], AssignFunctor< RangeType * >( values + pt*numShapeFunctions ) );
        jacobianEach( points[ pt ], AssignFunctor< JacobianRangeType * >( jacobians + pt*numShapeFunctions ) );
      }

      valueCaches_[ id ] = values;
      jacobianCaches_[ id ] = jacobians;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
